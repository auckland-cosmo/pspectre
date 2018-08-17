/*
 * SpectRE - A Spectral Code for Reheating
 * Copyright (C) 2009-2010 Hal Finkel, Nathaniel Roth and Richard Easther
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
 * AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#define _XOPEN_SOURCE 600

#include "pow/pow.hpp"
#include "model.hpp"
#include "integrator.hpp"
#include "verlet.hpp"
#include "rk4.hpp"
#include "initializer.hpp"
#include "le_style_initializer.hpp"
#include "defrost_style_initializer.hpp"
#include "grid_funcs.hpp"
#include "spectra_outputter.hpp"
#include "twoptcorr_outputter.hpp"
#include "stats_outputter.hpp"
#include "energy_outputter.hpp"

#include <cstdlib>
#include <cstring>

#include <ctime>

#include <algorithm>

#include <iostream>
#include <fstream>
#include <iomanip>

#include <unistd.h>

#include <sys/stat.h>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// Forward declaration of integrator classes...
template <typename R> class integrator;
template <typename R> class verlet;
template <typename R> class rk4;

template <typename R>
const char *precision_name()
{
	return "unknown";
}

template <>
const char *precision_name<double>()
{
	return "double";
}

template <>
const char *precision_name<long double>()
{
	return "extended";
}

template <typename R>
void model<R>::set_output_directory(const char *uodn)
{
	time_t curt = time(NULL);
	tm *curtm = localtime(&curt);

	char outdirname[2048];
	strftime(outdirname, 2048, "output-%Y%m%d%H%M%S", curtm);
	if (!uodn || uodn[0] == '\0') {
		uodn = outdirname;
	}

	mkdir(uodn, 0744);
	chdir(uodn);	
}

/**
 * @page running Running
 *
 * @section Command-Line Parameters
 * SpectRE Usage:
 * @code
 * ./pspectre [-h]
 * ./pspectre [-r] [-l [-B <real>]] [-V] [-H <name>[,<name>]*] [-O] [-N <int>] [-P <int>] [-L <real>] [-R <int>] [-o <dir name>] [-t <real>[:<real>]] [-T <real>] [-A <real>] [-p <name>=<value>[,<name>=<value>]*] [-s <name>[,<name>]*] [-S <name>[=<value>][,<name>[=<value>]]*] [-I <name>=<value>[,<name>=<value>]*] [--long] [@<file name>]
 * @endcode
 *
 * @li -h: Display usage information and exit
 * @li -r: Use the RK4 integrator (default is the Verlet integrator)
 * @li -l: Use LatticeEasy-style initial conditions (default is DEFROST-style initial conditions)
 * @li -B: The base length scale (default is 1.0 to match LatticeEasy)
 * @li -V: Allow the field variance to change with L
 * @li -H \<name\>[,\<name\>]*: Use homogeneous (zero variance) initial conditions. Field names are:
 * @code
 *      phi
 *      chi
 * @endcode
 * @li -O: Use out-of-place transforms
 * @li -N \<int\>: The number of grid points per side of the box
 * @li -P \<int\>: The padding factor used for position-space integration
 * @li -L \<real\>: The physical size of the box
 * @li -R \<int\>: The random seed
 * @li -o \<dir name\>: Set the output directory name
 * @li -t \<real\>[:\<real\>]: Set dt with an optional start time in program units
 * @li -T \<real\>: The final time in program units
 * @li -A \<real\>: The final scale factor
 * @li -p \<name\>=\<value\>[,\<name\>=\<value\>]*: Set a parameter value. Valid parameters are:
 * @code
 *	gamma_phi
 *	gamma_chi
 *	lambda_phi
 *	lambda_chi
 *	g
 *	m_phi
 *	m_chi
 *	phi0
 *	chi0
 *	phidot0
 *	chidot0
 *	ics_scale
 * @endcode
 * @li -s \<name\>[,\<name\>]*: Enable slice output of a variable. Valid variables are:
 * @code
 *	phi
 *	chi
 *	V
 *	T_phi
 *	T_chi
 *	G_phi
 *	G_chi
 *	rho
 *	p
 *	gpot
 * @endcode
 * @li -S \<name\>[=\<value\>][,\<name\>[=\<value\>]]*: Set a slice output option value. Valid options are:
 * @code
 *	dim
 *	length
 *	skip
 *	avg
 *	(avg does not take a value)
 * @endcode
 * @li -I \<name\>=\<value\>[,\<name\>=\<value\>]*: Set an output interval. Valid intervals are:
 * @code
 *	scale
 *	energy
 *	spectra
 *	twoptcorr
 *	screen
 *	slice
 *	stats
 *	all
 *	(intervals are specified as a number of iterations)
 * @endcode
 * @li --long: Run using long-double (extended) precision (this must be the *last* command-line option argument)
 * @li \@\<file name\>: The name of a parameters file. The parameters file has the same syntax as the command
 * line except that it may be divided among multiple lines and may contain comment lines which begin with
 * a \# character.
 *
 * @par
 * The default parameters model a situation generally similar to the default model provided with DEFROST version 1.1.
 *
 * @section rexamples Examples
 *
 * The following runs the model with the default parameters except that it sets a 128^3 grid with dt = 0.0005. Also,
 * -r selects the RK4 integrator (Verlet is default). -l selects LE-style initial conditions. -I all=1 
 * sets all output intervals to 1 time step (the default is 25).
 *
 * @code
 * ./pspectre -N 128 -t 0.0005 -r -l -I all=1
 * @endcode
 *
 * The following runs the model with the default parameters and has binary slice outputs for the energy density, pressure
 * and gravitational potential. The slices to have a length of 32 points per side and were constructed by averaging (not skipping)
 * over every eight-point cube (since the dimension is 3). -P 2 causes the integration over the potential energy to use a (2N)^3 grid.
 *
 * @code
 * ./pspectre -P 2 -s rho,p,gpot -S dim=3,length=32,skip=1,avg
 * @endcode
 */

template <typename R>
model<R>::model(int argc, char *argv[])
	: fs(64), use_verlet(true), le_init(false), homo_ic_phi(false), homo_ic_chi(false), seed(1), tf(200.0),
	scale_interval(25), energy_interval(25), spectra_interval(25),
		screen_interval(25), slice_interval(25), stats_interval(25), twoptcorr_interval(25),
		phi("phi"), phidot("phidot"), chi("chi"), chidot("chidot"), gc(0), som(0), ics_scale(1), len0(1.0), vvwl(false), af(0.0)
#ifdef HAVE_PRIVATE
		, priv(phi, chi)
#endif
{
	char *subopts, *value;
	int opt;

	extern char *optarg;
	extern int optopt;
		
	const char *param_names[] = {
		"gamma_phi", "gamma_chi",
		"lambda_phi", "lambda_chi",
		"g", "m_phi", "m_chi",
		"phi0", "chi0", "phidot0", "chidot0",
		"ics_scale", 0
	};

	const char *interval_names[] = {
		"scale", "energy", "spectra",
		"screen", "slice", "stats",
		"twoptcorr", "all", 0
	};

	const char *slice_opt_names[] = {
		"dim", "length", "skip",
		"avg", 0
	};
	
	const char *slice_names[] = {
		"phi", "chi", "V", "T_phi", "T_chi",
		"G_phi", "G_chi", "rho", "p", "gpot", 0
	};

	const char *field_names[] = {
		"phi", "chi", 0
	};
	
	bool show_usage = false, help_requested = false;

	bool slice_phi = false, slice_chi = false, slice_V = false,
		slice_T_phi = false, slice_T_chi = false,
		slice_G_phi = false, slice_G_chi = false,
		slice_rho = false, slice_p = false, slice_gpot = false;

	int slicedim = 3, slicelength = 0, sliceskip = 1;
	bool sliceaverage = false;

	bool oop = false;
	string odn;

	while ((opt = getopt(argc, argv, ":rlVB:hH:ON:P:L:R:p:o:t:T:A:s:S:I:z:")) != -1) {
		switch (opt) {
		case 'h':
			help_requested = true;
			show_usage = true;
			break;
		case 'r':
			use_verlet = false;
			break;
		case 'l':
			le_init = true;
			break;
		case 'B':
			len0 = atof(optarg);
			break;
		case 'V':
			vvwl = true;
			break;
		case 'H':
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) field_names, &value);
				if (index < 0) {
					cerr << "Invalid field specification: " << value << endl;
					show_usage = true;
				}			
				else if (!strcmp(field_names[index], "phi")) {
					homo_ic_phi = true;
				}
				else if (!strcmp(field_names[index], "chi")) {
					homo_ic_chi = true;
				}
			}
			break;
		case 'O':
			oop = true;
			break;
		case 'N':
			fs.n = atoi(optarg);
			break;
		case 'P':
			fs.n_pad_factor = atoi(optarg);
			break;
		case 'L':
			mp.len = atof(optarg);
			break;
		case 'R':
			seed = atoi(optarg);
			break;
		case 'o':
			odn = optarg;
			break;
		case 'p':
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) param_names, &value);
				if (index < 0) {
					cerr << "Invalid parameter specification: " << value << endl;
					show_usage = true;
				}
				else if (!value) {
					cerr << "No value specified for parameter: " << param_names[index] << endl;
					show_usage = true;
				}
				else if (!strcmp(param_names[index], "gamma_phi")) {
					mp.gamma_phi = atof(value);
				}
				else if (!strcmp(param_names[index], "gamma_chi")) {
					mp.gamma_chi = atof(value);
				}
				else if (!strcmp(param_names[index], "lambda_phi")) {
					mp.lambda_phi = atof(value);
				}
				else if (!strcmp(param_names[index], "lambda_chi")) {
					mp.lambda_chi = atof(value);
				}
				else if (!strcmp(param_names[index], "g")) {
					mp.g = atof(value);
				}
				else if (!strcmp(param_names[index], "m_phi")) {
					mp.m_phi = atof(value);
				}
				else if (!strcmp(param_names[index], "m_chi")) {
					mp.m_chi = atof(value);
				}
				else if (!strcmp(param_names[index], "phi0")) {
					mp.phi0 = atof(value);
				}
				else if (!strcmp(param_names[index], "chi0")) {
					mp.chi0 = atof(value);
				}
				else if (!strcmp(param_names[index], "phidot0")) {
					mp.phidot0 = atof(value);
				}
				else if (!strcmp(param_names[index], "chidot0")) {
					mp.chidot0 = atof(value);
				}
				else if (!strcmp(param_names[index], "ics_scale")) {
					ics_scale = atof(value);
				}
			}
			break;
		case 't':
			// Parse as dt or dt:start_time
			{
				char *en;
				R dt = strtod(optarg, &en);
				R start_time = 0.0;
				if (*en != '0') {
					start_time = atof(en + 1);
				}

				ts.add_dt(start_time, dt);
			}
			break;
		case 'T':
			tf = atof(optarg);
			break;
		case 'A':
			af = atof(optarg);
			break;
		case 's':
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) slice_names, &value);
				if (index < 0) {
					cerr << "Invalid slice output specification: " << value << endl;
					show_usage = true;
				}			
				else if (!strcmp(slice_names[index], "phi")) {
					slice_phi = true;
				}
				else if (!strcmp(slice_names[index], "chi")) {
					slice_chi = true;
				}
				else if (!strcmp(slice_names[index], "V")) {
					slice_V = true;
				}
				else if (!strcmp(slice_names[index], "T_phi")) {
					slice_T_phi = true;
				}
				else if (!strcmp(slice_names[index], "T_chi")) {
					slice_T_chi = true;
				}
				else if (!strcmp(slice_names[index], "G_phi")) {
					slice_G_phi = true;
				}
				else if (!strcmp(slice_names[index], "G_chi")) {
					slice_G_chi = true;
				}
				else if (!strcmp(slice_names[index], "rho")) {
					slice_rho = true;
				}
				else if (!strcmp(slice_names[index], "p")) {
					slice_p = true;
				}
				else if (!strcmp(slice_names[index], "gpot")) {
					slice_gpot = true;
				}
			}
			break;
		case 'S':
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) slice_opt_names, &value);
				if (index < 0) {
					cerr << "Invalid slice option specification: " << value << endl;
					show_usage = true;
				}
				else if (!strcmp(slice_opt_names[index], "avg")) {
					sliceaverage = true;
				}
				else if (!value) {
					cerr << "No value specified for slice option: " << slice_opt_names[index] << endl;
					show_usage = true;
				}
				else if (!strcmp(slice_opt_names[index], "dim")) {
					slicedim = atoi(value);
				}
				else if (!strcmp(slice_opt_names[index], "length")) {
					slicelength = atoi(value);
				}
				else if (!strcmp(slice_opt_names[index], "skip")) {
					sliceskip = atoi(value) + 1;
				}
			}
			break;
		case 'I':
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) interval_names, &value);
				if (index < 0) {
					cerr << "Invalid interval specification: " << value << endl;
					show_usage = true;
				}
				else if (!value) {
					cerr << "No value specified for interval: " << param_names[index] << endl;
					show_usage = true;
				}
				else if (!strcmp(interval_names[index], "scale")) {
					scale_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "energy")) {
					energy_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "spectra")) {
					spectra_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "screen")) {
					screen_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "slice")) {
					slice_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "stats")) {
					stats_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "twoptcorr")) {
					twoptcorr_interval = atoi(value);
				}
				else if (!strcmp(interval_names[index], "all")) {
					scale_interval = energy_interval = spectra_interval =
						screen_interval = slice_interval = stats_interval =
						twoptcorr_interval = atoi(value);
				}
			}
			break;
		case 'z':
#ifdef HAVE_PRIVATE
			subopts = optarg;
			while (*subopts != '\0') {
				int index = getsubopt(&subopts, (char**) priv.get_opt_names(), &value);
				if (index < 0) {
					cerr << "Invalid parameter specification: " << value << endl;
					show_usage = true;
				}
				else if (!priv.process_opt(priv.get_opt_names()[index], value)) {
					cerr << "No/Invalid value specified for parameter: " << priv.get_opt_names()[index] << endl;
					show_usage = true;
				}
			}
#endif
			break;
		case ':':
			cerr << "Missing operand for option " << (char) optopt << endl;
			show_usage = true;
			break;
		case '?':
			cerr << "Unrecognized option: " << (char) optopt << endl;
			show_usage = true;
			break;
		}
	}

	if (optind < argc) {
		cerr << "The following options could not be parsed: \"";
		for (int i = optind; i < argc; ++i) {
			cerr << argv[i] << (i < (argc-1) ? " " : "");
		}
		cerr << "\"" << endl;
		show_usage = true;
	}

	if (show_usage) {
		ostream &hout = help_requested ? cout : cerr;
		
		hout << "SpectRE Usage:" << endl;
		hout << argv[0] << " [-h]" << endl;
		hout << argv[0] << " [-r] [-l [-B <real>]] [-V] [-H <name>[,<name>]*] [-O] [-N <int>] [-P <int>] [-L <real>] [-R <int>] "
			"[-o <dir name>] [-t <real>[:<real>]] [-T <real>] [-A <real>] "
			"[-p <name>=<value>[,<name>=<value>]*] [-s <name>[,<name>]*] [-S <name>[=<value>][,<name>[=<value>]]*] "
			"[-I <name>=<value>[,<name>=<value>]*] "
#ifdef HAVE_PRIVATE
			"[-z <name>[=<value>][,<name>[=<value>]]*] "
#endif
#ifdef USE_LD
			"[--long] "
#endif
			"[@<file name>]"
			<< endl;
		hout << endl;
		
		hout << "\t-h: Display usage information and exit" << endl;
		hout << "\t-r: Use the RK4 integrator (default is the Verlet integrator)" << endl;
		hout << "\t-l: Use LatticeEasy-style initial conditions (default is DEFROST-style initial conditions)" << endl;
		hout << "\t-B: The base length scale (default is 1.0 to match LatticeEasy)" << endl;
		hout << "\t-V: Allow the field variance to change with L" << endl;
		hout << "\t-H: Use homogeneous (zero variance) initial conditions. Field names are:" << endl;

		for (int i = 0; i < (int) (sizeof(field_names)/sizeof(field_names[0])) - 1; ++i) {
			hout << "\t\t" << field_names[i] << endl;
		}

		hout << "\t-O: Use out-of-place transforms" << endl;
		hout << "\t-N <int>: The number of grid points per side of the box" << endl;
		hout << "\t-P <int>: The padding factor used for position-space integration" << endl;
		hout << "\t-L <real>: The physical size of the box" << endl;
		hout << "\t-R <int>: The random seed" << endl;
		hout << "\t-o <dir name>: Set the output directory name" << endl;
		hout << "\t-t <real>[:<real>]: Set dt with an optional start time in program units" << endl;
		hout << "\t-T <real>: The final time in program units" << endl;
		hout << "\t-A <real>: The final scale factor" << endl;
		hout << "\t-p <name>=<value>[,<name>=<value>]*: Set a parameter value. Valid parameters are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(param_names)/sizeof(param_names[0])) - 1; ++i) {
			hout << "\t\t" << param_names[i] << endl;
		}

		hout << "\t-s <name>[,<name>]*: Enable slice output of a variable. Valid variables are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(slice_names)/sizeof(slice_names[0])) - 1; ++i) {
			hout << "\t\t" << slice_names[i] << endl;
		}

		hout << "\t-S <name>[=<value>][,<name>[=<value>]]*: Set a slice output option value. Valid options are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(slice_opt_names)/sizeof(slice_opt_names[0])) - 1; ++i) {
			hout << "\t\t" << slice_opt_names[i] << endl;
		}
		hout << "\t\t(avg does not take a value)" << endl;

		hout << "\t-I <name>=<value>[,<name>=<value>]*: Set an output interval. Valid intervals are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(interval_names)/sizeof(interval_names[0])) - 1; ++i) {
			hout << "\t\t" << interval_names[i] << endl;
		}
		hout << "\t\t(intervals are specified as a number of iterations)" << endl;

#ifdef HAVE_PRIVATE
		hout << "\t-z <name>[=<value>][,<name>[=<value>]]*: Set a private option value. Valid options are:" << endl;
		
		for (int i = 0; priv.get_opt_names()[i]; ++i) {
			hout << "\t\t" << priv.get_opt_names()[i] << endl;
		}
		hout << "\t\t(some options may not take a value)" << endl;
#endif

#ifdef USE_LD
		hout << "\t--long: Run using long-double (extended) precision (this must be the *last* command-line option argument)" << endl;
#endif

		hout << "\t@<file name>: The name of a parameters file" << endl;

		exit(help_requested ? 0 : 1);
	}

#ifdef MKL_NO_DCT
	if (!le_init) {
		cout << "Using LatticeEasy-style initial conditions (cannot use DEFROST-style initial conditions with MKL 10)" << endl;
		le_init = true;
	}
#endif

	ts.finalize_dts();

	srand48(seed);
	
	cout << "Using " << precision_name<R>() << " precision." << endl;

	mp.calculate_derived_params(true);
	mp.phidot0 /= mp.rescale_B;
	mp.chidot0 /= mp.rescale_B;

	fs.calculate_size_totals();

	phi.construct(fs, oop);
	phidot.construct(fs, oop);
	chi.construct(fs, oop);
	chidot.construct(fs, oop);

	set_output_directory(odn.c_str());

	gc = new grad_computer<R>(fs, mp, phi, chi);
	gpotc = new gpot_computer<R>(fs, mp, ts, phi, chi, phidot, chidot, *gc);
	som = new slice_output_manager<R>(fs, mp, ts, phi, chi, phidot, chidot, *gc, *gpotc,
		slicedim, slicelength, sliceskip, sliceaverage);
	
	if (slice_phi) som->add_outputter("phi", grid_funcs<R>::compute_phi);
	if (slice_chi) som->add_outputter("chi", grid_funcs<R>::compute_chi);
	if (slice_V) som->add_outputter("V", grid_funcs<R>::compute_V);
	if (slice_T_phi) som->add_outputter("T_phi", grid_funcs<R>::compute_T_phi);
	if (slice_T_chi) som->add_outputter("T_chi", grid_funcs<R>::compute_T_chi);
	if (slice_G_phi) som->add_outputter("G_phi", grid_funcs<R>::compute_G_phi);
	if (slice_G_chi) som->add_outputter("G_chi", grid_funcs<R>::compute_G_chi);
	if (slice_rho) som->add_outputter("rho", grid_funcs<R>::compute_rho);
	if (slice_p) som->add_outputter("p", grid_funcs<R>::compute_p);
	if (slice_gpot) som->add_outputter("gpot", grid_funcs<R>::compute_gpot);

#ifdef HAVE_PRIVATE
	priv.initialize(fs);
#endif
}

template <typename R>
model<R>::~model()
{
	delete gc;
	delete gpotc;
	delete som;
}

template <typename R>
void model<R>::set_initial_conditions()
{
	// This is the initial value of adot, and since a = 1 at t = 0, this is the initial value of H.
	// See equation 6.91 in the LatticeEasy manual.

	R hf = 3. * pow<2>(mp.rescale_A)/(4. * M_PI);
	R h0norm = 1. / (hf - pow<2>(mp.rescale_r*mp.rescale_A) * (pow<2>(mp.phidot0) + pow<2>(mp.chidot0)));
	ts.adot = h0norm * (
		-mp.rescale_r*pow<2>(mp.rescale_A)*(mp.phidot0*mp.phi0 + mp.chidot0*mp.chi0) +
		sqrt(
			hf * pow<2>(mp.rescale_A)*(pow<2>(mp.phidot0) + pow<2>(mp.chidot0)) +
			2. * mp.V(mp.rescale_A * mp.phi0, mp.rescale_A * mp.chi0, 1.) *
				(hf - pow<2>(mp.rescale_r*mp.rescale_A)*(pow<2>(mp.phi0) + pow<2>(mp.chi0)))
		)
	);

	if (!homo_ic_phi || !homo_ic_chi) {
		initializer<R> *init = le_init ?
			(initializer<R> *) new le_style_initializer<R>(fs, mp, 
				phi, phidot, chi, chidot, ts.adot, len0)
			: (initializer<R> *) new defrost_style_initializer<R>(fs, mp, 
				phi, phidot, chi, chidot, ts.adot);

		init->initialize();

		if (vvwl) {
			const R cf = pow(mp.len/(le_init ? len0 : R(1.0)), R(1.5))/pow<3>(2 * M_PI/(mp.len));

			phi.switch_state(momentum);
			phidot.switch_state(momentum);

			phi.divby(cf);
			phidot.divby(cf);

			chi.switch_state(momentum);
			chidot.switch_state(momentum);

			chi.divby(cf);
			chidot.divby(cf);
		}

		if (homo_ic_phi) {
			phi.switch_state(momentum);
			phidot.switch_state(momentum);
			fill((R *) phi.mdata, ((R *) phi.mdata) + 2*fs.total_momentum_gridpoints, 0);
			fill((R *) phidot.mdata, ((R *) phidot.mdata) + 2*fs.total_momentum_gridpoints, 0);
		}

		if (homo_ic_chi) {
			chi.switch_state(momentum);
			chidot.switch_state(momentum);
			fill((R *) chi.mdata, ((R *) chi.mdata) + 2*fs.total_momentum_gridpoints, 0);
			fill((R *) chidot.mdata, ((R *) chidot.mdata) + 2*fs.total_momentum_gridpoints, 0);
		}
	}
	else {
		phi.switch_state(momentum);
		chi.switch_state(momentum);
		phidot.switch_state(momentum);
		chidot.switch_state(momentum);
	}

	phi.divby(ics_scale);
	chi.divby(ics_scale);
	phidot.divby(ics_scale);
	chidot.divby(ics_scale);
		
	// Note that the 0-mode in Fourier space is the sum over all points in position space.
	phi.mdata[0][0] = mp.rescale_A * fs.total_gridpoints * mp.phi0;
	phi.mdata[0][1] = 0.;
	phidot.mdata[0][0] = mp.rescale_A * fs.total_gridpoints * mp.phidot0;
	phidot.mdata[0][1] = 0.;
	chi.mdata[0][0] = mp.rescale_A * fs.total_gridpoints * mp.chi0;
	chi.mdata[0][1] = 0.;
	chidot.mdata[0][0] = mp.rescale_A * fs.total_gridpoints * mp.chidot0;
	chidot.mdata[0][1] = 0.;
}

/**
 * @page outputs Output Files
 * All output files generated by SpectRE are placed into a directory named output-YYYYMMDDHHMMSS
 * where YYYY is the current year, etc.
 *
 * @li @ref info_txt
 * @li @ref sf_tsv
 * @li @ref energy_tsv
 * @li @ref stats_tsv
 * @li @ref spectra_tsv
 * @li @ref twoptcorr_tsv
 * @li @ref slices
 *
 */

/**
 * @page info_txt info.txt
 * The info.txt contains a human-readable summary of the run parameters (both physical and numerical).
 */

/**
 * @page sf_tsv sf.tsv
 * sf.tsv is a tab serarated file with the following fields:
 * @li Program time
 * @li Physical time
 * @li a
 * @li H
 */

template <typename R>
void model<R>::evolve(integrator<R> *ig)
{
	int counter = 0;
	spectra_outputter<R> so(fs, mp, ts, phi, chi);
	twoptcorr_outputter<R> tpo(fs, mp, ts, phi, chi);
	stats_outputter<R> sto(fs, mp, ts, phi, chi);
	energy_outputter<R> eo(fs, mp, ts, phi, chi, phidot, chidot);
	ofstream scaleof("sf.tsv");
	scaleof << setprecision(30) << fixed;

	// Make H self-consistent.
	R adot1 = 0, adot_homo = ts.adot;
	const R avg_rho_thrsh = 1e-14;
	int avg_rho_iters = 0;
	do {
		eo.output(true);
		ts.adot = ts.a *
			sqrt( 8./(3. * pow<2>(mp.rescale_A) * pow(ts.a, 2. * mp.rescale_r)) * M_PI * eo.avg_rho_phys);
		if (!avg_rho_iters) adot1 = ts.adot;
		++avg_rho_iters;
	} while (fabs(eo.avg_rho - R(1.0)) > avg_rho_thrsh);
	cout << "Initial adot converged in " << avg_rho_iters << " iterations: <rho> = " << eo.avg_rho
		<< " (homo. delta = " << ts.adot - adot_homo << ", from " << adot_homo << ")" 
		<< " (1st iter delta = " << ts.adot - adot1 << ", from " << adot1 << ")" << endl;

	ig->initialize();

	while (ts.t <= tf) {
		if (af > 0.0 && ts.a > af) {
			break;
		}

		if (counter % scale_interval == 0) {
			scaleof << ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<
				ts.a << "\t" << ts.adot/ts.a << endl;
			scaleof.flush();
		}
		
		if (counter % energy_interval == 0) {
			eo.output();
		}

		if (counter % stats_interval == 0) {
			sto.output();
		}

		if (counter % spectra_interval == 0) {
			so.output();
		}

		if (counter % twoptcorr_interval == 0) {
			tpo.output();
		}

		if (counter % slice_interval == 0) {
			som->output();
		}

		ig->step();

#ifdef HAVE_PRIVATE
		priv.set_sf_info(mp, ts);
		priv.evolve(fs, mp, ts, counter);
#endif

		ts.advance();

		if (counter % screen_interval == 0) {
			cout << ts.t/tf * 100 << " %" << endl;
		}

		++counter;
	}
}

template <typename R>
void model<R>::run()
{
	write_info_file();

	set_initial_conditions();
	
	integrator<R> *ig = use_verlet ?
		(integrator<R> *) new verlet<R>(fs, mp, ts, phi, phidot, chi, chidot)
		: (integrator<R> *) new rk4<R>(fs, mp, ts, phi, phidot, chi, chidot);

	cout << "Beginning field evolution..." << endl;
	evolve(ig);
	delete ig;			
}

template <typename R>
void model<R>::write_info_file()
{
	ofstream info_file("info.txt");
	info_file << setprecision(30);
	info_file << scientific;
	
	info_file << "N: " << fs.n << endl;
	info_file << "final time: " << tf << endl;
	info_file << "gamma_phi: " << mp.gamma_phi << endl;
	info_file << "gamma_chi: " << mp.gamma_chi << endl;
	info_file << "lambda_phi: " << mp.lambda_phi << endl;
	info_file << "lambda_chi: " << mp.lambda_chi << endl;
	info_file << "m_phi: " << mp.m_phi << endl;
	info_file << "m_chi: " << mp.m_chi << endl;
	info_file << "g: " << mp.g << endl;
	info_file << "L: " << mp.len << endl;
	info_file << "phi0: " << mp.phi0 << endl;
	info_file << "chi0: " << mp.chi0 << endl;

	info_file << endl;

	info_file << "rescale A: " << mp.rescale_A << endl;
	info_file << "rescale B: " << mp.rescale_B << endl;
	info_file << "rescale r: " << mp.rescale_r << endl;
	info_file << "rescale s: " << mp.rescale_s << endl;

	info_file << endl;

	info_file << "precision: " << precision_name<R>() << endl;
	info_file << "integrator: " << (use_verlet ? "verlet" : "rk4") << endl;
	info_file << "initial conditions: " <<
		(
			(homo_ic_phi && homo_ic_chi) ? "homogeneous" : (le_init ? "latticeeasy" : "defrost")
		) << endl;
	if (le_init) {
		info_file << "base length scale: " << len0 << endl;
	}
	info_file << "homogeneous initial conditions for phi: " << (homo_ic_phi ? "yes" : "no") << endl;
	info_file << "homogeneous initial conditions for chi: " << (homo_ic_chi ? "yes" : "no") << endl;
	info_file << "initial conditions scale: " << ics_scale << endl;
	info_file << "vary field variance with L" << (vvwl ? "yes" : "no") << endl;

	info_file << "parallelize: " <<
#ifdef _OPENMP
		"yes (" << omp_get_max_threads() << " threads)"
#else
		"no"
#endif
		<< endl;

	info_file << "N pad factor: " << fs.n_pad_factor << endl;
	info_file << "random seed: " << seed << endl;
	
	info_file << endl;

	ts.dt_summary(info_file);

	info_file << endl;

#ifdef HAVE_PRIVATE
	priv.info_file_output(info_file);
	info_file << endl;
#endif
}

// Explicit instantiations
template class model<double>;
#ifdef USE_LD
template class model<long double>;
#endif
