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
#include <errno.h>

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
 * ./pspectre [-r] [-l [-B <real>]] [-V] [-H <name>[,<name>]*] [-O] [-N <int>] [-P <int>] [-L <real>] [-R <int>] [-o <dir name>] [-t <real>[:<real>]] [-T <real>] [-A <real>] [-p <name>=<value>[,<name>=<value>]*] [-e] [-s <name>[,<name>]*] [-S <name>[=<value>][,<name>[=<value>]]*] [-I <name>=<value>[,<name>=<value>]*] [--long] [@<file name>]
 * @endcode
 *
 * @li -h: Display usage information and exit
 * @li -r: Use the RK4 integrator (default is the Verlet integrator)
 * @li -l: Use LatticeEasy-style initial conditions (default is DEFROST-style initial conditions)
 * @li -B: The base length scale (default is 1.0 to match LatticeEasy)
 * @li -V: Allow the field variance to change with L
 * @li -e: Use power-law expansion
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
 *	monodromy_exp_phi
 *	monodromy_exp_chi
 *	monodromy_scale_phi
 *	monodromy_scale_chi
 *	H0
 *	phi0_slice
 *	chi0_slice
 *	phidot0_slice
 *	chidot0_slice
 *	ics_eff_size
 *	(a0 can be specified when H0 is specified by appending :\<a0\> to the H0 value;
 *	 Hdot0 can be similarly appended for use with power-law background expansion)
 *	(file paths provided for *_slice parameters cannot contain comma characters)
 *	(ics_eff_size is an integer <= N)
 * @endcode
 * @li -s \<name\>[,\<name\>]*: Enable slice output of a variable. Valid variables are:
 * @code
 *	phi
 *	chi
 *	phidot
 *	chidot
 *	V
 *	V_phys
 *	T_phi
 *	T_chi
 *	T_phi_phys
 *	T_chi_phys
 *	G_phi
 *	G_chi
 *	G_phi_phys
 *	G_chi_phys
 *	G_phi_x
 *	G_chi_x
 *	G_phi_phys_x
 *	G_chi_phys_x
 *	G_phi_y
 *	G_chi_y
 *	G_phi_phys_y
 *	G_chi_phys_y
 *	G_phi_z
 *	G_chi_z
 *	G_phi_phys_z
 *	G_chi_phys_z
 *	grad_phi_phys_x
 *	grad_chi_phys_x
 *	grad_phi_phys_y
 *	grad_chi_phys_y
 *	grad_phi_phys_z
 *	grad_chi_phys_z
 *	rho
 *	rho_phys
 *	p
 *	p_phys
 *	gpot
 * @endcode
 * @li -S \<name\>[=\<value\>][,\<name\>[=\<value\>]]*: Set a slice output option value. Valid options are:
 * @code
 *	dim
 *	length
 *	skip
 *	avg
 *	fullprec
 *	(avg and fullprec do not take a value)
 * @endcode
 * @li -I \<name\>=\<value\>[:\<real\>][,\<name\>=\<value\>[:\<real\>]]*: Set an output interval with an optional start time. Valid intervals are:
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
		scale_intervals(scale_interval, 0.0, scale_interval, "t", "scale-factor output interval"),
		energy_intervals(energy_interval, 0.0, energy_interval, "t", "energy output interval"),
		spectra_intervals(spectra_interval, 0.0, spectra_interval, "t", "spectra output interval"),
		screen_intervals(screen_interval, 0.0, screen_interval, "t", "screen output interval"),
		slice_intervals(slice_interval, 0.0, slice_interval, "t", "slice output interval"),
		stats_intervals(stats_interval, 0.0, stats_interval, "t", "stats output interval"),
		twoptcorr_intervals(twoptcorr_interval, 0.0, twoptcorr_interval, "t", "two-pt. corr. output interval"),
		phi("phi"), phidot("phidot"), chi("chi"), chidot("chidot"), gc(0), som(0), ics_scale(1), len0(1.0),
		vvwl(false), af(0.0), external_H0(false), ics_eff_size(0), phidot0pr(0.0), chidot0pr(0.0)
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
		"ics_scale", "monodromy_exp_phi", "monodromy_exp_chi",
		"monodromy_scale_phi", "monodromy_scale_chi",
		"H0", "phi0_slice", "chi0_slice",
		"phidot0_slice", "chidot0_slice", "ics_eff_size", 0
	};

	const char *interval_names[] = {
		"scale", "energy", "spectra",
		"screen", "slice", "stats",
		"twoptcorr", "all", 0
	};

	const char *slice_opt_names[] = {
		"dim", "length", "skip",
		"avg", "fullprec", 0
	};
	
	const char *slice_names[] = {
		"phi", "chi", "phidot", "chidot",
		"V", "V_phys",
		"T_phi", "T_chi",
		"T_phi_phys", "T_chi_phys",
		"G_phi", "G_chi",
		"G_phi_phys", "G_chi_phys",
		"G_phi_x", "G_chi_x",
		"G_phi_phys_x", "G_chi_phys_x",
		"G_phi_y", "G_chi_y",
		"G_phi_phys_y", "G_chi_phys_y",
		"G_phi_z", "G_chi_z",
		"G_phi_phys_z", "G_chi_phys_z",
		"grad_phi_phys_x", "grad_chi_phys_x",
		"grad_phi_phys_y", "grad_chi_phys_y",
		"grad_phi_phys_z", "grad_chi_phys_z",
		"rho", "rho_phys", "p", "p_phys",
		"gpot", 0
	};

	const char *field_names[] = {
		"phi", "chi", 0
	};
	
	bool show_usage = false, help_requested = false;

	bool slice_phi = false, slice_chi = false,
		slice_phidot = false, slice_chidot = false,
		slice_V = false, slice_V_phys = false,
		slice_T_phi = false, slice_T_chi = false,
		slice_T_phi_phys = false, slice_T_chi_phys = false,
		slice_G_phi = false, slice_G_chi = false,
		slice_G_phi_phys = false, slice_G_chi_phys = false,
		slice_G_phi_x = false, slice_G_chi_x = false,
		slice_G_phi_phys_x = false, slice_G_chi_phys_x = false,
		slice_G_phi_y = false, slice_G_chi_y = false,
		slice_G_phi_phys_y = false, slice_G_chi_phys_y = false,
		slice_G_phi_z = false, slice_G_chi_z = false,
		slice_G_phi_phys_z = false, slice_G_chi_phys_z = false,
		slice_grad_phi_phys_x = false, slice_grad_chi_phys_x = false,
		slice_grad_phi_phys_y = false, slice_grad_chi_phys_y = false,
		slice_grad_phi_phys_z = false, slice_grad_chi_phys_z = false,
		slice_rho = false, slice_p = false,
		slice_rho_phys = false, slice_p_phys = false,
		slice_gpot = false;

	int slicedim = 3, slicelength = 0, sliceskip = 1;
	bool sliceaverage = false, sliceflt = true;

	bool oop = false;
	string odn;

	while ((opt = getopt(argc, argv, ":rlVB:hH:ON:P:L:R:p:o:t:T:A:s:S:I:z:e")) != -1) {
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
		case 'e':
			mp.pwr_exp = true;
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
				else if (!strcmp(param_names[index], "monodromy_exp_phi")) {
					mp.md_e_phi = atof(value);
				}
				else if (!strcmp(param_names[index], "monodromy_exp_chi")) {
					mp.md_e_chi = atof(value);
				}
				else if (!strcmp(param_names[index], "monodromy_scale_phi")) {
					mp.md_s_phi = atof(value);
				}
				else if (!strcmp(param_names[index], "H0")) {
					// Parse as H0:a0, where a0 is 1 by default.
					char *en;
					R H0 = strtod(value, &en);
					R a0 = 1.0, Hdot0 = 0.0;
					if (*en != 0) {
						char *en2;
						a0 = strtod(en + 1, &en2);
						if (*en2 != 0) {
							Hdot0 = atof(en2 + 1);
						}
					}

					ts.a = a0;
					ts.adot = H0 * ts.a;
					ts.addot = Hdot0 * ts.a;
					external_H0 = true;
				}
				else if (!strcmp(param_names[index], "phi0_slice")) {
					phi0_slice = std::string(value);
				}
				else if (!strcmp(param_names[index], "chi0_slice")) {
					chi0_slice = std::string(value);
				}
				else if (!strcmp(param_names[index], "phidot0_slice")) {
					phidot0_slice = std::string(value);
				}
				else if (!strcmp(param_names[index], "chidot0_slice")) {
					chidot0_slice = std::string(value);
				}
				else if (!strcmp(param_names[index], "ics_eff_size")) {
					ics_eff_size = atoi(value);
				}
			}
			break;
		case 't':
			// Parse as dt or dt:start_time
			{
				char *en;
				R dt = strtod(optarg, &en);
				R start_time = 0.0;
				if (*en != 0) {
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
				else if (!strcmp(slice_names[index], "phidot")) {
					slice_phidot = true;
				}
				else if (!strcmp(slice_names[index], "chidot")) {
					slice_chidot = true;
				}
				else if (!strcmp(slice_names[index], "V")) {
					slice_V = true;
				}
				else if (!strcmp(slice_names[index], "V_phys")) {
					slice_V_phys = true;
				}
				else if (!strcmp(slice_names[index], "T_phi")) {
					slice_T_phi = true;
				}
				else if (!strcmp(slice_names[index], "T_chi")) {
					slice_T_chi = true;
				}
				else if (!strcmp(slice_names[index], "T_phi_phys")) {
					slice_T_phi_phys = true;
				}
				else if (!strcmp(slice_names[index], "T_chi_phys")) {
					slice_T_chi_phys = true;
				}
				else if (!strcmp(slice_names[index], "G_phi")) {
					slice_G_phi = true;
				}
				else if (!strcmp(slice_names[index], "G_chi")) {
					slice_G_chi = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_phys")) {
					slice_G_phi_phys = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_phys")) {
					slice_G_chi_phys = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_x")) {
					slice_G_phi_x = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_x")) {
					slice_G_chi_x = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_phys_x")) {
					slice_G_phi_phys_x = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_phys_x")) {
					slice_G_chi_phys_x = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_y")) {
					slice_G_phi_y = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_y")) {
					slice_G_chi_y = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_phys_y")) {
					slice_G_phi_phys_y = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_phys_y")) {
					slice_G_chi_phys_y = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_z")) {
					slice_G_phi_z = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_z")) {
					slice_G_chi_z = true;
				}
				else if (!strcmp(slice_names[index], "G_phi_phys_z")) {
					slice_G_phi_phys_z = true;
				}
				else if (!strcmp(slice_names[index], "G_chi_phys_z")) {
					slice_G_chi_phys_z = true;
				}
				else if (!strcmp(slice_names[index], "grad_phi_phys_x")) {
					slice_grad_phi_phys_x = true;
				}
				else if (!strcmp(slice_names[index], "grad_chi_phys_x")) {
					slice_grad_chi_phys_x = true;
				}
				else if (!strcmp(slice_names[index], "grad_phi_phys_y")) {
					slice_grad_phi_phys_y = true;
				}
				else if (!strcmp(slice_names[index], "grad_chi_phys_y")) {
					slice_grad_chi_phys_y = true;
				}
				else if (!strcmp(slice_names[index], "grad_phi_phys_z")) {
					slice_grad_phi_phys_z = true;
				}
				else if (!strcmp(slice_names[index], "grad_chi_phys_z")) {
					slice_grad_chi_phys_z = true;
				}
				else if (!strcmp(slice_names[index], "rho")) {
					slice_rho = true;
				}
				else if (!strcmp(slice_names[index], "rho_phys")) {
					slice_rho_phys = true;
				}
				else if (!strcmp(slice_names[index], "p")) {
					slice_p = true;
				}
				else if (!strcmp(slice_names[index], "p_phys")) {
					slice_p_phys = true;
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
				else if (!strcmp(slice_opt_names[index], "fullprec")) {
					sliceflt = false;
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
				else {
					char *en;
					int iv = (int) strtol(value, &en, 10);
					R start_time = 0.0;
					if (*en != '0') {
						start_time = atof(en + 1);
					}

					if (!strcmp(interval_names[index], "scale")) {
						scale_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "energy")) {
						energy_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "spectra")) {
						spectra_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "screen")) {
						screen_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "slice")) {
						slice_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "stats")) {
						stats_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "twoptcorr")) {
						twoptcorr_intervals.add_value(start_time, iv);
					}
					else if (!strcmp(interval_names[index], "all")) {
						scale_intervals.add_value(start_time, iv);
						energy_intervals.add_value(start_time, iv);
						spectra_intervals.add_value(start_time, iv);
						screen_intervals.add_value(start_time, iv);
						slice_intervals.add_value(start_time, iv);
						stats_intervals.add_value(start_time, iv);
						twoptcorr_intervals.add_value(start_time, iv);
					}
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
			"[-p <name>=<value>[,<name>=<value>]*] [-e] [-s <name>[,<name>]*] [-S <name>[=<value>][,<name>[=<value>]]*] "
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
		hout << "\t-e: Use power-law expansion" << endl;

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
		hout << "\t\t(a0 can be specified when H0 is specified by appending :<a0> to the H0 value" << endl;
		hout << "\t\t Hdot0 can be similarly appended for use with power-law background expansion)" << endl;
		hout << "\t\t(file paths provided for *_slice parameters cannot contain comma characters)" << endl;
		hout << "\t\t(ics_eff_size is an integer <= N)" << endl;

		hout << "\t-s <name>[,<name>]*: Enable slice output of a variable. Valid variables are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(slice_names)/sizeof(slice_names[0])) - 1; ++i) {
			hout << "\t\t" << slice_names[i] << endl;
		}

		hout << "\t-S <name>[=<value>][,<name>[=<value>]]*: Set a slice output option value. Valid options are:" << endl;
		
		for (int i = 0; i < (int) (sizeof(slice_opt_names)/sizeof(slice_opt_names[0])) - 1; ++i) {
			hout << "\t\t" << slice_opt_names[i] << endl;
		}
		hout << "\t\t(avg and fullprec do not take a value)" << endl;

		hout << "\t-I <name>=<value>[:<real>][,<name>=<value>[:<real>]]*: Set an output interval with an optional start time. Valid intervals are:" << endl;
		
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

	scale_intervals.finalize_values();
	energy_intervals.finalize_values();
	spectra_intervals.finalize_values();
	screen_intervals.finalize_values();
	slice_intervals.finalize_values();
	stats_intervals.finalize_values();
	twoptcorr_intervals.finalize_values();

	srand48(seed);
	
	cout << "Using " << precision_name<R>() << " precision." << endl;

	mp.calculate_derived_params(true);
#ifdef DOT0_IN_PLANCK
	mp.phidot0 *= mp.rescale_B;
	mp.chidot0 *= mp.rescale_B;
#endif

	fs.calculate_size_totals();

	phi.construct(fs, oop);
	phidot.construct(fs, oop);
	chi.construct(fs, oop);
	chidot.construct(fs, oop);

	char *swd = 0; std::size_t swdl = 1024;
	do {
		delete [] swd;
		swdl *= 2;
		swd = new char[swdl];
	} while (!getcwd(swd, swdl) && errno == ERANGE);
	start_wd = swd;
	delete [] swd;

	set_output_directory(odn.c_str());

	gc = new grad_computer<R>(fs, mp, phi, chi);
	gpotc = new gpot_computer<R>(fs, mp, ts, phi, chi, phidot, chidot, *gc);
	som = new slice_output_manager<R>(fs, mp, ts, phi, chi, phidot, chidot, *gc, *gpotc,
		slicedim, slicelength, sliceskip, sliceaverage, sliceflt);
	
	if (slice_phi) som->add_outputter("phi", grid_funcs<R>::compute_phi);
	if (slice_chi) som->add_outputter("chi", grid_funcs<R>::compute_chi);
	if (slice_phidot) som->add_outputter("phidot", grid_funcs<R>::compute_phidot);
	if (slice_chidot) som->add_outputter("chidot", grid_funcs<R>::compute_chidot);
	if (slice_V) som->add_outputter("V", grid_funcs<R>::compute_V);
	if (slice_V_phys) som->add_outputter("V_phys", grid_funcs<R>::compute_V_phys);
	if (slice_T_phi) som->add_outputter("T_phi", grid_funcs<R>::compute_T_phi);
	if (slice_T_chi) som->add_outputter("T_chi", grid_funcs<R>::compute_T_chi);
	if (slice_T_phi_phys) som->add_outputter("T_phi_phys", grid_funcs<R>::compute_T_phi_phys);
	if (slice_T_chi_phys) som->add_outputter("T_chi_phys", grid_funcs<R>::compute_T_chi_phys);
	if (slice_G_phi) som->add_outputter("G_phi", grid_funcs<R>::compute_G_phi);
	if (slice_G_chi) som->add_outputter("G_chi", grid_funcs<R>::compute_G_chi);
	if (slice_G_phi_phys) som->add_outputter("G_phi_phys", grid_funcs<R>::compute_G_phi_phys);
	if (slice_G_chi_phys) som->add_outputter("G_chi_phys", grid_funcs<R>::compute_G_chi_phys);
	if (slice_G_phi_x) som->add_outputter("G_phi_x", grid_funcs<R>::compute_G_phi_x);
	if (slice_G_chi_x) som->add_outputter("G_chi_x", grid_funcs<R>::compute_G_chi_x);
	if (slice_G_phi_phys_x) som->add_outputter("G_phi_phys_x", grid_funcs<R>::compute_G_phi_phys_x);
	if (slice_G_chi_phys_x) som->add_outputter("G_chi_phys_x", grid_funcs<R>::compute_G_chi_phys_x);
	if (slice_G_phi_y) som->add_outputter("G_phi_y", grid_funcs<R>::compute_G_phi_y);
	if (slice_G_chi_y) som->add_outputter("G_chi_y", grid_funcs<R>::compute_G_chi_y);
	if (slice_G_phi_phys_y) som->add_outputter("G_phi_phys_y", grid_funcs<R>::compute_G_phi_phys_y);
	if (slice_G_chi_phys_y) som->add_outputter("G_chi_phys_y", grid_funcs<R>::compute_G_chi_phys_y);
	if (slice_G_phi_z) som->add_outputter("G_phi_z", grid_funcs<R>::compute_G_phi_z);
	if (slice_G_chi_z) som->add_outputter("G_chi_z", grid_funcs<R>::compute_G_chi_z);
	if (slice_G_phi_phys_z) som->add_outputter("G_phi_phys_z", grid_funcs<R>::compute_G_phi_phys_z);
	if (slice_G_chi_phys_z) som->add_outputter("G_chi_phys_z", grid_funcs<R>::compute_G_chi_phys_z);
	if (slice_grad_phi_phys_x) som->add_outputter("grad_phi_phys_x", grid_funcs<R>::compute_grad_phi_phys_x);
	if (slice_grad_chi_phys_x) som->add_outputter("grad_chi_phys_x", grid_funcs<R>::compute_grad_chi_phys_x);
	if (slice_grad_phi_phys_y) som->add_outputter("grad_phi_phys_y", grid_funcs<R>::compute_grad_phi_phys_y);
	if (slice_grad_chi_phys_y) som->add_outputter("grad_chi_phys_y", grid_funcs<R>::compute_grad_chi_phys_y);
	if (slice_grad_phi_phys_z) som->add_outputter("grad_phi_phys_z", grid_funcs<R>::compute_grad_phi_phys_z);
	if (slice_grad_chi_phys_z) som->add_outputter("grad_chi_phys_z", grid_funcs<R>::compute_grad_chi_phys_z);
	if (slice_rho) som->add_outputter("rho", grid_funcs<R>::compute_rho);
	if (slice_rho_phys) som->add_outputter("rho_phys", grid_funcs<R>::compute_rho_phys);
	if (slice_p) som->add_outputter("p", grid_funcs<R>::compute_p);
	if (slice_p_phys) som->add_outputter("p_phys", grid_funcs<R>::compute_p_phys);
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

	// Note that the relationship between phidot0 in physical and program units is bit complicated:
	// f_pr = A a^r f => f'_pr = d/dt_pr f_pr = d/dt_pr A a^r f = 1/B a^{-s} d/dt ( A a^r f ) =>
	// f'_pr = A/B a^{-s} d/dt ( a^r f ) = A/B a^{-s} [ a^r f' + r a^{r-1} a' f ] =>
	// f'_pr = A/B [ a^{r-s} f' + r a^{r-1-s} a' f ]
	// So setting a' depends on f'_pr and vice versa, so we'll iterate to convergence...

	if (!external_H0) {
		phidot0pr = mp.rescale_A*mp.phidot0;
		chidot0pr = mp.rescale_A*mp.chidot0;

		const R adot_thrsh = 1e-14;
		R adot_prev;
		int adot_iters = 0;

		ts.adot = 0.0;
		do {
			adot_prev = ts.adot;

			R hf = 3. * pow<2>(mp.rescale_A)/(4. * M_PI);
			R h0norm = 1. / (hf - pow<2>(mp.rescale_r*mp.rescale_A) * (pow<2>(mp.phi0) + pow<2>(mp.chi0)));
			for (int s = -1; s <= 1; s += 2) {
				ts.adot = h0norm * (
					-mp.rescale_r*pow<2>(mp.rescale_A)*((phidot0pr/mp.rescale_A)*mp.phi0 + (chidot0pr/mp.rescale_A)*mp.chi0) +
					s*sqrt(
						hf * (pow<2>(phidot0pr) + pow<2>(chidot0pr)) +
						2. * mp.V(mp.rescale_A * mp.phi0, mp.rescale_A * mp.chi0, 1.) *
							(hf - pow<2>(mp.rescale_r*mp.rescale_A)*(pow<2>(mp.phi0) + pow<2>(mp.chi0)))
					)
				);

				if (ts.adot >= 0) {
					break;				
				}
			}

			// Assuming here that a = 1.
			phidot0pr = mp.rescale_A*(mp.phidot0 + mp.rescale_r*ts.adot*mp.phi0);
			chidot0pr = mp.rescale_A*(mp.chidot0 + mp.rescale_r*ts.adot*mp.chi0);

			++adot_iters;
		} while (adot_iters < 2 || fabs(ts.adot - adot_prev) > ts.adot*adot_thrsh);

		cout << "Initial homogeneous adot (to be corrected later) = " << ts.adot << " (converged in " << adot_iters << " iteration(s))" << endl;
	}
	else {
		phidot0pr = mp.rescale_A*(
			pow(ts.a, mp.rescale_r - mp.rescale_s)*mp.phidot0 +
			mp.rescale_r*pow(ts.a, mp.rescale_r-mp.rescale_s-1)*ts.adot*mp.phi0
		);
		chidot0pr = mp.rescale_A*(
			pow(ts.a, mp.rescale_r - mp.rescale_s)*mp.chidot0 +
			mp.rescale_r*pow(ts.a, mp.rescale_r-mp.rescale_s-1)*ts.adot*mp.chi0
		);
	}

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
	phi.mdata[0][0] = mp.rescale_A * pow(ts.a, mp.rescale_r) * fs.total_gridpoints * mp.phi0;
	phi.mdata[0][1] = 0.;
	phidot.mdata[0][0] = fs.total_gridpoints * phidot0pr;
	phidot.mdata[0][1] = 0.;
	chi.mdata[0][0] = mp.rescale_A * pow(ts.a, mp.rescale_r) * fs.total_gridpoints * mp.chi0;
	chi.mdata[0][1] = 0.;
	chidot.mdata[0][0] = fs.total_gridpoints * chidot0pr;
	chidot.mdata[0][1] = 0.;

	load_initial_slice_file(phi0_slice, phi, pow(ts.a, mp.rescale_r)*mp.rescale_A);
	load_initial_slice_file(chi0_slice, chi, pow(ts.a, mp.rescale_r)*mp.rescale_A);
	load_initial_slice_file(phidot0_slice, phidot, pow(ts.a, mp.rescale_r - mp.rescale_s)*mp.rescale_A);
	load_initial_slice_file(chidot0_slice, chidot, pow(ts.a, mp.rescale_r - mp.rescale_s)*mp.rescale_A);

	if (ics_eff_size > 0) {
		phi.switch_state(momentum);
		chi.switch_state(momentum);
		phidot.switch_state(momentum);
		chidot.switch_state(momentum);

		// Note that F_{x,y,z} = F*_{N-x,N-y,N-z}
		// What does ics_eff_size mean? In means that all momentum modes will be zeroed out which
		// would not be on a grid of size ics_eff_size.
		int effpmax = ics_eff_size/2;
		int effpmin = -effpmax+1;

#ifdef _OPENMP
#pragma omp parallel for
#endif

		for (int x = 0; x < fs.n; ++x) {
			int px = (x <= fs.n/2 ? x : x - fs.n);
			for (int y = 0; y < fs.n; ++y) {
				int py = (y <= fs.n/2 ? y : y - fs.n);
				for (int z = 0; z < fs.n/2+1; ++z) {
					int pz = z;
					if (
						px > effpmax || py > effpmax || pz > effpmax ||
						px < effpmin || py < effpmin || pz < effpmin
					) {
						int idx = z + (fs.n/2+1)*(y + fs.n * x);
						phi.mdata[idx][0] = phi.mdata[idx][1] = R(0);
						chi.mdata[idx][0] = chi.mdata[idx][1] = R(0);
						phidot.mdata[idx][0] = phidot.mdata[idx][1] = R(0);
						chidot.mdata[idx][0] = chidot.mdata[idx][1] = R(0);
					}
				}
			}
		}
	}
}

template <typename R>
void model<R>::load_initial_slice_file(std::string &ifn, field<R> &fld, R pf)
{
	using namespace std;

	if (!ifn.empty()) {
		cout << "Replacing " << fld.name << " field data with slice: " << ifn << endl;
		ifstream ifs((start_wd + "/" + ifn).c_str(), ios::in | ios::binary);
		if (!ifs) {
			cout << "ERROR: unable to open " << ifn << endl;
		}
		else {
			ifs.seekg(0, ios::end);
			std::size_t sz = ifs.tellg();
			ifs.seekg(0);

			if (sz == fld.fs.total_gridpoints*sizeof(R)) {
				fld.switch_state(position);
				ifs.read((char *) fld.data, sz);
				cout << "\tread " << sz << " bytes into the position-space grid" << endl;
				fld.divby(1./pf);
			}
			else if (sz == fld.fs.total_padded_gridpoints*sizeof(R)) {
				fld.switch_state(padded_position);
				ifs.read((char *) fld.data, sz);
				cout << "\tread " << sz << " bytes into the padded position-space grid" << endl;
				fld.divby(1./pf);
			}
			else if (sz < fld.fs.total_gridpoints*sizeof(R)) {
				// Need to interpolate up to the current size...
				int p = 1;
				for (; p < fld.fs.n; ++p) {
					if (sz == fld.fs.total_gridpoints*sizeof(R)/pow<3>(p)) {
						field_size fs(fld.fs.n/p, p);
						field<R> ifld(fs);

						fld.switch_state(position);
						ifld.switch_state(position);
						ifs.read((char *) ifld.data, sz);
						ifld.switch_state(padded_position);
						std::copy(ifld.data, ifld.data + fld.fs.total_gridpoints, fld.data);

						cout << "\tread " << sz << " bytes, interpolating up by a factor of " << p << endl;
						fld.divby(1./pf);

						break;
					}
				}
				if (p == fld.fs.n) {
					cout << "ERROR: no padding factor found for file size: " << sz << endl;
				}
			}
			else {
				// Need to interpolate down to the current size...
				int p = 1;
				for (; p < fld.fs.n; ++p) {
					if (sz == fld.fs.total_gridpoints*sizeof(R)*pow<3>(p)) {
						field_size fs(fld.fs.n/p, p);
						field<R> ifld(fs);

						fld.switch_state(position);
						ifld.switch_state(padded_position);
						ifs.read((char *) ifld.data, sz);
						ifld.switch_state(position);
						std::copy(ifld.data, ifld.data + fld.fs.total_gridpoints, fld.data);

						cout << "\tread " << sz << " bytes, interpolating down by a factor of " << p << endl;
						fld.divby(1./pf);

						break;
					}
				}
				if (p == fld.fs.n) {
					cout << "ERROR: no padding factor found for file size: " << sz << endl;
				}
			}
		}
	}
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
	stats_outputter<R> sto(fs, mp, ts, phi, chi, phidot, chidot);
	energy_outputter<R> eo(fs, mp, ts, phi, chi, phidot, chidot);
	ofstream scaleof("sf.tsv");
	scaleof << setprecision(30) << fixed;

	if (!external_H0) {
		// Make H self-consistent.
		R adot1 = 0, adot_homo = ts.adot;
		const R avg_rho_thrsh = 1e-14;
		int avg_rho_iters = 0;
		do {
			eo.output(true);
			// cout << "avg_rho_phys: " << eo.avg_rho_phys << endl;
			ts.adot = ts.a *
				sqrt( 8./(3. * pow<2>(mp.rescale_A) * pow(ts.a, 2. * mp.rescale_r)) * M_PI * eo.avg_rho_phys);
			if (!avg_rho_iters) adot1 = ts.adot;
			++avg_rho_iters;
		} while (fabs(eo.avg_rho - R(1.0)) > avg_rho_thrsh);
		cout << "Initial adot converged in " << avg_rho_iters << " iteration(s) to " << ts.adot << ": <rho> = " << eo.avg_rho
			<< " (homo. delta = " << ts.adot - adot_homo << ", from " << adot_homo << ")" 
			<< " (1st iter delta = " << ts.adot - adot1 << ", from " << adot1 << ")" << endl;
	}
	else {
		cout << "User-provided adot = " << ts.adot << ", a = " << ts.a << endl;
	}

	if (mp.pwr_exp) {
		cout << "Using power-law background expansion..." << endl;
		// addot = (G-1)/G 1/a adot^2 1/f^2, f = 1 =>
		// a * addot/adot^2 = 1-1/G =>
		// 1 - a * addot/adot^2 = 1/G
		mp.pwr_exp_G = 1./(1. - ts.a*ts.addot/pow<2>(ts.adot));
		cout << "\tG = " << mp.pwr_exp_G << endl;

		// G = gamma/(gamma*s + 1) =>
		// 1/G = (gamma*s + 1)/gamma =>
		// 1/G = s + 1/gamma =>
		// 1/gamma = 1/G - s
		R gamma = 1./(1./mp.pwr_exp_G - mp.rescale_s);
		cout << "Effective power-law exponent: " << gamma << endl;

		// gamma = 2/(3(1+alpha)) =>
		// 2/(3*gamma) = 1 + alpha =>
		// alpha = 2/(3*gamma) - 1
		R alpha = 2./(3.*gamma) - R(1);
		cout << "Effective E.o.S.: p = " << alpha << "*rho" << endl;
	}

	ig->initialize();

	while (ts.t <= tf) {
		if (af > 0.0 && ts.a > af) {
			cout << "Exiting because the scale factor is now " << ts.a << endl;
			break;
		}

		if (counter % scale_interval == 0) {
			scaleof << ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<
				ts.a << "\t" << ts.adot/ts.a << "\t" << ts.addot/ts.a << endl;
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
		scale_intervals.advance(ts.t);
		energy_intervals.advance(ts.t);
		spectra_intervals.advance(ts.t);
		screen_intervals.advance(ts.t);
		slice_intervals.advance(ts.t);
		stats_intervals.advance(ts.t);
		twoptcorr_intervals.advance(ts.t);

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
	info_file << "monodromy_exp_phi: " << mp.md_e_phi << endl;
	info_file << "monodromy_exp_chi: " << mp.md_e_chi << endl;
	info_file << "monodromy_scale_phi: " << mp.md_s_phi << endl;
	info_file << "monodromy_scale_chi: " << mp.md_s_chi << endl;
	info_file << "L: " << mp.len << endl;
	info_file << "phi0: " << mp.phi0 << endl;
	info_file << "chi0: " << mp.chi0 << endl;
	info_file << "phidot0: " << mp.phidot0 << endl;
	info_file << "chidot0: " << mp.chidot0 << endl;

	info_file << endl;

	info_file << "rescale A: " << mp.rescale_A << endl;
	info_file << "rescale B: " << mp.rescale_B << endl;
	info_file << "rescale r: " << mp.rescale_r << endl;
	info_file << "rescale s: " << mp.rescale_s << endl;
	if (mp.pwr_exp) {
		info_file << "power-law expansion: yes" <<  endl;
	}

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
	info_file << "effective-size cutoff: "; if (ics_eff_size > 0) info_file << ics_eff_size; else info_file << "none"; info_file << endl;
	info_file << "vary field variance with L: " << (vvwl ? "yes" : "no") << endl;

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
