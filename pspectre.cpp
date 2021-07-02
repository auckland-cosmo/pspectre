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

/**
 * @mainpage SpectRE - A Spectral Code for Reheating
 * SpectRE is a pseudo-spectral code for simulating a pair of interacting scalar fields
 * in a self-consistently expanding background. These fields are named phi and chi.
 *
 * @par
 * The time-dependent variable-rescaling scheme from LatticeEasy is used to eliminate the
 * first order term from the equations of motion. The fields can be initialized using either
 * the scheme from LatticeEasy or the scheme from Defrost.
 *
 * @li @ref building
 * @li @ref running
 * @li @ref outputs
 *
 * @section refs References
 * @li Gary Felder and Igor Tkachev. LATTICEEASY: A Program for Lattice Simulations of Scalar
 * Fields in an Expanding Universe. arXiv:hep-ph/0011159v1.
 * http://www.science.smith.edu/departments/Physics/fstaff/gfelder/latticeeasy/
 * @li Andrei V. Frolov. DEFROST: A New Code for Simulating Preheating after Inflation.
 * arXiv:0809.4904v2 [hep-ph]. http://www.sfu.ca/physics/cosmology/defrost/
 */

/**
 * @page building Building
 *
 * @section make Make
 * Building SpectRE requires GNU make. On systems where GNU make is not the
 * system's default make command, GNU make is often called gmake. 
 *
 * @section reqs Requirements
 * SpectRE should build and run on any POSIX-style operating system, and uses OpenMP for
 * shared-memory parallelism. It requires:
 * @li FFTW 3 or Intel's MKL version 10+.
 * @li G++ (the GNU C++ compiler version 4+) or ICC (the Intel C++ compiler).
 *
 * @section targets Targets
 * The following (phony) targets are defined:
 * @li rel - Build the release (optimized) spectre executable. This is the default target.
 * @li profile - Build the optimized profiling executable spectre-pg.
 * @li debug - Build the debug spectre-dbg executable.
 * @li debug-mudflap - Build the mudflap-enabled debug executable spectre-dbg-mf.
 * @li doc - Build the documentation (doxygen and dot required).
 * @li clean - Remove all generated files (including executables) except for the documentation.
 * @li clean-doc - Remove the documentation files.
 * @li clean-all - A combination of clean and clean-doc.
 * @li dist - A combination of clean-all and doc followed by the creation of a source archive.
 *
 * @section vars Variables
 * The make file recognizes the following variables which can be specified on the command line
 * prior to or after the target name(s):
 * @li USE_ICC=yes - Use the Intel C++ compiler instead of the GNU C++ compiler.
 * @li USE_MKL=yes - Use the Intel Math Kernel Libraries intead of FFTW. The MKL FFTW wrapper
 * library is used, which is provided in source form with the MKL installation, and so the MKLROOT
 * environmental variable must be set appropriately.
 * @li USE_LD=yes - Enable long-double support (not supported when using MKL). If the fftwl-wisdom
 * utility exists in a directory in the current search path, then long double support is active by
 * default.
 *
 * @section bexamples Examples
 * To build spectre using g++ and FFTW:
 * @code
 * make
 * @endcode
 *
 * To build spectre using icc and the MKL:
 * @code
 * make USE_ICC=yes USE_MKL=yes
 * @endcode
 *
 * To build spectre-dbg using icc and FFTW:
 * @code
 * make USE_ICC=yes debug
 * @endcode
 *
 * @section compiler Compiler Selection
 * The name of the compiler used can be overridden by setting the GXX variable.
 * By default, this variable has the value g++ or icc. If an executable called
 * g++-4 is found in the current search path, then it is used in preference to g++.
 */

#define _XOPEN_SOURCE 600

#include "field.hpp"
#include "integrator.hpp"
#include "model.hpp"	

#include <cstdlib>
#include <cstring>

#include <vector>

#include <iostream>

#include <climits>
#include <cfloat>
#include <cmath>

#include <fenv.h>
#if defined(__i386__) && defined(__SSE__)
#include <xmmintrin.h>
#endif

#include <unistd.h>
#include <wordexp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

int main(int argc, char *argv[])
{
#if defined(FE_NOMASK_ENV) && !defined(__INTEL_COMPILER)
	fesetenv(FE_NOMASK_ENV);
	fedisableexcept(/* FE_OVERFLOW | */ FE_UNDERFLOW | FE_INEXACT);
#elif defined(__i386__) && defined(__SSE__)
	_MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~(_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO));
#endif

#ifdef _OPENMP
	fftw_init_threads(); //Initialization to use threads for FFTW
	fftw_plan_with_nthreads(omp_get_max_threads());

#ifdef USE_LD
	fftwl_init_threads();
	fftwl_plan_with_nthreads(omp_get_max_threads());
#endif

	cout << "Initialized FFT library with " << omp_get_max_threads() << " threads." << endl;
#endif

	vector<char *> args(argv, argv + argc);
	bool first_line = true;
	wordexp_t we;
	memset(&we, 0, sizeof(we));

	// read in the parameters file and append to the argv array...
	if (argc > 1 && argv[argc-1][0] == '@') {
		args.pop_back();

		ifstream pfl(argv[argc-1]+1);
		if (!pfl) {
			cerr << "Unable to open parameters file: " << (argv[argc-1]+1) << endl;
			exit(1);
		}

		string ws = " \t";
		while (pfl) {
			string line;
			getline(pfl, line);

			if (line.length() < 1) {
				continue;
			}

			size_t wse = line.find_first_not_of(ws);
			if (wse != string::npos) {
				line = line.substr(wse);
			}

			if (line.length() < 1) {
				continue;
			}

			if (line[0] == '#') {
				continue;
			}

			if (wordexp(line.c_str(), &we, (first_line ? 0 : WRDE_APPEND) | WRDE_SHOWERR) != 0) {
				cerr << "Error parsing line: " << line << endl;
				exit(1);
			}

			first_line = false;
		}

		if (we.we_wordc) {
			char **as = we.we_wordv, **ae = we.we_wordv + we.we_wordc;
			while (as[0][0] == '\0') ++as;
			args.insert(++args.begin(), as, ae);
		}
	}

#ifdef USE_LD
	// If the *last* argument is --long, then use long doubles
	if (args.size() > 1 && !strcmp(args[args.size()-1], "--long")) {
		args.pop_back();
		model<long double> mdl(args.size(), &args[0]);
		mdl.run();
	}
	else {
#endif
		model<double> mdl(args.size(), &args[0]);
		mdl.run();
#ifdef USE_LD
	}
#endif

	if (!first_line) {
		wordfree(&we);
	}

	return 0;
}

// Explicit instantiations
template struct model_params<double>;
template struct time_state<double>;
#ifdef USE_LD
template struct model_params<long double>;
template struct time_state<long double>;
#endif

