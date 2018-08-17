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

#include "v_integrator.hpp"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

template <typename R>
v_integrator<R>::v_integrator(field_size &fs_, model_params<R> &mp_)
	: fs(fs_), mp(mp_)
{
#ifdef _OPENMP
	int sz = fs.n_pad_factor*fs.n*omp_get_max_threads();
#else
	int sz = fs.n_pad_factor*fs.n;
#endif

	y_integral = new R[sz];
	z_integral = new R[sz];
}

template <typename R>
v_integrator<R>::~v_integrator()
{
	delete [] y_integral;
	delete [] z_integral;
}

// Integrate the potential using Simpson's rule, and assume periodic boundaries.
// Returns the average value.
template <typename R>
R v_integrator<R>::integrate(field<R> &phi, field<R> &chi, R a_t)
{
	R total_V = 0.0;

	phi.switch_state(padded_position, true);
	chi.switch_state(padded_position, true);
	
#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_V)
#endif

	for (int x = 0; x < fs.n_pad_factor*fs.n; ++x) {
#ifdef _OPENMP
		int tsi = fs.n_pad_factor*fs.n*omp_get_thread_num();
#else
		int tsi = 0;
#endif

		y_integral[tsi+x] = 0.0;

		for (int y = 0; y < fs.n_pad_factor*fs.n; ++y) {
			
			z_integral[tsi+y] = 0.0;
			
			for (int z = 0; z < fs.n_pad_factor*fs.n; ++z) {
				int idx = z + phi.pldl*(y + fs.n_pad_factor*fs.n*x);				

				z_integral[tsi+y] += (2 + 2*(z % 2))*mp.V(
					phi.data[idx], chi.data[idx], a_t
				);
			}
			
			y_integral[tsi+x] += (2 + 2*(y % 2))*z_integral[tsi+y];
		}
		
		total_V += (2 + 2*(x % 2))*y_integral[tsi+x];
	}

	// The normalizing factor for Simpson's rule iterated over 3 dimensions.
	return (total_V * 1./(3.*3.*3.)) / fs.total_padded_gridpoints;
}

// Explicit instantiations
template class v_integrator<double>;
#ifdef USE_LD
template class v_integrator<long double>;
#endif
