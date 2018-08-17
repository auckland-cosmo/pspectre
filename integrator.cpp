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

#include "pow/pow.hpp"
#include "integrator.hpp"

using namespace std;

template <typename R>
void integrator<R>::avg_gradients(field_size &fs, model_params<R> &mp,
	field<R> &phi, field<R> &chi,
	R &avg_gradient_phi, R &avg_gradient_chi)
{
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	R total_gradient_phi = 0.0, total_gradient_chi = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_gradient_phi,total_gradient_chi)
#endif
	
	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n*x);

				R mom2 = pow<2>(mp.dp)*(pow<2>(px) + pow<2>(py) + pow<2>(pz));
				mom2 *= (z == 0 || z == fs.n/2) ? 1 : 2;
				
				total_gradient_phi += mom2*(pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]));
				total_gradient_chi += mom2*(pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]));
			}
		}
	}

	// Divide by total_gridpoints again to get *average* squared gradient and *average* potential energy.
	avg_gradient_phi = total_gradient_phi/pow<2, R>(fs.total_gridpoints);
	avg_gradient_chi = total_gradient_chi/pow<2, R>(fs.total_gridpoints);
}

// Explicit instantiations
template class integrator<double>;
#ifdef USE_LD
template class integrator<long double>;
#endif
