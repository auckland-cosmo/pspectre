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
#include "gpot_computer.hpp"
#include "grid_funcs.hpp"

#include <cmath>

using namespace std;

template <typename R>
void gpot_computer<R>::compute(field_state final_state, bool grad_computed)
{
	if (!grad_computed) gc.compute();

	gc.phigradx.switch_state(position);
	gc.phigrady.switch_state(position);
	gc.phigradz.switch_state(position);
	gc.chigradx.switch_state(position);
	gc.chigrady.switch_state(position);
	gc.chigradz.switch_state(position);

	phi.switch_state(position, true);
	chi.switch_state(position, true);
	
	phidot.switch_state(position, true);
	chidot.switch_state(position, true);
	
	gpot.switch_state(uninitialized);
	gpot.switch_state(position);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n; ++x) {
		for (int y = 0; y < fs.n; ++y) {
			for (int z = 0; z < fs.n; ++z) {
				int fdx = z + phi.ldl*(y + fs.n*x);
				int idx = z + gpot.ldl*(y + fs.n*x);
				gpot.data[idx] = grid_funcs<R>::compute_rho_phys(
					fs, mp, ts, phi.data[fdx], chi.data[fdx],
					phidot.data[fdx], chidot.data[fdx],
					gc.phigradx.data[idx], gc.chigradx.data[idx],
					gc.phigrady.data[idx], gc.chigrady.data[idx],
					gc.phigradz.data[idx], gc.chigradz.data[idx], 0.0);
			}
		}
	}
	
	gpot.switch_state(momentum);
	
	R w = 2.0*M_PI/fs.n;

	R c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -3.0, cc = 0.5;
	// R c3 = 0.0, c2 = 1.0, c1 = 1.0, c0 = -6.0, cc = 1.5;
	// R c3 = 1.0, c2 = 0.0, c1 = 2.0, c0 = -7.0, cc = 1.5;
	// R c3 = 2.0, c2 = 3.0, c1 = 7.0, c0 = -32.0, cc = 7.5;

	R dx = mp.len/fs.n;
	R c = cc * 4.0 * M_PI * pow<2>(dx);

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n; ++x) {
		R ii = cos(w*x);
		for (int y = 0; y < fs.n; ++y) {
			R jj = cos(w*y);
			for (int z = 0; z < fs.n/2+1; ++z) {
				R kk = cos(w*z);
				int idx = z+(fs.n/2+1)*(y + fs.n*x);

				R ffd = c0 + c1*(ii+jj+kk) + c2*(ii*jj+ii*kk+jj*kk) + c3*ii*jj*kk;
				R ff = (ffd == 0.0) ? 0.0 : c/ffd;

				gpot.mdata[idx][0] *= ff;
				gpot.mdata[idx][1] *= ff;
			}
		}
	}

	gpot.switch_state(final_state);
}

// Explicit instantiations
template class gpot_computer<double>;
#ifdef USE_LD
template class gpot_computer<long double>;
#endif
