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

#include "grad_computer.hpp"

using namespace std;

template <typename R>
void grad_computer<R>::compute(field_state final_state)
{
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	phigradx.switch_state(uninitialized);
	phigradx.switch_state(momentum);
	
	chigradx.switch_state(uninitialized);
	chigradx.switch_state(momentum);
	
	phigrady.switch_state(uninitialized);
	phigrady.switch_state(momentum);
	
	chigrady.switch_state(uninitialized);
	chigrady.switch_state(momentum);

	phigradz.switch_state(uninitialized);
	phigradz.switch_state(momentum);
	
	chigradz.switch_state(uninitialized);
	chigradz.switch_state(momentum);
	
#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n * x);

				// The Fourier transform conventions require the normalization of the position-space values by 1/N^3.
				phigradx.mdata[idx][0] = -1. * px * mp.dp * phi.mdata[idx][1];
				phigradx.mdata[idx][1] = px * mp.dp * phi.mdata[idx][0];

				chigradx.mdata[idx][0] = -1. * px * mp.dp * chi.mdata[idx][1];
				chigradx.mdata[idx][1] = px * mp.dp * chi.mdata[idx][0];

				phigrady.mdata[idx][0] = -1. * py * mp.dp * phi.mdata[idx][1];
				phigrady.mdata[idx][1] = py * mp.dp * phi.mdata[idx][0];

				chigrady.mdata[idx][0] = -1. * py * mp.dp * chi.mdata[idx][1];
				chigrady.mdata[idx][1] = py * mp.dp * chi.mdata[idx][0];

				phigradz.mdata[idx][0] = -1. * pz * mp.dp * phi.mdata[idx][1];
				phigradz.mdata[idx][1] = pz * mp.dp * phi.mdata[idx][0];

				chigradz.mdata[idx][0] = -1. * pz * mp.dp * chi.mdata[idx][1];
				chigradz.mdata[idx][1] = pz * mp.dp * chi.mdata[idx][0];
			}
		}
	}

	phigradx.switch_state(final_state);
	chigradx.switch_state(final_state);
	
	phigrady.switch_state(final_state);
	chigrady.switch_state(final_state);

	phigradz.switch_state(final_state);
	chigradz.switch_state(final_state);	
}

// Explicit instantiations
template class grad_computer<double>;
#ifdef USE_LD
template class grad_computer<long double>;
#endif
