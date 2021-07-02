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
#include "fft.hpp"
#include "defrost_style_initializer.hpp"

#include <cstdlib>
#include <cmath>

using namespace std;

template <typename R>
void defrost_style_initializer<R>::sample_grf(field<R> &fld, R gamma, R m2eff)
{
        const int os = 16, nos = fs.n*pow<2>(os), nn = fs.n/2+1;
        const R dx = mp.len/((R) fs.n), dk = 2.0*M_PI/(fs.n*dx);
#ifdef DONT_FIX_DEFROST_ICS
#define KCUT_PF 1.0
#else
#define KCUT_PF 2.0
#endif
        const R dxos = dx/os, dkos = dk/(2.0*os), kcut = KCUT_PF*nn*dk/2.0;
	R norm = 0.5/(sqrt(2*M_PI*pow<3>(dk)) * (2.e5/sqrt(8*M_PI)))*(dkos/dxos);

	if (m2eff < 0) return;

        R *ker = fft_malloc<R>(sizeof(R)*nos);

        for (int k = 0; k < nos; ++k) {
                R kk = (k+0.5)*dkos;
                ker[k] = kk*pow(pow<2>(kk) + m2eff, gamma) * exp(-pow<2>(kk/kcut));
        }

	fft_r2r_1d_plan<R> plan(nos, ker, ker, rodft10);
	plan.execute();

        for (int k = 0; k < nos; ++k) {
                ker[k] *= norm/(k+1);
        }

	fld.switch_state(position);

        for (int x = 0; x < fs.n; ++x)
        for (int y = 0; y < fs.n; ++y)
        for (int z = 0; z < fs.n; ++z) {
		int idx = z + fld.ldl*(y + fs.n*x);
                R kk = sqrt(pow<2>(x+1-nn) + pow<2>(y+1-nn) + pow<2>(z+1-nn)) * os;

                int kkl = (int) floor(kk);

                if (kkl > 0) {
                        fld.data[idx] = ker[kkl-1] + (kk-kkl)*(ker[kkl]-ker[kkl-1]);
                }
                else {
                        fld.data[idx] = (4.0*ker[0]-ker[1])/3.0;
                }
        }

        fft_free<R>(ker);
	fld.switch_state(momentum);

        for (int idx = 0; idx < fld.fs.total_momentum_gridpoints; ++idx) {
                R a = drand48(), p = drand48();
                R rr = mp.rescale_A*sqrt(-2.0*log(a))/(8.0*M_PI);

                typename field<R>::complex_t val;
                val[0] = fld.mdata[idx][0];
                val[1] = fld.mdata[idx][1];

                fld.mdata[idx][0] = rr*(cos(2.0*M_PI*p)*val[0] - sin(2.0*M_PI*p)*val[1]);
                fld.mdata[idx][1] = rr*(cos(2.0*M_PI*p)*val[1] + sin(2.0*M_PI*p)*val[0]);
        }
}

template <typename R>
void defrost_style_initializer<R>::initialize()
{
	// TODO: update effective mass for other potential parameters.
	sample_grf(phi, -0.25, (pow<2>(mp.g*mp.chi0) + pow<2>(mp.m_phi))/pow<2>(mp.rescale_B) - 2.25*pow<2>(adot));
	sample_grf(chi, -0.25, (pow<2>(mp.g*mp.phi0) + pow<2>(mp.m_chi))/pow<2>(mp.rescale_B) - 2.25*pow<2>(adot));

	sample_grf(phidot, 0.25, (pow<2>(mp.g*mp.chi0) + pow<2>(mp.m_phi))/pow<2>(mp.rescale_B) - 2.25*pow<2>(adot));
	sample_grf(chidot, 0.25, (pow<2>(mp.g*mp.phi0) + pow<2>(mp.m_chi))/pow<2>(mp.rescale_B) - 2.25*pow<2>(adot));
}

template class defrost_style_initializer<double>;
#ifdef USE_LD
template class defrost_style_initializer<long double>;
#endif
