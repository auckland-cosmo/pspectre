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

#include "nonlinear_transformer.hpp"

using namespace std;

template <typename R>
void nonlinear_transformer<R>::transform(field<R> &phi, field<R> &chi, R a_t, field_state final_state)
{
	phi.switch_state(position, true);
	chi.switch_state(position, true);

	phi2chi.switch_state(uninitialized);
	phi2chi.switch_state(position);
	
	chi2phi.switch_state(uninitialized);
	chi2phi.switch_state(position);
	
	if (mp.lambda_phi != 0.0) {
		phi3.switch_state(uninitialized);
		phi3.switch_state(position);
	}
	
	if (mp.lambda_chi != 0.0) {
		chi3.switch_state(uninitialized);
		chi3.switch_state(position);
	}

	if (mp.gamma_phi != 0.0) {
		phi5.switch_state(uninitialized);
		phi5.switch_state(position);
	}
	
	if (mp.gamma_chi != 0.0) {
		chi5.switch_state(uninitialized);
		chi5.switch_state(position);
	}

	if (mp.md_e_phi != 0.0) {
		phi_md.switch_state(uninitialized);
		phi_md.switch_state(position);
	}

	if (mp.md_e_chi != 0.0) {
		chi_md.switch_state(uninitialized);
		chi_md.switch_state(position);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n; ++x) {
		for (int y = 0; y < fs.n; ++y) {
			for (int z = 0; z < fs.n; ++z) {
				int fdx = z + phi.ldl*(y + fs.n*x);
				int idx = z + phi2chi.ldl*(y + fs.n*x);

				R p = phi.data[fdx];
				R c = chi.data[fdx];
				
				phi2chi.data[idx] = pow<2>(p)*c;
				chi2phi.data[idx] = pow<2>(c)*p;
				
				if (mp.lambda_phi != 0.0) {
					phi3.data[idx] = pow<3>(p);
				}

				if (mp.lambda_chi != 0.0) {
					chi3.data[idx] = pow<3>(c);
				}

				if (mp.gamma_phi != 0.0) {
					phi5.data[idx] = pow<5>(p);
				}

				if (mp.gamma_chi != 0.0) {
					chi5.data[idx] = pow<5>(c);
				}

				if (mp.md_e_phi != 0.0) {
					phi_md.data[idx] = 2.0*mp.md_c_phi*mp.md_e_phi*p *
						pow(1.0 + pow(a_t, -2. * mp.rescale_r)*pow<2>(p/mp.rescale_A)/pow<2>(mp.md_s_phi), mp.md_e_phi - 1.0);
				}

				if (mp.md_e_chi != 0.0) {
					chi_md.data[idx] = 2.0*mp.md_c_chi*mp.md_e_chi*p *
						pow(1.0 + pow(a_t, -2. * mp.rescale_r)*pow<2>(p/mp.rescale_A)/pow<2>(mp.md_s_chi), mp.md_e_chi - 1.0);
				}
			}
		}
	}
	
	phi2chi.switch_state(final_state);
	chi2phi.switch_state(final_state);
	
	if (mp.lambda_phi != 0.0) {
		phi3.switch_state(final_state);
	}
	
	if (mp.lambda_chi != 0.0) {
		chi3.switch_state(final_state);
	}

	if (mp.gamma_phi != 0.0) {
		phi5.switch_state(final_state);
	}
	
	if (mp.gamma_chi != 0.0) {
		chi5.switch_state(final_state);
	}

	if (mp.md_e_phi != 0.0) {
		phi_md.switch_state(final_state);
	}

	if (mp.md_e_chi != 0.0) {
		chi_md.switch_state(final_state);
	}
}

// Explicit instantiations
template class nonlinear_transformer<double>;
#ifdef USE_LD
template class nonlinear_transformer<long double>;
#endif
