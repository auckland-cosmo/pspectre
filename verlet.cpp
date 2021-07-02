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
#include "verlet.hpp"

using namespace std;

template <typename R>
void verlet<R>::initialize()
{
	R avg_gradient_phi = 0.0, avg_gradient_chi = 0.0;
	
	integrator<R>::avg_gradients(fs, mp, phi, chi,
		avg_gradient_phi, avg_gradient_chi);

	R avg_V = vi.integrate(phi, chi, ts.a);

	ts.addot = addot = mp.adoubledot(ts.t, ts.a, ts.adot, avg_gradient_phi, avg_gradient_chi, avg_V);
	ddptdt = -mp.rescale_s/mp.rescale_B * pow(ts.a, -mp.rescale_s - 1.) * ts.adot;

	nlt.transform(phi, chi, ts.a);

	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);
	phidot.switch_state(momentum, true);
	chidot.switch_state(momentum, true);
	
	phiddot.switch_state(momentum);
	chiddot.switch_state(momentum);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n*x);

				R dphidt, dchidt, dphidotdt, dchidotdt;
				R mom2 = pow<2>(mp.dp)*(pow<2>(px) + pow<2>(py) + pow<2>(pz));

				for (int c = 0; c < 2; ++c) {
					mp.derivs(phi.mdata[idx][c], chi.mdata[idx][c],
						phidot.mdata[idx][c], chidot.mdata[idx][c],
						nlt.chi2phi.mdata[idx][c], nlt.phi2chi.mdata[idx][c],
						mp.lambda_phi != 0 ? nlt.phi3.mdata[idx][c] : 0.0,
						mp.lambda_chi != 0 ? nlt.chi3.mdata[idx][c] : 0.0,
						mp.gamma_phi != 0 ? nlt.phi5.mdata[idx][c] : 0.0,
						mp.gamma_chi != 0 ? nlt.chi5.mdata[idx][c] : 0.0,
						mp.md_e_phi != 0 ? nlt.phi_md.mdata[idx][c] : 0.0,
						mp.md_e_chi != 0 ? nlt.chi_md.mdata[idx][c] : 0.0,
						ts.a, ts.adot, addot, mom2,
						dphidt, dchidt, dphidotdt, dchidotdt);
				
					phiddot.mdata[idx][c] = dphidotdt;
					chiddot.mdata[idx][c] = dchidotdt;
				}
			}
		}
	}
}

template <typename R>
void verlet<R>::step()
{
	adot_staggered = ts.adot + 0.5 * addot * ts.dt;
	dptdt_staggered = dptdt + 0.5 * ddptdt * ts.dt;

	ts.a += ts.adot * ts.dt + 0.5 * addot * pow<2>(ts.dt);
	ts.physical_time += dptdt * ts.dt + 0.5 * ddptdt * pow<2>(ts.dt);

	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);
	phidot.switch_state(momentum, true);
	chidot.switch_state(momentum, true);
	
	phiddot.switch_state(momentum);
	chiddot.switch_state(momentum);
	phidot_staggered.switch_state(momentum);
	chidot_staggered.switch_state(momentum);

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

				for (int c = 0; c < 2; ++c) {

					phidot_staggered.mdata[idx][c] = phidot.mdata[idx][c] + 0.5 * phiddot.mdata[idx][c] * ts.dt;
					chidot_staggered.mdata[idx][c] = chidot.mdata[idx][c] + 0.5 * chiddot.mdata[idx][c] * ts.dt;

					phi.mdata[idx][c] += phidot_staggered.mdata[idx][c]*ts.dt;
					chi.mdata[idx][c] += chidot_staggered.mdata[idx][c]*ts.dt;
				}

				mom2 *= (z == 0 || z == fs.n/2) ? 1 : 2;

				total_gradient_phi += mom2*(pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]));
				total_gradient_chi += mom2*(pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]));
			}
		}
	}

	R avg_gradient_phi = total_gradient_phi/pow<2, R>(fs.total_gridpoints);
	R avg_gradient_chi = total_gradient_chi/pow<2, R>(fs.total_gridpoints);

	R avg_V = vi.integrate(phi, chi, ts.a);

	ts.addot = addot = mp.adoubledot_staggered(ts.t, ts.dt, ts.a, adot_staggered, avg_gradient_phi, avg_gradient_chi, avg_V);
	ts.adot = adot_staggered + 0.5 * addot * ts.dt;

	ddptdt = -mp.rescale_s / mp.rescale_B * pow(ts.a, -mp.rescale_s - 1) * ts.adot;
	dptdt = dptdt_staggered + 0.5 * ddptdt * ts.dt;

	nlt.transform(phi, chi, ts.a);

	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n*x);

				R dphidt, dchidt, dphidotdt, dchidotdt;
				R mom2 = pow<2>(mp.dp)*(pow<2>(px) + pow<2>(py) + pow<2>(pz));

				for (int c = 0; c < 2; ++c) {
					mp.derivs(phi.mdata[idx][c], chi.mdata[idx][c],
						phidot.mdata[idx][c], chidot.mdata[idx][c],
						nlt.chi2phi.mdata[idx][c], nlt.phi2chi.mdata[idx][c],
						mp.lambda_phi != 0 ? nlt.phi3.mdata[idx][c] : 0.0,
						mp.lambda_chi != 0 ? nlt.chi3.mdata[idx][c] : 0.0,
						mp.gamma_phi != 0 ? nlt.phi5.mdata[idx][c] : 0.0,
						mp.gamma_chi != 0 ? nlt.chi5.mdata[idx][c] : 0.0,
						mp.md_e_phi != 0 ? nlt.phi_md.mdata[idx][c] : 0.0,
						mp.md_e_chi != 0 ? nlt.chi_md.mdata[idx][c] : 0.0,
						ts.a, ts.adot /* adot_staggered */, addot, mom2,
						/* Verlet assumes that \ddot{x} is independent of \dot{x},
 						 * using adot here instead of adot_staggered turns out to be better.
 						 */
						dphidt, dchidt, dphidotdt, dchidotdt);
				
					phiddot.mdata[idx][c] = dphidotdt;
					chiddot.mdata[idx][c] = dchidotdt;

					phidot.mdata[idx][c] = phidot_staggered.mdata[idx][c] + 0.5 * phiddot.mdata[idx][c] * ts.dt;
					chidot.mdata[idx][c] = chidot_staggered.mdata[idx][c] + 0.5 * chiddot.mdata[idx][c] * ts.dt;
				}
			}
		}
	}
}

// Explicit instantiations
template class verlet<double>;
#ifdef USE_LD
template class verlet<long double>;
#endif
