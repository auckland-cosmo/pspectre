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
#include "rk4.hpp"

using namespace std;

template <typename R>
void rk4<R>::initialize()
{
	// empty
}

// Uses the fourth order Runge-Kutta scheme to evolve field values, the field derivatives, 
// and the scale factor.
template <typename R>
void rk4<R>::step()
{	
	R avg_gradient_phi = 0.0, avg_gradient_chi = 0.0;
	
	integrator<R>::avg_gradients(fs, mp, phi, chi,
		avg_gradient_phi, avg_gradient_chi);

	substep(0.5, phi, chi, phidot, chidot,
		phi1, chi1, phidot1, chidot1,
		ts.a, ts.adot, ts.physical_time, a1, adot1, pt1,
		avg_gradient_phi, avg_gradient_chi);
	substep(0.5, phi1, chi1, phidot1, chidot1,
		phi2, chi2, phidot2, chidot2,
		a1, adot1, pt1, a2, adot2, pt2,
		avg_gradient_phi, avg_gradient_chi);
	substep(1.0, phi2, chi2, phidot2, chidot2,
		phi3, chi3, phidot3, chidot3,
		a2, adot2, pt2, a3, adot3, pt3,
		avg_gradient_phi, avg_gradient_chi);
	
	R da4, dadot4, dpt4;
	substep_scale(1.0, phi3, chi3,
		a3, adot3, pt3, a4, adot4, pt4,
		da4, dadot4, dpt4,
		avg_gradient_phi, avg_gradient_chi);

	nlt.transform(phi3, chi3, a3);

	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);
	phidot.switch_state(momentum, true);
	chidot.switch_state(momentum, true);
	
	phi1.switch_state(momentum, true);
	chi1.switch_state(momentum, true);
	phidot1.switch_state(momentum, true);
	chidot1.switch_state(momentum, true);

	phi2.switch_state(momentum, true);
	chi2.switch_state(momentum, true);
	phidot2.switch_state(momentum, true);
	chidot2.switch_state(momentum, true);

	phi3.switch_state(momentum, true);
	chi3.switch_state(momentum, true);
	phidot3.switch_state(momentum, true);
	chidot3.switch_state(momentum, true);

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

				R dphi4, dchi4, dphidot4, dchidot4;
				R mom2 = pow<2>(mp.dp)*(pow<2>(px) + pow<2>(py) + pow<2>(pz));

				for (int c = 0; c < 2; ++c) {
					mp.derivs(phi3.mdata[idx][c], chi3.mdata[idx][c],
						phidot3.mdata[idx][c], chidot3.mdata[idx][c],
						nlt.chi2phi.mdata[idx][c], nlt.phi2chi.mdata[idx][c],
						mp.lambda_phi != 0 ? nlt.phi3.mdata[idx][c] : 0.0,
						mp.lambda_chi != 0 ? nlt.chi3.mdata[idx][c] : 0.0,
						mp.gamma_phi != 0 ? nlt.phi5.mdata[idx][c] : 0.0,
						mp.gamma_chi != 0 ? nlt.chi5.mdata[idx][c] : 0.0,
						mp.md_e_phi != 0 ? nlt.phi_md.mdata[idx][c] : 0.0,
						mp.md_e_chi != 0 ? nlt.chi_md.mdata[idx][c] : 0.0,
						a3, adot3, dadot4/ts.dt, mom2,
						dphi4, dchi4, dphidot4, dchidot4);
	
					dphi4 *= ts.dt, dchi4 *= ts.dt, dphidot4 *= ts.dt, dchidot4 *= ts.dt;

					R dphi1 = 2.0*(phi1.mdata[idx][c] - phi.mdata[idx][c]);
					R dchi1 = 2.0*(chi1.mdata[idx][c] - chi.mdata[idx][c]);
					
					R dphidot1 = 2.0*(phidot1.mdata[idx][c] - phidot.mdata[idx][c]);
					R dchidot1 = 2.0*(chidot1.mdata[idx][c] - chidot.mdata[idx][c]);

					R dphi2 = 2.0*(phi2.mdata[idx][c] - phi.mdata[idx][c]);
					R dchi2 = 2.0*(chi2.mdata[idx][c] - chi.mdata[idx][c]);
					
					R dphidot2 = 2.0*(phidot2.mdata[idx][c] - phidot.mdata[idx][c]);
					R dchidot2 = 2.0*(chidot2.mdata[idx][c] - chidot.mdata[idx][c]);
				
					R dphi3 = phi3.mdata[idx][c] - phi.mdata[idx][c];
					R dchi3 = chi3.mdata[idx][c] - chi.mdata[idx][c];
					
					R dphidot3 = phidot3.mdata[idx][c] - phidot.mdata[idx][c];
					R dchidot3 = chidot3.mdata[idx][c] - chidot.mdata[idx][c];
					
					phi.mdata[idx][c] += (0.5*dphi1 + dphi2 + dphi3 + 0.5*dphi4)/3.0;
					chi.mdata[idx][c] += (0.5*dchi1 + dchi2 + dchi3 + 0.5*dchi4)/3.0;

					phidot.mdata[idx][c] += (0.5*dphidot1 + dphidot2 + dphidot3 + 0.5*dphidot4)/3.0;
					chidot.mdata[idx][c] += (0.5*dchidot1 + dchidot2 + dchidot3 + 0.5*dchidot4)/3.0;
				}
			}
		}
	}

	// Update scale factor, etc.

	R da1 = 2.0*(a1 - ts.a), da2 = 2.0*(a2 - ts.a),
		da3 = a3 - ts.a;
	R dadot1 = 2.0*(adot1 - ts.adot), dadot2 = 2.0*(adot2 - ts.adot),
		dadot3 = adot3 - ts.adot;
	R dpt1 = 2.0*(pt1 - ts.physical_time), dpt2 = 2.0*(pt2 - ts.physical_time),
		dpt3 = pt3 - ts.physical_time;

	ts.a +=(da1 + 2.*da2 + 2.*da3 + da4)/6.;
	ts.adot +=(dadot1 + 2.*dadot2 + 2.*dadot3 + dadot4)/6.;
	ts.physical_time += (dpt1 + 2.*dpt2 + 2.*dpt3 + dpt4)/6.;
}

template <typename R>
void rk4<R>::substep_scale(R fac, field<R> &phip, field<R> &chip,
	R ap, R adotp, R ptp, R &an, R &adotn, R &ptn, R &dan, R &dadotn, R &dptn,
	R &avg_gradient_phi, R &avg_gradient_chi)
{
	R avg_V = vi.integrate(phip, chip, ap);

	dan = ts.dt*adotp;
	ts.addot = mp.adoubledot(ts.t + (fac - 0.5)*ts.dt, ap, adotp, avg_gradient_phi, avg_gradient_chi, avg_V);
	dadotn = ts.dt*ts.addot;
	dptn = pow(ap, -mp.rescale_s)/mp.rescale_B * ts.dt;

	an = ts.a + fac*dan;
	adotn = ts.adot + fac*dadotn;
	ptn = ts.physical_time + fac*dptn;	
}

template <typename R>
void rk4<R>::substep(R fac, field<R> &phip, field<R> &chip, field<R> &phidotp, field<R> &chidotp,
	field<R> &phin, field<R> &chin, field<R> &phidotn, field<R> &chidotn,
	R ap, R adotp, R ptp, R &an, R &adotn, R &ptn, R &avg_gradient_phi, R &avg_gradient_chi)
{
	R dan, dadotn, dptn;
	substep_scale(fac, phip, chip,
		ap, adotp, ptp, an, adotn, ptn, dan, dadotn, dptn,
		avg_gradient_phi, avg_gradient_chi);

	nlt.transform(phip, chip, ap);

	phip.switch_state(momentum, true);
	chip.switch_state(momentum, true);
	phidotp.switch_state(momentum, true);
	chidotp.switch_state(momentum, true);

	phin.switch_state(momentum, true);
	chin.switch_state(momentum, true);
	phidotn.switch_state(momentum, true);
	chidotn.switch_state(momentum, true);

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

				R dphi, dchi, dphidot, dchidot;
				R mom2 = pow<2>(mp.dp)*(pow<2>(px) + pow<2>(py) + pow<2>(pz));

				for (int c = 0; c < 2; ++c) {
					mp.derivs(phip.mdata[idx][c], chip.mdata[idx][c],
						phidotp.mdata[idx][c], chidotp.mdata[idx][c],
						nlt.chi2phi.mdata[idx][c], nlt.phi2chi.mdata[idx][c],
						mp.lambda_phi != 0 ? nlt.phi3.mdata[idx][c] : 0.0,
						mp.lambda_chi != 0 ? nlt.chi3.mdata[idx][c] : 0.0,
						mp.gamma_phi != 0 ? nlt.phi5.mdata[idx][c] : 0.0,
						mp.gamma_chi != 0 ? nlt.chi5.mdata[idx][c] : 0.0,
						mp.md_e_phi != 0 ? nlt.phi_md.mdata[idx][c] : 0.0,
						mp.md_e_chi != 0 ? nlt.chi_md.mdata[idx][c] : 0.0,
						ap, adotp, dadotn/ts.dt, mom2,
						dphi, dchi, dphidot, dchidot);
	
					dphi *= ts.dt, dchi *= ts.dt, dphidot *= ts.dt, dchidot *= ts.dt;

					phin.mdata[idx][c] = phi.mdata[idx][c] + fac*dphi;
					chin.mdata[idx][c] = chi.mdata[idx][c] + fac*dchi;
					
					phidotn.mdata[idx][c] = phidot.mdata[idx][c] + fac*dphidot;
					chidotn.mdata[idx][c] = chidot.mdata[idx][c] + fac*dchidot;
				}

				mom2 *= (z == 0 || z == fs.n/2) ? 1 : 2;

				total_gradient_phi += mom2*(pow<2>(phin.mdata[idx][0]) + pow<2>(phin.mdata[idx][1]));
				total_gradient_chi += mom2*(pow<2>(chin.mdata[idx][0]) + pow<2>(chin.mdata[idx][1]));
			}
		}
	}

	avg_gradient_phi = total_gradient_phi/pow<2, R>(fs.total_gridpoints);
	avg_gradient_chi = total_gradient_chi/pow<2, R>(fs.total_gridpoints);
}

// Explicit instantiations
template class rk4<double>;
#ifdef USE_LD
template class rk4<long double>;
#endif
