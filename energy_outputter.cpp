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
#include "energy_outputter.hpp"
#include "grid_funcs.hpp"

#include <sstream>
#include <iomanip>

using namespace std;

template <typename R>
energy_outputter<R>::energy_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
	field<R> &phi_, field<R> &chi_, field<R> &phidot_, field<R> &chidot_)
	: fs(fs_), mp(mp_), ts(ts_), phi(phi_), chi(chi_), phidot(phidot_), chidot(chidot_),
	vi(fs, mp), avg_rho_phys(0.0), avg_rho(0.0)
{
	of.open("energy.tsv");
	of << setprecision(30) << scientific;
}

template <typename R>
void energy_outputter<R>::output(bool no_output)
{
	R avg_V = vi.integrate(phi, chi, ts.a);
	
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	phidot.switch_state(momentum, true);
	chidot.switch_state(momentum, true);
	
	R avg_phi_squared = 0.0, avg_chi_squared = 0.0,
		avg_phidot_squared = 0.0, avg_chidot_squared = 0.0,
		avg_gradient_phi_x = 0.0, avg_gradient_chi_x = 0.0,
		avg_gradient_phi_y = 0.0, avg_gradient_chi_y = 0.0,
		avg_gradient_phi_z = 0.0, avg_gradient_chi_z = 0.0,
		avg_ffd_phi = 0.0, avg_ffd_chi = 0.0;

		int tc = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tc,avg_phi_squared,avg_chi_squared,avg_phidot_squared,avg_chidot_squared,avg_gradient_phi_x,avg_gradient_chi_x,avg_gradient_phi_y,avg_gradient_chi_y,avg_gradient_phi_z,avg_gradient_chi_z,avg_ffd_phi,avg_ffd_chi)
#endif

	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1) * (y + fs.n * x);
				int cnt = (z == 0 || z == fs.n/2) ? 1 : 2;

				tc += cnt;
				R phi_squared = pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]);
				R chi_squared = pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]);
				avg_phi_squared += cnt * phi_squared;
				avg_chi_squared += cnt * chi_squared;
				// cout << x << " " << y << " " << z << " " << idx << " " << phi_squared << " " << cnt * phi_squared << " " << avg_phi_squared << endl;
				
				avg_phidot_squared += cnt * (pow<2>(phidot.mdata[idx][0]) + pow<2>(phidot.mdata[idx][1]));
				avg_chidot_squared += cnt * (pow<2>(chidot.mdata[idx][0]) + pow<2>(chidot.mdata[idx][1]));
				
				avg_ffd_phi += cnt * (
					phi.mdata[idx][0] * phidot.mdata[idx][0] +
					phi.mdata[idx][1] * phidot.mdata[idx][1]
				);
				avg_ffd_chi += cnt * (
					chi.mdata[idx][0] * chidot.mdata[idx][0] +
					chi.mdata[idx][1] * chidot.mdata[idx][1]
				);
				
				R mom2x = pow<2>(mp.dp) * pow<2>(px);
				R mom2y = pow<2>(mp.dp) * pow<2>(py);
				R mom2z = pow<2>(mp.dp) * pow<2>(pz);
				avg_gradient_phi_x += cnt * mom2x * phi_squared;
				avg_gradient_chi_x += cnt * mom2x * chi_squared;
				avg_gradient_phi_y += cnt * mom2y * phi_squared;
				avg_gradient_chi_y += cnt * mom2y * chi_squared;
				avg_gradient_phi_z += cnt * mom2z * phi_squared;
				avg_gradient_chi_z += cnt * mom2z * chi_squared;
			}
		}
	}

	// The first factor of 1./N^3 comes from Parseval's theorem.
	avg_phidot_squared /= 2*pow<2, R>(fs.total_gridpoints);
	avg_chidot_squared /= 2*pow<2, R>(fs.total_gridpoints);

	R fld_fac = 0.5 * pow<2>(mp.rescale_r) * pow<2>(ts.adot/ts.a);
	avg_phi_squared *= fld_fac / pow<2, R>(fs.total_gridpoints);
	avg_chi_squared *= fld_fac / pow<2, R>(fs.total_gridpoints);

	R grad_fac = 0.5 * pow(ts.a, -2. * mp.rescale_s - 2.);
	avg_gradient_phi_x *= grad_fac / pow<2, R>(fs.total_gridpoints);
	avg_gradient_chi_x *= grad_fac / pow<2, R>(fs.total_gridpoints);
	avg_gradient_phi_y *= grad_fac / pow<2, R>(fs.total_gridpoints);
	avg_gradient_chi_y *= grad_fac / pow<2, R>(fs.total_gridpoints);
	avg_gradient_phi_z *= grad_fac / pow<2, R>(fs.total_gridpoints);
	avg_gradient_chi_z *= grad_fac / pow<2, R>(fs.total_gridpoints);
	
	R ffd_fac = -mp.rescale_r * ts.adot/ts.a;
	avg_ffd_phi *= ffd_fac / pow<2, R>(fs.total_gridpoints);
	avg_ffd_chi *= ffd_fac / pow<2, R>(fs.total_gridpoints);

	// This is the *average* energy per gridpoint.
	avg_rho_phys = avg_V + avg_phi_squared + avg_chi_squared +
		avg_phidot_squared + avg_chidot_squared +
		avg_gradient_phi_x + avg_gradient_chi_x +
		avg_gradient_phi_y + avg_gradient_chi_y +
		avg_gradient_phi_z + avg_gradient_chi_z +
		avg_ffd_phi + avg_ffd_chi;

	R avg_p_phys = -avg_V + avg_phi_squared + avg_chi_squared +
		avg_phidot_squared + avg_chidot_squared -
		(
			avg_gradient_phi_x + avg_gradient_chi_x +
			avg_gradient_phi_y + avg_gradient_chi_y +
			avg_gradient_phi_z + avg_gradient_chi_z
		)/3 +
		avg_ffd_phi + avg_ffd_chi;
	R avg_w = avg_p_phys/avg_rho_phys;

	const R es = grid_funcs<R>::compute_energy_scaling(mp, ts);
	avg_rho = es * avg_rho_phys;

	if (!no_output) {
		of << setw(10) << setfill('0') <<
			ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<
			avg_rho_phys << "\t" << avg_rho << "\t" <<
			es * avg_phidot_squared << "\t" <<
			es * avg_chidot_squared << "\t" <<
			es * avg_ffd_phi << "\t" <<
			es * avg_ffd_chi << "\t" <<
			es * avg_phi_squared << "\t" <<
			es * avg_chi_squared << "\t" <<
			es * (avg_gradient_phi_x + avg_gradient_phi_y + avg_gradient_phi_z) << "\t" <<
			es * (avg_gradient_chi_x + avg_gradient_chi_y + avg_gradient_chi_z) << "\t" <<
			es * avg_V << "\t" <<
			avg_phidot_squared << "\t" <<
			avg_chidot_squared << "\t" <<
			avg_ffd_phi << "\t" <<
			avg_ffd_chi << "\t" <<
			avg_phi_squared << "\t" <<
			avg_chi_squared << "\t" <<
			avg_gradient_phi_x + avg_gradient_phi_y + avg_gradient_phi_z << "\t" <<
			avg_gradient_chi_x + avg_gradient_chi_y + avg_gradient_chi_z << "\t" <<
			avg_V << "\t" <<

			avg_gradient_phi_x << "\t" <<
			avg_gradient_chi_x << "\t" <<
			avg_gradient_phi_y << "\t" <<
			avg_gradient_chi_y << "\t" <<
			avg_gradient_phi_z << "\t" <<
			avg_gradient_chi_z << "\t" <<

			avg_p_phys << "\t" <<
			avg_w <<

			endl;
		of.flush();
	}
}

/**
 * @page energy_tsv energy.tsv
 * energy.tsv is a tab serarated file with the following fields:
 * @li Program time
 * @li Physical time
 * @li Average physical energy (w.r.t. the rescaled length)
 * @li Average energy normalized by the Friedmann equation
 * @li Average normalized @f$ \phi'^2 @f$ energy contribution
 * @li Average normalized @f$ \chi'^2 @f$ energy contribution
 * @li Average normalized @f$ \phi\phi' @f$ energy contribution
 * @li Average normalized @f$ \chi\chi' @f$ energy contribution
 * @li Average normalized @f$ \phi^2 @f$ energy contribution
 * @li Average normalized @f$ \chi^2 @f$ energy contribution
 * @li Average normalized @f$ \nabla \phi @f$ energy contribution
 * @li Average normalized @f$ \nabla \chi @f$ energy contribution
 * @li Average normalized potential-energy contribution
 * @li Average physical @f$ \phi'^2 @f$ energy contribution
 * @li Average physical @f$ \chi'^2 @f$ energy contribution
 * @li Average physical @f$ \phi\phi' @f$ energy contribution
 * @li Average physical @f$ \chi\chi' @f$ energy contribution
 * @li Average physical @f$ \phi^2 @f$ energy contribution
 * @li Average physical @f$ \chi^2 @f$ energy contribution
 * @li Average physical @f$ \nabla \phi @f$ energy contribution
 * @li Average physical @f$ \nabla \chi @f$ energy contribution
 * @li Average physical potential-energy contribution
 * @li Average physical @f$ \nabla \phi @f$ x-direction energy contribution
 * @li Average physical @f$ \nabla \chi @f$ x-direction energy contribution
 * @li Average physical @f$ \nabla \phi @f$ y-direction energy contribution
 * @li Average physical @f$ \nabla \chi @f$ y-direction energy contribution
 * @li Average physical @f$ \nabla \phi @f$ z-direction energy contribution
 * @li Average physical @f$ \nabla \chi @f$ z-direction energy contribution
 * @li Average physical pressure
 * @li Average w (the e.o.s. parameter)
 */

// Explicit instantiations
template class energy_outputter<double>;
#ifdef USE_LD
template class energy_outputter<long double>;
#endif
