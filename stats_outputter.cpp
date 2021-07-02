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
#include "stats_outputter.hpp"

#include <sstream>
#include <iomanip>

using namespace std;

template <typename R>
stats_outputter<R>::stats_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
	field<R> &phi_, field<R> &chi_, field<R> &phidot_, field<R> &chidot_)
	: fs(fs_), mp(mp_), ts(ts_), phi(phi_), chi(chi_), phidot(phidot_), chidot(chidot_)
{
	of.open("stats.tsv");
	of << setprecision(30) << scientific;
}

template <typename R>
void stats_outputter<R>::compute(field<R> &fld1, field<R> &fld2, R &fld1_mean, R &fld2_mean,
	R &fld1_var, R &fld2_var)
{
	fld1.switch_state(momentum, true);
	fld2.switch_state(momentum, true);

	R total_fld1_squared = 0.0, total_fld2_squared = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_fld1_squared,total_fld2_squared)
#endif

	for (int x = 0; x < fs.n; ++x) {
		for (int y = 0; y < fs.n; ++y) {
			for (int z = 0; z < fs.n/2+1; ++z) {
				int idx = z + (fs.n/2 + 1) * (y + fs.n * x);
				int cnt = (z == 0 || z == fs.n/2) ? 1 : 2;

				total_fld1_squared += cnt * (pow<2>(fld1.mdata[idx][0]) + pow<2>(fld1.mdata[idx][1]));
				total_fld2_squared += cnt * (pow<2>(fld2.mdata[idx][0]) + pow<2>(fld2.mdata[idx][1]));
			}
		}
	}

	// The factor of 1./N^3 comes from Parseval's theorem.
	total_fld1_squared /= fs.total_gridpoints;
	total_fld2_squared /= fs.total_gridpoints;

	// Keep in mind that the zero-mode (DC mode) in momentum-space is also the average value of the field.
	fld1_mean = fld1.mdata[0][0]/fs.total_gridpoints;
	fld2_mean = fld2.mdata[0][0]/fs.total_gridpoints;

	// Compute the field variance using the formula Var(X) = <X^2> - <X>^2.
	fld1_var = total_fld1_squared/fs.total_gridpoints - pow<2>(fld1_mean);
	fld2_var = total_fld2_squared/fs.total_gridpoints - pow<2>(fld2_mean);
}

template <typename R>
void stats_outputter<R>::compute_cov(field<R> &fld, field<R> &flddot, R &fld_flddot_cov)
{
	// Uses the Plancherel theorem (the generalization of Parseval's theorem) to compute the covariance:
	// \sum_{i=0}^{N-1} x_i y_i^* = 1/N \sum_{k=0}^{N-1} X_k Y_k^*
	// Note that Cov(X, Y) = E[XY] - E[X]*E[Y];

	fld.switch_state(momentum, true);
	flddot.switch_state(momentum, true);

	// Because the answer *must* be real, keep only the real part...
	R total_ffd_mul = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_ffd_mul)
#endif

	for (int x = 0; x < fs.n; ++x) {
		for (int y = 0; y < fs.n; ++y) {
			for (int z = 0; z < fs.n/2+1; ++z) {
				int idx = z + (fs.n/2 + 1) * (y + fs.n * x);
				int cnt = (z == 0 || z == fs.n/2) ? 1 : 2;

				total_ffd_mul += cnt * (
					fld.mdata[idx][0]*flddot.mdata[idx][0] - fld.mdata[idx][1]*flddot.mdata[idx][1]
				);
			}
		}
	}

	// The factor of 1./N^3 comes from the Plancherel theorem.
	total_ffd_mul /= fs.total_gridpoints;

	R fld_mean = fld.mdata[0][0]/fs.total_gridpoints;
	R flddot_mean = flddot.mdata[0][0]/fs.total_gridpoints;
	fld_flddot_cov = total_ffd_mul/fs.total_gridpoints - fld_mean*flddot_mean;
}

template <typename R>
void stats_outputter<R>::output()
{
	R phi_avg, phi_var, chi_avg, chi_var, phidot_avg, phidot_var, chidot_avg, chidot_var;
	compute(phi, chi, phi_avg, chi_avg, phi_var, chi_var);
	compute(phidot, chidot, phidot_avg, chidot_avg, phidot_var, chidot_var);

	R phi_phidot_cov, chi_chidot_cov;
	compute_cov(phi, phidot, phi_phidot_cov);
	compute_cov(chi, chidot, chi_chidot_cov);

	R tophys = mp.rescale_A * pow(ts.a, mp.rescale_r);

	// Note: f'/B = a^{s-r} f'_pr/A - r a^{-1-r} a'/B f_pr/A
	R phidot_avg_pt = pow(ts.a, mp.rescale_s-mp.rescale_r)*phidot_avg/mp.rescale_A
		- mp.rescale_r*pow(ts.a, -mp.rescale_r-1)*ts.adot*phi_avg/mp.rescale_A;
	R chidot_avg_pt = pow(ts.a, mp.rescale_s-mp.rescale_r)*chidot_avg/mp.rescale_A
		- mp.rescale_r*pow(ts.a, -mp.rescale_r-1)*ts.adot*chi_avg/mp.rescale_A;

	// Note: Var(aX + bY) = a^2 Var(X) + b^2 Var(Y) + 2ab Cov(X, Y)
	R var_c1 = pow(ts.a, mp.rescale_s-mp.rescale_r)/mp.rescale_A;
	R var_c2 = -mp.rescale_r*ts.adot*pow(ts.a, -mp.rescale_r-1)/mp.rescale_A;
	R phidot_var_pt = pow<2>(var_c1)*phidot_var + pow<2>(var_c2)
		+ 2*var_c1*var_c2*phi_phidot_cov;
	R chidot_var_pt = pow<2>(var_c1)*chidot_var + pow<2>(var_c2)
		+ 2*var_c1*var_c2*chi_chidot_cov;

	of << setw(10) << setfill('0') <<
		ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<

		phi_avg << "\t" << phi_var << "\t" <<
		chi_avg << "\t" << chi_var << "\t" <<
		phi_avg/tophys << "\t" << phi_var/pow<2>(tophys) << "\t" <<
		chi_avg/tophys << "\t" << chi_var/pow<2>(tophys) << "\t" <<

		phidot_avg << "\t" << phidot_var << "\t" <<
		chidot_avg << "\t" << chidot_var << "\t" <<

		phidot_avg_pt << "\t" <<
		phidot_var_pt << "\t" <<
		chidot_avg_pt << "\t" <<
		chidot_var_pt << "\t" <<

		mp.rescale_B*phidot_avg_pt << "\t" <<
		pow<2>(mp.rescale_B)*phidot_var_pt << "\t" << 
		mp.rescale_B*chidot_avg_pt << "\t" <<
		pow<2>(mp.rescale_B)*chidot_var_pt << "\t" << 

		endl;
	of.flush(); 
}

/**
 * @page stats_tsv stats.tsv
 * stats.tsv is a tab serarated file with the following fields:
 * @li Program time
 * @li Physical time
 * @li Mean of phi in program units
 * @li Variance of phi in program units
 * @li Mean of chi in program units
 * @li Variance of chi in program units
 * @li Mean of phi
 * @li Variance of phi
 * @li Mean of chi
 * @li Variance of chi
 * @li Mean of phidot in program units
 * @li Variance of phidot in program units
 * @li Mean of chidot in program units
 * @li Variance of chidot in program units
 * @li Mean of phidot (rescaled time)
 * @li Variance of phidot (rescaled time)
 * @li Mean of chidot (rescaled time)
 * @li Variance of chidot (rescaled time)
 * @li Mean of phidot
 * @li Variance of phidot
 * @li Mean of chidot
 * @li Variance of chidot
 */

// Explicit instantiations
template class stats_outputter<double>;
#ifdef USE_LD
template class stats_outputter<long double>;
#endif
