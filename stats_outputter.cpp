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
	field<R> &phi_, field<R> &chi_)
	: fs(fs_), mp(mp_), ts(ts_), phi(phi_), chi(chi_)
{
	of.open("stats.tsv");
	of << setprecision(30) << scientific;
}

template <typename R>
void stats_outputter<R>::output()
{
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	R total_phi_squared = 0.0, total_chi_squared = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+:total_phi_squared,total_chi_squared)
#endif

	for (int x = 0; x < fs.n; ++x) {
		for (int y = 0; y < fs.n; ++y) {
			for (int z = 0; z < fs.n/2+1; ++z) {
				int idx = z + (fs.n/2 + 1) * (y + fs.n * x);
				int cnt = (z == 0 || z == fs.n/2) ? 1 : 2;

				total_phi_squared += cnt * (pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]));
				total_chi_squared += cnt * (pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]));
			}
		}
	}

	// The factor of 1./N^3 comes from Parseval's theorem.
	total_phi_squared /= fs.total_gridpoints;
	total_chi_squared /= fs.total_gridpoints;

	// Compute the field variance using the formula Var(X) = <X^2> - <X>^2.
	// Keep in mind that the zero-mode (DC mode) in momentum-space is also the average value of the field.
	R phi_variance = total_phi_squared/fs.total_gridpoints - pow<2>(phi.mdata[0][0]/fs.total_gridpoints);
	R chi_variance = total_chi_squared/fs.total_gridpoints - pow<2>(chi.mdata[0][0]/fs.total_gridpoints);

	R tophys = mp.rescale_A * pow(ts.a, mp.rescale_r);

	of << setw(10) << setfill('0') <<
		ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" << phi.mdata[0][0]/fs.total_gridpoints << "\t" <<
		phi_variance << "\t" << chi.mdata[0][0]/fs.total_gridpoints << "\t" << chi_variance
		<< "\t" << (phi.mdata[0][0]/fs.total_gridpoints)/tophys << "\t" <<
		phi_variance/pow<2>(tophys) << "\t" << (chi.mdata[0][0]/fs.total_gridpoints)/tophys << "\t" << chi_variance/pow<2>(tophys)
		<< endl;
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
 */

// Explicit instantiations
template class stats_outputter<double>;
#ifdef USE_LD
template class stats_outputter<long double>;
#endif
