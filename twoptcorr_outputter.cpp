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
#include "twoptcorr_outputter.hpp"

#include <sstream>
#include <iomanip>

using namespace std;

template <typename R>
twoptcorr_outputter<R>::twoptcorr_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
	field<R> &phi_, field<R> &chi_)
	: fs(fs_), upfs(fs_.n), mp(mp_), ts(ts_), phi(phi_), chi(chi_), corr("corr")
{
	of.open("twoptcorr.tsv");
	of << setprecision(30) << scientific;

	corr.construct(upfs);

	dmax = int(sqrt(3.0)*(fs.n-1));
	
	phi_total = new R[dmax];
	chi_total = new R[dmax];
	
	counts = new int[dmax];
}

template <typename R>
twoptcorr_outputter<R>::~twoptcorr_outputter()
{
	delete [] phi_total;
	delete [] chi_total;
	
	delete [] counts;
}

template <typename R>
void twoptcorr_outputter<R>::output()
{
	for (int i = 0; i < dmax; ++i) {
		phi_total[i] = chi_total[i] = 0.0;
		counts[i] = 0;
	}
	
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	corr.switch_state(uninitialized);
	corr.switch_state(momentum);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for (int x = 0; x < fs.n; ++x)
	for (int y = 0; y < fs.n; ++y)
	for (int z = 0; z < fs.n/2+1; ++z) {
		int idx = z + (fs.n/2 + 1) * (y + fs.n * x);
		corr.mdata[idx][0] = pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]);
		corr.mdata[idx][1] = 0.0;
	}

	corr.mdata[0][0] = 0.0;

	corr.switch_state(position);

	for (int x = 0; x < fs.n; ++x)
	for (int y = 0; y < fs.n; ++y)
	for (int z = 0; z < fs.n; ++z) {
		int idx = z + corr.ldl*(y + fs.n*x);
		int d = int(sqrt(pow<2, R>(x) + pow<2, R>(y) + pow<2, R>(z)));
		if (d > 0) phi_total[d-1] += corr.data[idx];
		if (d > 0) counts[d-1] += 1;
	}

	corr.switch_state(uninitialized);
	corr.switch_state(momentum);

#ifdef _OPENMP
#pragma omp parallel for
#endif
	
	for (int x = 0; x < fs.n; ++x)
	for (int y = 0; y < fs.n; ++y)
	for (int z = 0; z < fs.n/2+1; ++z) {
		int idx = z + (fs.n/2 + 1) * (y + fs.n * x);
		corr.mdata[idx][0] = pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]);
		corr.mdata[idx][1] = 0.0;
	}

	corr.mdata[0][0] = 0.0;

	corr.switch_state(position);

	for (int x = 0; x < fs.n; ++x)
	for (int y = 0; y < fs.n; ++y)
	for (int z = 0; z < fs.n; ++z) {
		int idx = z + corr.ldl*(y + fs.n*x);
		int d = int(sqrt(pow<2, R>(x) + pow<2, R>(y) + pow<2, R>(z)));
		if (d > 0) chi_total[d-1] += corr.data[idx];
	}

	R dx = mp.len/fs.n;
	R pdx = dx/mp.rescale_B;
	for (int i = 0; i < dmax; ++i) {
		phi_total[i] /= counts[i]*pow<2, R>(fs.total_gridpoints);
		chi_total[i] /= counts[i]*pow<2, R>(fs.total_gridpoints);
		
		R tophys = pow<2>(mp.rescale_A * pow(ts.a, mp.rescale_r));
		R phi_total_phys = phi_total[i]/tophys;
		R chi_total_phys = chi_total[i]/tophys;
		
		of << setw(10) << setfill('0') << ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<
			i << "\t" << counts[i] << "\t" << dx*(i+1) << "\t" << pdx*(i+1) << "\t" <<
			phi_total[i] << "\t" << chi_total[i] << "\t" <<
			phi_total_phys << "\t" << chi_total_phys << endl;
	}
	
	of << endl << endl;
	of.flush(); 
}

/**
 * @page twoptcorr_tsv twoptcorr.tsv
 * twoptcorr.tsv is a tab serarated file with the following fields:
 * @li Program time
 * @li Physical time
 * @li Length-bin number
 * @li Length grid points
 * @li Program-unit length
 * @li Physical length
 * @li Program-unit phi two-point correlation
 * @li Program-unit chi two-point correlation
 * @li Physical phi two-point correlation
 * @li Physical chi two-point correlation
 */

// Explicit instantiations
template class twoptcorr_outputter<double>;
#ifdef USE_LD
template class twoptcorr_outputter<long double>;
#endif
