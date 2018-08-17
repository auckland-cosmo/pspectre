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
#include "spectra_outputter.hpp"

#include <sstream>
#include <iomanip>

using namespace std;

template <typename R>
spectra_outputter<R>::spectra_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
	field<R> &phi_, field<R> &chi_)
	: fs(fs_), mp(mp_), ts(ts_), phi(phi_), chi(chi_)
{
	of.open("spectra.tsv");
	of << setprecision(30) << scientific;
	
	phi_total = new R[fs.power_length];
	chi_total = new R[fs.power_length];
	
	counts = new int[fs.power_length];
}

template <typename R>
spectra_outputter<R>::~spectra_outputter()
{
	delete [] phi_total;
	delete [] chi_total;
	
	delete [] counts;
}

template <typename R>
void spectra_outputter<R>::output()
{
	for (int i = 0; i < fs.power_length; ++i) {
		phi_total[i] = chi_total[i] = 0.0;
		counts[i] = 0;
	}
	
	phi.switch_state(momentum, true);
	chi.switch_state(momentum, true);

	for (int x = 0; x < fs.n; ++x) {
		int px = (x <= fs.n/2 ? x : x - fs.n);
		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);
			for (int z = 0; z < fs.n/2+1; ++z) {
				int pz = z;
				int idx = z + (fs.n/2 + 1) * (y + fs.n * x);

				int cnt = (z == 0 || z == fs.n/2) ? 1 : 2;
				int i = (int) sqrt(pow<2>(px) + pow<2>(py) + pow<2>(pz));

				phi_total[i] += cnt * (pow<2>(phi.mdata[idx][0]) + pow<2>(phi.mdata[idx][1]));
				chi_total[i] += cnt * (pow<2>(chi.mdata[idx][0]) + pow<2>(chi.mdata[idx][1]));

				counts[i] += cnt;
			}
		}
	}

	for (int i = 0; i < fs.power_length; ++i) {
		phi_total[i] /= counts[i]*pow<2, R>(fs.total_gridpoints);
		chi_total[i] /= counts[i]*pow<2, R>(fs.total_gridpoints);
		
		R tophys = pow<2>(mp.rescale_A * pow(ts.a, mp.rescale_r));
		R phi_total_phys = phi_total[i]/tophys;
		R chi_total_phys = chi_total[i]/tophys;
		
		of << setw(10) << setfill('0') << ts.t << "\t" << mp.rescale_B * ts.physical_time << "\t" <<
			i << "\t" << counts[i] << "\t" << mp.dp*i << "\t" <<
			phi_total[i] << "\t" << chi_total[i] << "\t" <<
			phi_total_phys << "\t" << chi_total_phys << endl;
	}
	
	of << endl << endl;
	of.flush(); 
}

/**
 * @page spectra_tsv spectra.tsv
 * spectra.tsv is a tab serarated file with the following fields:
 * @li Program time
 * @li Physical time
 * @li Bin number
 * @li Grid points per bin
 * @li Physical bin momentum
 * @li Total program-unit phi momentum
 * @li Total program-unit chi momentum
 * @li Total physical phi momentum
 * @li Total physical chi momentum
 */

// Explicit instantiations
template class spectra_outputter<double>;
#ifdef USE_LD
template class spectra_outputter<long double>;
#endif
