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

#include "slice_outputter.hpp"

#include <sstream>
#include <iomanip>

using namespace std;

template <typename R>
slice_outputter<R>::slice_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
	int slicelength_, std::string varname_, var_func vf_, bool flt_)
	: fs(fs_), mp(mp_), ts(ts_), slicelength(slicelength_), varname(varname_),
	vf(vf_), cp(0), cn(0), flt(flt_)
{
	buffer = new R[slicelength];
	if (flt) bufferf = new float[slicelength];
	
	cout << "Slice output enabled for: " << varname << endl;
}

template <typename R>
slice_outputter<R>::~slice_outputter()
{
	delete [] buffer;
	if (flt) delete [] bufferf;
}

template <typename R>
void slice_outputter<R>::begin(int bin_idx)
{
	stringstream ss;
	ss << varname << "_" <<
		setfill('0') << setw(5) << bin_idx << ".bin";
	
	if (of.is_open()) {
		of.close();	
	}
	
	of.open(ss.str().c_str(), ios::binary);
	
	cp = cn = 0;
	buffer[cp] = 0.0;
}

template <typename R>
void slice_outputter<R>::flush()
{
	if (flt) {
		for (int i = 0; i < cp; ++i) {
			bufferf[i] = (float) buffer[i];
		}
		
		of.write((char *) bufferf, cp*sizeof(float));
	}
	else {
		of.write((char *) buffer, cp*sizeof(R));
	}

	cp = cn = 0;
	buffer[cp] = 0.0;
}

template <typename R>
void slice_outputter<R>::advance()
{
	buffer[cp++] /= cn;
	buffer[cp] = 0.0;
	cn = 0;
}

template <typename R>
void slice_outputter<R>::accumulate(R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
	R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	buffer[cp] += vf(fs, mp, ts, phi, chi, phidot, chidot, phigradx, chigradx,
		phigrady, chigrady, phigradz, chigradz, gpot);
	++cn;
}

// Explicit instantiations
template class slice_outputter<double>;
#ifdef USE_LD
template class slice_outputter<long double>;
#endif
