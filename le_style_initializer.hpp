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

/**
 * @file
 * @brief LatticeEasy-style initialization.
 */

#ifndef LE_STYLE_INITIALIZER_HPP
#define LE_STYLE_INITIALIZER_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "field.hpp"
#include "initializer.hpp"

#include <cmath>

template <typename R>
class le_style_initializer : public initializer<R>
{
public:
	le_style_initializer(field_size &fs_, model_params<R> &mp_,
		field<R> &phi_, field<R> &phidot_, field<R> &chi_, field<R> &chidot_, R adot_, R len0)
		: fs(fs_), mp(mp_), phi(phi_), phidot(phidot_), chi(chi_), chidot(chidot_), adot(adot_)
	{
		using namespace std;
		
		fluctuation_amplitude = mp_.rescale_A * mp_.rescale_B * fs_.total_gridpoints /
				(pow(mp_.len/len0, (R) 1.5) * sqrt(2.));
	}

public:
	virtual void initialize();

protected:
	void set_mode(field<R> &fld, field<R> &flddot, R m_fld_eff, int px, int py, int pz, int idx, bool real = false);
	void initialize_field(field<R> &fld, field<R> &flddot, R m_fld_eff);

protected:
	field_size &fs;
	model_params<R> &mp;
	field<R> &phi, &phidot;
	field<R> &chi, &chidot;
	R adot;
	R fluctuation_amplitude;
};

#endif // LE_STYLE_INITIALIZER_HPP
