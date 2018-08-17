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
 * @brief DEFROST-style initial conditions.
 */

#ifndef DEFROST_STYLE_INITIALIZER_HPP
#define DEFROST_STYLE_INITIALIZER_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "field.hpp"
#include "initializer.hpp"

/** 
 * @brief DEFROST-style initial conditions.
 */

template <typename R>
class defrost_style_initializer : public initializer<R>
{
public:
	defrost_style_initializer(field_size &fs_, model_params<R> &mp_,
		field<R> &phi_, field<R> &phidot_, field<R> &chi_, field<R> &chidot_, R adot_)
		: fs(fs_), mp(mp_), phi(phi_), phidot(phidot_), chi(chi_), chidot(chidot_), adot(adot_) {}

public:	
	/** 
	 * @brief Initialize the phi, phidot, chi and chidot fields.
	 */

	virtual void initialize();

protected:
	/** 
	 * @brief Sample a Gaussian random field.
	 *
	 * Random Gaussian-mode amplitudes @f$b_k@f$ are chosen such that @f$<b_k b^*_{k'}> = \delta(k - k')@f$
	 * using the Box-Muller transformation. The kernel function is defined as:
	 * @f[
	 *	\zeta(r) = \frac{1}{\sqrt{\pi}} \int dk \, k^2 (k^2 + m_{\mbox{eff}})^\gamma \frac{\sin(kr)}{kr} e^{-k^2/q^2}
	 * @f]
	 *
	 * @f$q@f$ is chosen to be some scale below the Nyquist frequency.
	 * 
	 * @param fld The field into which to store the random field sample.
	 * @param gamma The @f$(k^2 + m^2)@f$ exponent.
	 * @param m2eff The effective mass.
	 */

	void sample_grf(field<R> &fld, R gamma, R m2eff);

protected:
	field_size &fs;
	model_params<R> &mp;
	field<R> &phi, &phidot;
	field<R> &chi, &chidot;
	R adot;
};

#endif // DEFROST_STYLE_INITIALIZER_HPP

