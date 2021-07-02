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
 * @brief Fourth-order Runge–Kutta (RK4) integrator.
 */

#ifndef RK4_HPP
#define RK4_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "time_state.hpp"
#include "field.hpp"
#include "integrator.hpp"
#include "nonlinear_transformer.hpp"
#include "v_integrator.hpp"

template <typename R>
class rk4 : public integrator<R>
{
public:
	rk4(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
		field<R> &phi_, field<R> &phidot_, field<R> &chi_, field<R> &chidot_)
		: fs(fs_), upfs(fs_.n), mp(mp_), ts(ts_),
		phi(phi_), phidot(phidot_), chi(chi_), chidot(chidot_),
		phi1("phi1"), phidot1("phidot1"),
		chi1("chi1"), chidot1("chidot1"),
		phi2("phi2"), phidot2("phidot2"),
		chi2("chi2"), chidot2("chidot2"),
		phi3("phi3"), phidot3("phidot3"),
		chi3("chi3"), chidot3("chidot3"),
		nlt(fs_, mp_, ts_), vi(fs_, mp_),
		a1(0), a2(0), a3(0), a4(0), adot1(0), adot2(0), adot3(0), adot4(0),
		pt1(0), pt2(0), pt3(0), pt4(0)
		{
			phi1.construct(fs);
			phidot1.construct(upfs);

			phi2.construct(fs);
			phidot2.construct(upfs);

			phi3.construct(fs);
			phidot3.construct(upfs);

			chi1.construct(fs);
			chidot1.construct(upfs);

			chi2.construct(fs);
			chidot2.construct(upfs);

			chi3.construct(fs);
			chidot3.construct(upfs);
		}

public:	
	virtual void step();
	virtual void initialize();

protected:
	void substep_scale(R fac, field<R> &phip, field<R> &chip,
		R ap, R adotp, R ptp, R &an, R &adotn, R &ptn, R &dan, R &dadotn, R &dptn,
		R &avg_gradient_phi, R &avg_gradient_chi);

	void substep(R fac, field<R> &phip, field<R> &chip, field<R> &phidotp, field<R> &chidotp,
		field<R> &phin, field<R> &chin, field<R> &phidotn, field<R> &chidotn,
		R ap, R adotp, R ptp, R &an, R &adotn, R &ptn, R &avg_gradient_phi, R &avg_gradient_chi);

protected:
	field_size &fs, upfs;
	model_params<R> &mp;
	time_state<R> &ts;
	field<R> &phi, &phidot;
	field<R> &chi, &chidot;
	field<R> phi1, phidot1;
	field<R> chi1, chidot1;
	field<R> phi2, phidot2;
	field<R> chi2, chidot2;
	field<R> phi3, phidot3;
	field<R> chi3, chidot3;
	nonlinear_transformer<R> nlt;
	v_integrator<R> vi;
	R a1, a2, a3, a4; 
	R adot1, adot2, adot3, adot4;
	R pt1, pt2, pt3, pt4;
};

#endif // RK4_HPP

