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
 * @brief Computation of the gradient in Fourier space.
 */

#ifndef GRAD_COMPUTER_HPP
#define GRAD_COMPUTER_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "field.hpp"

template <typename R>
class grad_computer {
public:
	grad_computer(field_size &fs_, model_params<R> &mp_, field<R> &phi_, field<R> &chi_)
		: fs(fs_), upfs(fs_.n), mp(mp_), phi(phi_), chi(chi_),
		phigradx("phigradx"), chigradx("chigradx"),
		phigrady("phigrady"), chigrady("chigrady"),
		phigradz("phigradz"), chigradz("chigradz")
	{
		phigradx.construct(upfs);
		chigradx.construct(upfs);
		
		phigrady.construct(upfs);
		chigrady.construct(upfs);
		
		phigradz.construct(upfs);
		chigradz.construct(upfs);		
	}

public:
	void compute(field_state final_state = position);

protected:
	field_size &fs, upfs;
	model_params<R> &mp;
	field<R> &phi, &chi;
	
public:
	field<R> phigradx, chigradx;
	field<R> phigrady, chigrady;
	field<R> phigradz, chigradz;	
};

#endif // GRAD_COMPUTER_HPP

