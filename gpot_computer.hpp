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
 * @brief Gravitational-potential computations.
 */

#ifndef GPOT_COMPUTER_HPP
#define GPOT_COMPUTER_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "time_state.hpp"
#include "field.hpp"
#include "grad_computer.hpp"

/** 
 * @brief Computer of the gravitational potential from the energy density of the phi and chi fields.
 */

template <typename R>
class gpot_computer {
public:
	gpot_computer(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
		field<R> &phi_, field<R> &chi_, field<R> &phidot_, field<R> &chidot_, grad_computer<R> &gc_)
		: fs(fs_), upfs(fs_.n), mp(mp_), ts(ts_), phi(phi_), chi(chi_), phidot(phidot_), chidot(chidot_), gc(gc_),
		gpot("gpot")
	{
		gpot.construct(upfs);
	}

public:
	/** 
	 * @brief Compute gpot.
	 * 
	 * @param final_state The final state of gpot.
	 * @param grad_computed True if the gradient fields have already been computed (otherwise gc.compute() is called).
	 */

	void compute(field_state final_state = position, bool grad_computed = false);

protected:
	field_size &fs, upfs;
	model_params<R> &mp;
	time_state<R> &ts;
	field<R> &phi, &chi;
	field<R> &phidot, &chidot;
	grad_computer<R> &gc;

public:

	/** 
	 * @brief The gravitational potential field.
	 */

	field<R> gpot;
};

#endif // GPOT_COMPUTER_HPP
