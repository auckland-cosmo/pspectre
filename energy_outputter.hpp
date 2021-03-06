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
 * @brief Outputter for the energy TSV file.
 */

#ifndef ENERGY_OUTPUTTER_HPP
#define ENERGY_OUTPUTTER_HPP

#include "field.hpp"
#include "model_params.hpp"
#include "time_state.hpp"
#include "v_integrator.hpp"

#include <iostream>
#include <fstream>

/** 
 * @brief Outputter for the energy TSV file.
 */

template <typename R>
class energy_outputter {
public:
	energy_outputter(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
		field<R> &phi_, field<R> &chi_, field<R> &phidot_, field<R> &chidot_);

public:
	/** 
	 * @brief Compute the integrated energy density and write it to the file.
	 * 
	 * @param no_output If true the result is not output to the file.
	 */

	void output(bool no_output = false);

protected:
	field_size &fs;
	model_params<R> &mp;
	time_state<R> &ts;
	field<R> &phi, &chi;
	field<R> &phidot, &chidot;
	v_integrator<R> vi;
	std::ofstream of;

public:
	/** 
	 * @brief The average energy density in physical units.
	 */

	R avg_rho_phys;

	/** 
	 * @brief The average energy density normalized by the Friedmann equation.
	 */

	R avg_rho;
};

#endif // ENERGY_OUTPUTTER_HPP

