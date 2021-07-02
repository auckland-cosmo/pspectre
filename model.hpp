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
 * @brief A particular simulated situation.
 */

#ifndef MODEL_HPP
#define MODEL_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "time_state.hpp"
#include "field.hpp"
#include "integrator.hpp"
#include "slice_output_manager.hpp"
#include "grad_computer.hpp"
#include "gpot_computer.hpp"

#ifdef HAVE_PRIVATE
#include "private/private_globals.hpp"
#endif

template <typename R>
class model
{
public:
	model(int argc, char *argv[]);
	~model();

public:
	void run();
	
protected:
	void set_output_directory(const char *uodn);
	void write_info_file();
	void set_initial_conditions();
	void evolve(integrator<R> *ig);
	void load_initial_slice_file(std::string &ifn, field<R> &fld, R pf);

protected:
	void private_allocate();
	void private_set_sf_info();
	void private_evolve(int counter);
	void private_info_file_output(std::ofstream &info_file);

protected:
	field_size fs;
	model_params<R> mp;
	time_state<R> ts;
	bool use_verlet, le_init, homo_ic_phi, homo_ic_chi;
	int seed;
	R tf;
	int scale_interval, energy_interval, spectra_interval,
		screen_interval, slice_interval, stats_interval,
		twoptcorr_interval;
	keyed_value<R, int> scale_intervals, energy_intervals, spectra_intervals,
                screen_intervals, slice_intervals, stats_intervals,
                twoptcorr_intervals;

protected:
	field<R> phi, phidot;
	field<R> chi, chidot;
	grad_computer<R> *gc;
	gpot_computer<R> *gpotc;
	slice_output_manager<R> *som;
	R ics_scale;
	R len0;
	bool vvwl;
	R af;
	bool external_H0;
	std::string phi0_slice, chi0_slice;
	std::string phidot0_slice, chidot0_slice;
	std::string start_wd;
	int ics_eff_size;
	R phidot0pr, chidot0pr;

#ifdef HAVE_PRIVATE
	private_globals<R> priv;
#endif
};

#endif // MODEL_HPP
