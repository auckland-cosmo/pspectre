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
 * @brief Field slice output manager.
 */

#ifndef SLICE_OUTPUT_MANAGER_HPP
#define SLICE_OUTPUT_MANAGER_HPP

#include "field_size.hpp"
#include "model_params.hpp"
#include "time_state.hpp"
#include "field.hpp"
#include "slice_outputter.hpp"
#include "grad_computer.hpp"
#include "gpot_computer.hpp"

#include <string>
#include <vector>

template <typename R>
class slice_output_manager
{
public:
	typedef typename slice_outputter<R>::var_func var_func;	

public:
	slice_output_manager(field_size &fs_, model_params<R> &mp_, time_state<R> &ts_,
		field<R> &phi_, field<R> &chi_, field<R> &phidot_, field<R> &chidot_,
		grad_computer<R> &gc_, gpot_computer<R> &gpotc_,
		int slicedim_ = 3, int slicelength_ = 0, int sliceskip_ = 1,
		bool sliceaverage_ = true, bool sliceflt_ = true)
		: fs(fs_), upfs(fs_.n), mp(mp_), ts(ts_), phi(phi_), chi(chi_),
		phidot(phidot_), chidot(chidot_), gc(gc_), gpotc(gpotc_),
		slicedim(slicedim_), slicelength(slicelength_ ? slicelength_ : fs_.n),
		sliceskip(sliceskip_), sliceaverage(sliceaverage_), sliceflt(sliceflt_), bin_idx(0) {}
	
public:
	void add_outputter(std::string varname, var_func vf);
	void output();

protected:
	field_size &fs, upfs;
	model_params<R> &mp;
	time_state<R> &ts;
	field<R> &phi, &chi;
	field<R> &phidot, &chidot;
	grad_computer<R> &gc;
	gpot_computer<R> &gpotc;
	int slicedim, slicelength, sliceskip;
	bool sliceaverage, sliceflt;
	int bin_idx;
	std::vector< slice_outputter<R> * > outputters;
};

#endif // SLICE_OUTPUT_MANAGER_HPP

