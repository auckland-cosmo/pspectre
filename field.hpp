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
 * @brief Three-dimensional scalar fields.
 */

#ifndef FIELD_HPP
#define FIELD_HPP

#include "fft.hpp"
#include "field_size.hpp"

#include <cmath>

enum field_state
{
	uninitialized,
	position,
	momentum,
	padded_position,
	padded_momentum
};

/** 
 * @brief A three-dimensional scalar field in both position and momentum space.
 */

template <typename R>
class field
{
public:
	typedef typename fft_dft_r2c_3d_plan<R>::complex_t complex_t;

public:
	field(field_size &fs_, bool oop = false, const char *name_ = 0)
		: state(uninitialized), name(name_)
	{		
		construct(fs_, oop);
	}
	
	field(const char *name_ = 0)
		: data(0), ldl(0), pldl(0), mdata(0), state(uninitialized), mdata_saved(0), name(name_) {};

	~field();

public:
	void construct(field_size &fs_, bool oop = false);
	void divby(R v);
	void switch_state(field_state state_, bool mmo = false);

public:
	bool is_in_place()
	{
		return (data == ((R *) mdata));
	}

protected:
	void pad_momentum_grid();
	void unpad_momentum_grid();

public:
	field_size fs;

	/** 
	 * @brief The position-space data.
	 *
	 * @note The inner (z) dimension is padded to a size of 2*(fs.n/2+1).
	 */

	R *data;

	/**
	 * @brief The length of the last dimension of the data array.
	 */

	int ldl;

	/**
	 * @brief The length of the last dimension of the padded data array.
	 */

	int pldl;

	/** 
	 * @brief The momentum-space data.
	 */

	typename fft_dft_c2r_3d_plan<R>::complex_t *mdata;
	
protected:
	field_state state;
	fft_dft_r2c_3d_plan<R> p2m_plan;
	fft_dft_c2r_3d_plan<R> m2p_plan;
	fft_dft_r2c_3d_plan<R> padded_p2m_plan;
	fft_dft_c2r_3d_plan<R> padded_m2p_plan;
	typename fft_dft_c2r_3d_plan<R>::complex_t *mdata_saved;

public:
	const char *name;
};

#endif // FIELD_HPP
