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
 * @brief Field grid size and derived size-related quantities.
 */

#ifndef FIELD_SIZE_HPP
#define FIELD_SIZE_HPP

#include <cmath>

struct field_size
{
	field_size(int n_ = 0, int n_pad_factor_ = 1)
		: n(n_), n_pad_factor(n_pad_factor_)
	{
		calculate_size_totals();
	}
	
	void calculate_size_totals()
	{
		using namespace std;
		
		total_gridpoints = n*n*n;
		total_padded_gridpoints = total_gridpoints*n_pad_factor*n_pad_factor*n_pad_factor;
		
		total_momentum_gridpoints = n*n*(n/2+1);
		total_padded_momentum_gridpoints = n_pad_factor*n*n_pad_factor*n*((n_pad_factor*n)/2+1);
		
		power_length = int(sqrt(3)*0.5*n) + 1;
	}

	int n, n_pad_factor;
	int total_gridpoints, total_padded_gridpoints;
	int total_momentum_gridpoints, total_padded_momentum_gridpoints;
	int power_length; // Length of a power spectrum array.
};

#endif // FIELD_SIZE_HPP
