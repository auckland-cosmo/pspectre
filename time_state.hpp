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
 * @brief Time-varying model parameters.
 */

#ifndef TIME_STATE
#define TIME_STATE

#include <vector>
#include <utility>

#include <iostream>

template <typename R>
struct time_state
{
	time_state()
	: t(0), physical_time(0), a(1), adot(0), addot(0), dt(0.0) {}
	
	R t, physical_time;
	R a, adot, addot;
	
	R dt;

	static R default_dt()
	{
		return R(0.005);
	}

	void advance()
	{
		t += dt;
		
		if (dts.size() > 0) {
			if (dts[0].first < t) {
				dt = dts[0].second;
				dts.erase(dts.begin());
				
				std::cout << "dt is now: " << dt << std::endl;
			}
		}
	}

	void add_dt(R start_time, R dt_)
	{
		if (dts.size() < 1) {
			if (start_time == 0.0) {
				dt = dt_;
			}
			else {
				add_dt(0.0, default_dt());
			}
		}
		else {
			for (typename std::vector< std::pair<R, R> >::iterator it = dts.begin(); it != dts.end(); ++it) {
				if (it->first == start_time) {
					it->second = dt_;
					break;
				}
			}

			// Ignore if the start time is less than one previously specified.
			if (start_time <= dts[dts.size()-1].first) {
				return;
			}
		}

		dts.push_back(std::make_pair(start_time, dt_));		
	}

	void finalize_dts()
	{
		if (dts.size() < 1) {
			add_dt(0.0, default_dt());
		}
	}

	void dt_summary(std::ostream& os)
	{
		for (typename std::vector< std::pair<R, R> >::iterator it = dts.begin(); it != dts.end(); ++it) {
			os << "dt: " << it->second << " (starting at t = " << it->first << ")" << std::endl;
		}
	}
	
protected:
	std::vector< std::pair<R, R> > dts;
};

#endif // TIME_STATE

