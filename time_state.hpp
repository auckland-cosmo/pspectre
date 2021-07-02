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

template <typename K, typename V>
struct keyed_value
{
	keyed_value(V &v, K ik, V dv, const char *kn, const char *vn)
	: value(v), initial_key(ik), default_value(dv), key_name(kn), value_name(vn) {}

	V &value;
	const K initial_key;
	const V default_value;

	void advance(K k)
	{
		if (values.size() > 0) {
			if (values[0].first < k) {
				value = values[0].second;
				values.erase(values.begin());
				
				std::cout << value_name << " is now: " << value << std::endl;
			}
		}
	}

	void add_value(K start_key, V value_)
	{
		if (values.size() < 1) {
			if (start_key == initial_key) {
				value = value_;
			}
			else {
				add_value(initial_key, default_value);
			}
		}
		else {
			for (typename std::vector< std::pair<K, V> >::iterator it = values.begin(); it != values.end(); ++it) {
				if (it->first == start_key) {
					it->second = value_;
					break;
				}
			}

			// Ignore if the start key is less than one previously specified.
			if (start_key <= values[values.size()-1].first) {
				return;
			}
		}

		values.push_back(std::make_pair(start_key, value_));		
	}

	void finalize_values()
	{
		if (values.size() < 1) {
			add_value(initial_key, default_value);
		}
	}

	void summary(std::ostream& os)
	{
		for (typename std::vector< std::pair<K, V> >::iterator it = values.begin(); it != values.end(); ++it) {
			os << value_name << ": " << it->second << " (starting at " << key_name << " = " << it->first << ")" << std::endl;
		}
	}

protected:
	const char *key_name, *value_name;
	std::vector< std::pair<K, V> > values;
};

template <typename R>
struct time_state
{
	time_state()
	: t(0), physical_time(0), a(1), adot(0), addot(0), dt(0.0),
	dtm(dt, 0.0, default_dt(), "t", "dt") {}
	
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
		dtm.advance(t);
	}

	void add_dt(R start_time, R dt_)
	{
		dtm.add_value(start_time, dt_);
	}

	void finalize_dts()
	{
		dtm.finalize_values();
	}

	void dt_summary(std::ostream& os)
	{
		dtm.summary(os);
	}
	
protected:
	keyed_value<R, R> dtm;
};

#endif // TIME_STATE

