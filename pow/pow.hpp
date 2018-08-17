//          Copyright (C) 2010 Hal Finkel.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

/**
 * @file
 * @brief Template function to compute the integer power of its argument.
 */

#ifndef POW_HPP
#define POW_HPP

// This is a modified version of the pow template from:
// C++ Meta<Programming> Concepts and Results
// Walter E. Brown, April 2001.

#ifndef DOXYGEN
template <unsigned N, typename R>
struct do_pow
{
	static inline R pow(R x)
	{
		return do_pow<N % 2, R>::pow(x) * do_pow<N/2, R>::pow(x*x);
	}
};

template <typename R>
struct do_pow<1u, R>
{
	static inline R pow(R x)
	{
		return x;
	}
};

template <typename R>
struct do_pow<0u, R>
{
	static inline R pow(R x)
	{
		return R(1.0);
	}
};
#endif // DOXYGEN

template <unsigned N, typename R>
static inline R pow(R x)
{
	return do_pow<N, R>::pow(x);
}

#endif // POW_HPP
