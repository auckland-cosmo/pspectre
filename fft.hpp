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
 * @brief FFT wrappers.
 */

#ifndef FFT_HPP
#define FFT_HPP

#ifdef USE_MKL
#include <fftw/fftw3.h>

#ifdef HAS_FFTW3_MKL_H
#include <fftw/fftw3_mkl.h>
#endif

#ifdef USE_LD
#error MKL does not support long double precision.
#endif

#else
#include <fftw3.h>
#endif

template <typename R>
inline R* fft_malloc(size_t sz)
{
	return 0;
}

template <>
inline double *fft_malloc<double>(size_t sz)
{
	return (double *) fftw_malloc(sz);
}

#ifdef USE_LD
template <>
inline long double *fft_malloc<long double>(size_t sz)
{
	return (long double *) fftwl_malloc(sz);
}
#endif

template <typename R>
inline void fft_free(R *ptr) {}

template <>
inline void fft_free<double>(double *ptr)
{
	return fftw_free(ptr);
}

#ifdef USE_LD
template <>
inline void fft_free<long double>(long double *ptr)
{
	return fftwl_free(ptr);
}
#endif

enum fft_r2r_kind
{
	r2hc = FFTW_R2HC,
	hc2r = FFTW_HC2R,
	dht = FFTW_DHT,
	redft00 = FFTW_REDFT00,
	redft10 = FFTW_REDFT10,
	redft01 = FFTW_REDFT01,
	redft11 = FFTW_REDFT11,
	rodft00 = FFTW_RODFT00,
	rodft10 = FFTW_RODFT10,
	rodft01 = FFTW_RODFT01,
	rodft11 = FFTW_RODFT11
};

template <typename R>
class fft_r2r_1d_plan {};

template <>
class fft_r2r_1d_plan<double>
{
public:
	fft_r2r_1d_plan(int n, double *in, double *out, fft_r2r_kind kind, bool estimate = true)
	{
		construct(n, in, out, kind, estimate);
	}
	
	fft_r2r_1d_plan()
		: plan(0) {}
	
	~fft_r2r_1d_plan()
	{
		fftw_destroy_plan(plan);
	}
	
public:
	void construct(int n, double *in, double *out, fft_r2r_kind kind, bool estimate = true)
	{
		plan = fftw_plan_r2r_1d(n, in, out, (fftw_r2r_kind) kind, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}	

	void execute()
	{
		fftw_execute(plan);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftw_plan plan;
};

#ifdef USE_LD
template <>
class fft_r2r_1d_plan<long double>
{
public:
	fft_r2r_1d_plan(int n, long double *in, long double *out, fft_r2r_kind kind, bool estimate = true)
	{
		construct(n, in, out, kind, estimate);
	}
	
	fft_r2r_1d_plan()
		: plan(0) {}
	
	~fft_r2r_1d_plan()
	{
		fftwl_destroy_plan(plan);
	}
	
public:
	void construct(int n, long double *in, long double *out, fft_r2r_kind kind, bool estimate = true)
	{
		plan = fftwl_plan_r2r_1d(n, in, out, (fftwl_r2r_kind) kind, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}	

	void execute()
	{
		fftwl_execute(plan);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftwl_plan plan;
};
#endif

template <typename R>
class fft_dft_c2r_3d_plan {};

template <>
class fft_dft_c2r_3d_plan<double>
{
public:
	typedef fftw_complex complex_t;
	
public:
	fft_dft_c2r_3d_plan(int n0, int n1, int n2, complex_t *in, double *out, bool estimate = true)
	{
		construct(n0, n1, n2, in, out, estimate);
	}
	
	fft_dft_c2r_3d_plan()
		: plan(0) {}
	
	~fft_dft_c2r_3d_plan()
	{
		fftw_destroy_plan(plan);
	}
	
public:
	void construct(int n0, int n1, int n2, complex_t *in, double *out, bool estimate = true)
	{
		plan = fftw_plan_dft_c2r_3d(n0, n1, n2, in, out, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}	

	void execute()
	{
		fftw_execute(plan);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftw_plan plan;
};

#ifdef USE_LD
template <>
class fft_dft_c2r_3d_plan<long double>
{
public:
	typedef fftwl_complex complex_t;

public:
	fft_dft_c2r_3d_plan(int n0, int n1, int n2, complex_t *in, long double *out, bool estimate = true)
	{
		construct(n0, n1, n2, in, out, estimate);
	}
	
	fft_dft_c2r_3d_plan()
		: plan(0) {}

	~fft_dft_c2r_3d_plan()
	{
		fftwl_destroy_plan(plan);
	}

public:
	void construct(int n0, int n1, int n2, complex_t *in, long double *out, bool estimate = true)
	{
		plan = fftwl_plan_dft_c2r_3d(n0, n1, n2, in, out, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}

	void execute()
	{
		fftwl_execute(plan);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftwl_plan plan;
};
#endif

template <typename R>
class fft_dft_r2c_3d_plan {};

template <>
class fft_dft_r2c_3d_plan<double>
{
public:
	typedef fftw_complex complex_t;
	
public:
	fft_dft_r2c_3d_plan(int n0, int n1, int n2, double *in, complex_t *out, bool estimate = true)
	{
		construct(n0, n1, n2, in, out, estimate);
	}
	
	fft_dft_r2c_3d_plan()
		: plan(0) {}
	
	~fft_dft_r2c_3d_plan()
	{
		fftw_destroy_plan(plan);
	}

public:
	void execute()
	{
		fftw_execute(plan);
	}
	
	void construct(int n0, int n1, int n2, double *in, complex_t *out, bool estimate = true)
	{
		plan = fftw_plan_dft_r2c_3d(n0, n1, n2, in, out, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftw_plan plan;
};

#ifdef USE_LD
template <>
class fft_dft_r2c_3d_plan<long double>
{
public:
	typedef fftwl_complex complex_t;

public:
	fft_dft_r2c_3d_plan(int n0, int n1, int n2, long double *in, complex_t *out, bool estimate = true)
	{
		construct(n0, n1, n2, in, out, estimate);
	}
	
	fft_dft_r2c_3d_plan()
		: plan(0) {}

	~fft_dft_r2c_3d_plan()
	{
		fftwl_destroy_plan(plan);
	}

public:
	void execute()
	{
		fftwl_execute(plan);
	}
	
	void construct(int n0, int n1, int n2, long double *in, complex_t *out, bool estimate = true)
	{
		plan = fftwl_plan_dft_r2c_3d(n0, n1, n2, in, out, estimate ? FFTW_ESTIMATE : FFTW_MEASURE);
	}

	bool constructed() {
		return plan == 0;
	}

protected:
	fftwl_plan plan;
};
#endif

#endif // FFT_HPP
