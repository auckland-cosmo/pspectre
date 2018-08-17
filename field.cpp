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

#define _XOPEN_SOURCE 600

#include "pow/pow.hpp"
#include "field.hpp"

#include <cstdlib>
#include <cstring>
#include <cmath>

#include <iostream>

using namespace std;

template <typename R>
void field<R>::construct(field_size &fs_, bool oop)
{
	fs = fs_;
	ldl = oop ? fs.n : 2*(fs.n/2+1);
	pldl = oop ? fs.n_pad_factor*fs.n : 2*((fs.n_pad_factor*fs.n)/2+1);

	size_t alloc_size = sizeof(R) * fs.total_padded_gridpoints;
	size_t c_alloc_size = 2*sizeof(R) * fs.total_padded_momentum_gridpoints;

	mdata = (typename fft_dft_c2r_3d_plan<R>::complex_t *) fft_malloc<R>(c_alloc_size);
	data = oop ? fft_malloc<R>(alloc_size) : (R *) mdata;

	// Out-of-place c2r implies FFTW_DESTROY_INPUT, data is saved to this array.
	mdata_saved = oop ? (typename fft_dft_c2r_3d_plan<R>::complex_t *) fft_malloc<R>(c_alloc_size) : 0;

	m2p_plan.construct(fs.n, fs.n, fs.n, mdata_saved ? mdata_saved : mdata, data, false);		
	p2m_plan.construct(fs.n, fs.n, fs.n, data, mdata, false);

	padded_m2p_plan.construct(fs.n_pad_factor*fs.n, fs.n_pad_factor*fs.n, fs.n_pad_factor*fs.n,
		mdata_saved ? mdata_saved : mdata, data, false);		
	padded_p2m_plan.construct(fs.n_pad_factor*fs.n, fs.n_pad_factor*fs.n, fs.n_pad_factor*fs.n,
		data, mdata, false);

	if (oop) memset(data, 0, alloc_size);
	memset(mdata, 0, c_alloc_size);
}

template <typename R>
field<R>::~field()
{
	if (!is_in_place()) {
		fft_free<R>((R *) mdata);
		fft_free<R>((R *) mdata_saved);
	}

	fft_free<R>(data);
}

template <typename R>
void field<R>::divby(R v)
{
	if (state == momentum) {

#ifdef _OPENMP
#pragma omp parallel for
#endif

		for (int idx = 0; idx < fs.total_momentum_gridpoints; ++idx) {
			mdata[idx][0] /= v;
			mdata[idx][1] /= v;
		}
	}
	else if (state == position) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int x = 0; x < fs.n; ++x)
		for (int y = 0; y < fs.n; ++y)
		for (int z = 0; z < fs.n; ++z) {
			int idx = z + ldl*(y + fs.n*x);
			data[idx] /= v;
		}
	}
	else if (state == padded_momentum) {

#ifdef _OPENMP
#pragma omp parallel for
#endif

		for (int idx = 0; idx < fs.total_padded_momentum_gridpoints; ++idx) {
			mdata[idx][0] /= v;
			mdata[idx][1] /= v;
		}		
	}
	else if (state == padded_position) {

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int x = 0; x < fs.n_pad_factor*fs.n; ++x)
		for (int y = 0; y < fs.n_pad_factor*fs.n; ++y)
		for (int z = 0; z < fs.n_pad_factor*fs.n; ++z) {
			int idx = z + pldl*(y + fs.n_pad_factor*fs.n*x);
			data[idx] /= v;
		}
	}
}

template <typename R>
void field<R>::switch_state(field_state state_, bool mmo)
{
	bool do_p2m = is_in_place() || !mmo;

	if (state_ == uninitialized) {
		state = uninitialized;
		return;
	}
	else if (state == uninitialized) {
		state = state_;
		return;
	}
	else if (state_ == state) {
		return;
	}
	else if (fs.n_pad_factor == 1 && (
		(state_ == padded_position && state == position) ||
		(state_ == position && state == padded_position) ||
		(state_ == padded_momentum && state == momentum) ||
		(state_ == momentum && state == padded_momentum))) {
		state = state_;
	}

	if (state == position) {
		state = momentum;
		if (do_p2m) p2m_plan.execute();
	}
	else if (state == padded_position) {
		state = padded_momentum;
		if (do_p2m) padded_p2m_plan.execute();
	}

switch_momentum_states:			
	if (state == momentum) {	
		if (state_ == padded_momentum || state_ == padded_position) {
			state = padded_momentum;
			pad_momentum_grid();
		}
		else if (state_ == position) {
			state = position;
			if (mdata_saved) memcpy(mdata_saved, mdata, 2*sizeof(R)*fs.total_momentum_gridpoints);
			m2p_plan.execute();
			divby(fs.total_gridpoints);
		}
	}

	if (state == padded_momentum) {
		if (state_ == momentum || state_ == position) {
			state = momentum;
			unpad_momentum_grid();
			
			if (state_ == position) {
				goto switch_momentum_states;
			}
		}
		else if (state_ == padded_position) {
			state = padded_position;
			if (mdata_saved) memcpy(mdata_saved, mdata, 2*sizeof(R)*fs.total_padded_momentum_gridpoints);
			padded_m2p_plan.execute();
			divby(fs.total_padded_gridpoints);
		}
	}
}

/*
 * The momentum-grid can be padded in place: The process works backwards.
 */

template <typename R>
void field<R>::pad_momentum_grid()
{
	if (fs.n_pad_factor < 2) {
		return;
	}

	for (int x = fs.n_pad_factor * fs.n - 1; x >= 0; --x) {
		int ux = -1;
		if (x <= fs.n/2) {
			ux = x;
		}
		else if (x > fs.n_pad_factor*fs.n - fs.n/2) {
			ux = x - (fs.n_pad_factor - 1) * fs.n;
		}

		for (int y = fs.n_pad_factor * fs.n - 1; y >= 0; --y) {
			int uy = -1;
			if (y <= fs.n/2) {
				uy = y;
			}
			else if (y > fs.n_pad_factor*fs.n - fs.n/2) {
				uy = y - (fs.n_pad_factor - 1) * fs.n;
			}

			for (int z = (fs.n_pad_factor*fs.n)/2; z >= 0; --z) {
				int idx = z+((fs.n_pad_factor*fs.n)/2 + 1)*(y + fs.n_pad_factor*fs.n*x);
				if (ux >= 0 && uy >= 0 && z <= fs.n/2) {
					int uidx = z + (fs.n/2 + 1)*(uy + fs.n*ux);
					mdata[idx][0] = pow<3>(fs.n_pad_factor) * mdata[uidx][0];
					mdata[idx][1] = pow<3>(fs.n_pad_factor) * mdata[uidx][1];

					if (z == fs.n/2) {
						mdata[idx][0] *= 0.5;
						mdata[idx][1] *= 0.5;
					}

					//cout << x << " " << y << " " << z << " <- " << ux << " " << uy << " " << z << " (" << idx << " <- " << uidx << ")" << endl;
				}
				else {
					mdata[idx][0] = mdata[idx][1] = 0;
					//cout << x << " " << y << " " << z << " (" << idx << ") = 0" << endl;
				}
			}
		}
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n_pad_factor*fs.n; ++x) {
		int x_conj = (x == 0 ? 0 : fs.n_pad_factor*fs.n - x);

		int ux = -1;
		if (x <= fs.n/2) {
			ux = x;
		}
		else if (x > fs.n_pad_factor*fs.n - fs.n/2) {
			ux = x - (fs.n_pad_factor - 1) * fs.n;
		}

		int ux_conj = -1;
		if (x_conj <= fs.n/2) {
			ux_conj = x_conj;
		}
		else if (x_conj > fs.n_pad_factor*fs.n - fs.n/2) {
			ux_conj = x_conj - (fs.n_pad_factor - 1) * fs.n;
		}

		for (int y = 0; y < fs.n_pad_factor*fs.n; ++y) {
			int y_conj = (y == 0 ? 0 : fs.n_pad_factor*fs.n - y);

			int uy = -1;
			if (y <= fs.n/2) {
				uy = y;
			}
			else if (y > fs.n_pad_factor*fs.n - fs.n/2) {
				uy = y - (fs.n_pad_factor - 1) * fs.n;
			}

			int uy_conj = -1;
			if (y_conj <= fs.n/2) {
				uy_conj = y_conj;
			}
			else if (y_conj > fs.n_pad_factor*fs.n - fs.n/2) {
				uy_conj = y_conj - (fs.n_pad_factor - 1) * fs.n;
			}

			// Set the z = 0 and z = N/2 modes by conj. symmetry.
			for (int z = 0; z < fs.n_pad_factor*fs.n; z += (fs.n_pad_factor*fs.n)/2) {
				if (ux >= 0 && uy >= 0 && (ux_conj < 0 || uy_conj < 0) && z <= fs.n/2) {
					int idx = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y + fs.n_pad_factor*fs.n*x);
					int idx_conj = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y_conj + fs.n_pad_factor*fs.n*x_conj);

					if (idx != idx_conj) {
						//cout << "mul " << x << " " << y << " " << z << " (" << idx << ") by 0.5, set " << idx_conj << " = conj idx" << endl;
						//cout << mdata[idx][0] << " " << mdata[idx][1] << " " << mdata[idx_conj][0] << " " << mdata[idx_conj][1] << endl;

						mdata[idx][0] *= 0.5;
						mdata[idx][1] *= 0.5;

						mdata[idx_conj][0] = mdata[idx][0];
						mdata[idx_conj][1] = -mdata[idx][1];

					}
				}
				else if (ux >= 0 && uy >= 0 && ux_conj >= 0 && uy_conj >= 0 && z <= fs.n/2) {
					int idx = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y + fs.n_pad_factor*fs.n*x);
					int idx_conj = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y_conj + fs.n_pad_factor*fs.n*x_conj);

					if (mdata[idx_conj][0] != mdata[idx][0] || mdata[idx_conj][1] != -mdata[idx][1]) {
						//cout << "NOT CONJ: " << idx << " " << idx_conj << " " << mdata[idx][0] << " " << mdata[idx][1] << " " << mdata[idx_conj][0] << " " << mdata[idx_conj][1] << endl;
					}
				}
			}
		}
	}
}

template <typename R>
void field<R>::unpad_momentum_grid()
{
	if (fs.n_pad_factor < 2) {
		return;
	}

	// Modes which are halved during the padding must be re-doubled here...
#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int x = 0; x < fs.n_pad_factor*fs.n; ++x) {
		int x_conj = (x == 0 ? 0 : fs.n_pad_factor*fs.n - x);

		int ux = -1;
		if (x <= fs.n/2) {
			ux = x;
		}
		else if (x > fs.n_pad_factor*fs.n - fs.n/2) {
			ux = x - (fs.n_pad_factor - 1) * fs.n;
		}

		int ux_conj = -1;
		if (x_conj <= fs.n/2) {
			ux_conj = x_conj;
		}
		else if (x_conj > fs.n_pad_factor*fs.n - fs.n/2) {
			ux_conj = x_conj - (fs.n_pad_factor - 1) * fs.n;
		}

		for (int y = 0; y < fs.n_pad_factor*fs.n; ++y) {
			int y_conj = (y == 0 ? 0 : fs.n_pad_factor*fs.n - y);

			int uy = -1;
			if (y <= fs.n/2) {
				uy = y;
			}
			else if (y > fs.n_pad_factor*fs.n - fs.n/2) {
				uy = y - (fs.n_pad_factor - 1) * fs.n;
			}

			int uy_conj = -1;
			if (y_conj <= fs.n/2) {
				uy_conj = y_conj;
			}
			else if (y_conj > fs.n_pad_factor*fs.n - fs.n/2) {
				uy_conj = y_conj - (fs.n_pad_factor - 1) * fs.n;
			}

			// Set the z = 0 and z = N/2 modes by conj. symmetry.
			for (int z = 0; z < fs.n_pad_factor*fs.n; z += (fs.n_pad_factor*fs.n)/2) {
				if (ux >= 0 && uy >= 0 && (ux_conj < 0 || uy_conj < 0) && z <= fs.n/2) {
					int idx = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y + fs.n_pad_factor*fs.n*x);
					int idx_conj = z + ((fs.n_pad_factor*fs.n)/2 + 1) * (y_conj + fs.n_pad_factor*fs.n*x_conj);

					if (idx != idx_conj) {
						mdata[idx][0] *= 2.0;
						mdata[idx][1] *= 2.0;

						mdata[idx_conj][0] = 0;
						mdata[idx_conj][1] = 0;

						//cout << "mul " << x << " " << y << " " << z << " (" << idx << ") by 2, set " << x_conj << " " << y_conj << " " << z << " (" << idx_conj << ") = 0" << endl;
					}
				}
			}
		}
	}

	for (int x = 0; x < fs.n_pad_factor * fs.n; ++x) {
		int ux = -1;
		if (x <= fs.n/2) {
			ux = x;
		}
		else if (x > fs.n_pad_factor*fs.n - fs.n/2) {
			ux = x - (fs.n_pad_factor - 1) * fs.n;
		}

		for (int y = 0; y < fs.n_pad_factor * fs.n; ++y) {
			int uy = -1;
			if (y <= fs.n/2) {
				uy = y;
			}
			else if (y > fs.n_pad_factor*fs.n - fs.n/2) {
				uy = y - (fs.n_pad_factor - 1) * fs.n;
			}

			for (int z = 0; z < (fs.n_pad_factor*fs.n)/2 + 1; ++z) {
				int idx = z+((fs.n_pad_factor*fs.n)/2 + 1)*(y + fs.n_pad_factor*fs.n*x);
				if (ux >= 0 && uy >= 0 && z <= fs.n/2) {
					int uidx = z + (fs.n/2 + 1)*(uy + fs.n*ux);
					mdata[uidx][0] = mdata[idx][0] / pow<3>(fs.n_pad_factor);
					mdata[uidx][1] = mdata[idx][1] / pow<3>(fs.n_pad_factor);

					if (z == fs.n/2) {
						mdata[uidx][0] *= 2;
						mdata[uidx][1] *= 2;
					}

					//cout << x << " " << y << " " << z << " -> " << ux << " " << uy << " " << z << " (" << idx << " -> " << uidx << ")" << endl;
				}
			}
		}
	}
}

// Explicit instantiations
template class field<double>;
#ifdef USE_LD
template class field<long double>;
#endif
