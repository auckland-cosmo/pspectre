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
#include "fft.hpp"
#include "le_style_initializer.hpp"

#include <cstdlib>
#include <cmath>

using namespace std;

// Sets the initial amplitude of a momentum mode with indices px, py, and pz.
// Refer to the LatticeEasy manual sections 6.3.2 through 6.3.4
template <typename R>
void le_style_initializer<R>::set_mode(int px, int py, int pz, int idx, bool real)
{
	R omega; // Effective mass of the fields
	R rms_amplitude, amplitude;
	R phase1, phase2;

	//Set phi mode
	omega = sqrt(pow<2>(mp.dp)* (pow<2>(px) + pow<2>(py) + pow<2>(pz)) + pow<2>(mp.m_phi/mp.rescale_B));
	if (omega > 0.) {
		rms_amplitude = fluctuation_amplitude / (sqrt(omega));
	}
	else {
		rms_amplitude = 0.;
	}
	
	// The extra factor of 1/sqrt(2) accounts for the fact that both
	// left-moving and right-moving waves are present.
	// See the section in the LatticeEasy manual on "preserving isotropy"
	amplitude = rms_amplitude / sqrt(2.) * sqrt(log(1./drand48())); 

	phase1 = 2 * M_PI * drand48();
	phase2 = 2 * M_PI * drand48();

	phi.mdata[idx][0] = amplitude * (cos(phase1) + cos(phase2));
	phi.mdata[idx][1] = amplitude * (sin(phase1) + sin(phase2));

	phidot.mdata[idx][0] = omega * amplitude * (sin(phase1) - sin(phase2)) +
		(mp.rescale_r - 1) * adot * phi.mdata[idx][0];
	phidot.mdata[idx][1] = -1. * omega * amplitude * (cos(phase1) - cos(phase2)) +
		(mp.rescale_r - 1) * adot * phi.mdata[idx][1];

	//Set chi mode
	omega = sqrt(pow<2>(mp.dp) * (pow<2>(px) + pow<2>(py) + pow<2>(pz)) + pow<2>(mp.m_chi/mp.rescale_B));
	if (omega > 0.) {
		rms_amplitude = fluctuation_amplitude / (sqrt(omega));
	}
	else {
		rms_amplitude = 0.;
	}

	amplitude = rms_amplitude / sqrt(2.) * sqrt(log(1./drand48()));

	phase1 = 2. * M_PI * drand48();
	phase2 = 2. * M_PI * drand48();

	chi.mdata[idx][0] = amplitude * (cos(phase1) + cos(phase2));
	chi.mdata[idx][1] = amplitude * (sin(phase1) + sin(phase2));

	chidot.mdata[idx][0] = omega * amplitude * (sin(phase1) - sin(phase2)) +
		(mp.rescale_r - 1) * adot * chi.mdata[idx][0];
	chidot.mdata[idx][1] = -1. * omega * amplitude * (cos(phase1) - cos(phase2)) +
		(mp.rescale_r - 1) * adot * chi.mdata[idx][1];

	if (real) // The 8 corners of the box are real-valued.
	{
		phi.mdata[idx][1] = 0.;
		chi.mdata[idx][1] = 0.;
		phidot.mdata[idx][1] = 0.;
		chidot.mdata[idx][1] = 0.;
	}
}

template <typename R>
void le_style_initializer<R>::initialize()
{
	phi.switch_state(momentum);
	phidot.switch_state(momentum);
	chi.switch_state(momentum);
	chidot.switch_state(momentum);

	for (int x = 0; x < fs.n; ++x) {
		int x_conj = (x == 0 ? 0 : fs.n - x);
		int px = (x <= fs.n/2 ? x : x - fs.n);

		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);

			for (int z = 1; z < fs.n/2; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n*x);

				set_mode(px, py, pz, idx, false);
			}
			
			// Different cases for the z = 0 and z = N/2 modes (making use of conj symmetry)
			if (y > fs.n/2 || (x > fs.n/2 && (y == 0 || y == fs.n/2))) {
				int y_conj = (y == 0 ? 0 : fs.n - y);

				//Set the z = 0 mode for these x and y indices
				int z = 0;
				int pz = z;
				int idx = z + (fs.n/2 + 1) * (y + fs.n*x);
				int idx_conj = z + (fs.n/2 + 1) * (y_conj + fs.n*x_conj);
					  
				set_mode(px, py, pz, idx, false);

				phi.mdata[idx_conj][0] = phi.mdata[idx][0];
				phi.mdata[idx_conj][1] = -1. * phi.mdata[idx][1];
				phidot.mdata[idx_conj][0] = phidot.mdata[idx][0];
				phidot.mdata[idx_conj][1] = -1. * phidot.mdata[idx][1];

				chi.mdata[idx_conj][0] = chi.mdata[idx][0];
				chi.mdata[idx_conj][1] = -1. * chi.mdata[idx][1];
				chidot.mdata[idx_conj][0] = chidot.mdata[idx][0];
				chidot.mdata[idx_conj][1] = -1. * chidot.mdata[idx][1];

				// Set the z = N / 2 mode
				z = fs.n/2;
				pz = z;
				idx = z + (fs.n/2 + 1) * (y + fs.n*x);
				idx_conj = z + (fs.n/2 + 1) * (y_conj + fs.n*x_conj);

				set_mode(px, py, pz, idx, false);

				phi.mdata[idx_conj][0] = phi.mdata[idx][0];
				phi.mdata[idx_conj][1] = -1. * phi.mdata[idx][1];
				phidot.mdata[idx_conj][0] = phidot.mdata[idx][0];
				phidot.mdata[idx_conj][1] = -1. * phidot.mdata[idx][1];

				chi.mdata[idx_conj][0] = chi.mdata[idx][0];
				chi.mdata[idx_conj][1] = -1. * chi.mdata[idx][1];
				chidot.mdata[idx_conj][0] = chidot.mdata[idx][0];
				chidot.mdata[idx_conj][1] = -1. * chidot.mdata[idx][1];
			}
			else if ((x == 0 || x == fs.n/2) && (y == 0 || y == fs.n/2)) {
				int z = 0;
				int pz = z;
				int idx = z + (fs.n/2 + 1) * (y + fs.n*x);

				set_mode(px, py, pz, idx, true);
				
				z = fs.n/2;
				pz = z;
				idx = z + (fs.n/2 + 1) * (y + fs.n*x);
				
				set_mode(px, py, pz, idx, true);
			}
		}
	}
}

template class le_style_initializer<double>;
#ifdef USE_LD
template class le_style_initializer<long double>;
#endif
