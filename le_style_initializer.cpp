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
void le_style_initializer<R>::set_mode(field<R> &fld, field<R> &flddot, R m_fld_eff, int px, int py, int pz, int idx, bool real)
{
	R omega; // Effective mass of the fields
	R rms_amplitude, amplitude;
	R phase1, phase2;

	omega = sqrt(pow<2>(mp.dp)* (pow<2>(px) + pow<2>(py) + pow<2>(pz)) + pow<2>(m_fld_eff));
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

	// cout << "setting mode using " << amplitude << " " << phase1 << " " << phase2 << endl;

	fld.mdata[idx][0] = amplitude * (cos(phase1) + cos(phase2));
	fld.mdata[idx][1] = amplitude * (sin(phase1) + sin(phase2));

	flddot.mdata[idx][0] = omega * amplitude * (sin(phase1) - sin(phase2)) +
		(mp.rescale_r - 1) * adot * fld.mdata[idx][0];
	flddot.mdata[idx][1] = -1. * omega * amplitude * (cos(phase1) - cos(phase2)) +
		(mp.rescale_r - 1) * adot * fld.mdata[idx][1];

	if (real) // The 8 corners of the box are real-valued.
	{
		fld.mdata[idx][1] = 0.;
		flddot.mdata[idx][1] = 0.;
	}
}

template <typename R>
void le_style_initializer<R>::initialize_field(field<R> &fld, field<R> &flddot, R m_fld_eff)
{
	fld.switch_state(momentum);
	flddot.switch_state(momentum);

	for (int x = 0; x < fs.n; ++x) {
		int x_conj = (x == 0 ? 0 : fs.n - x);
		int px = (x <= fs.n/2 ? x : x - fs.n);

		for (int y = 0; y < fs.n; ++y) {
			int py = (y <= fs.n/2 ? y : y - fs.n);

			for (int z = 1; z < fs.n/2; ++z) {
				int pz = z;
				int idx = z + (fs.n/2+1)*(y + fs.n*x);

				// cout << "\tmode: " << x << " " << y << " " << z << endl;
				set_mode(fld, flddot, m_fld_eff, px, py, pz, idx, false);
			}
			
			// Different cases for the z = 0 and z = N/2 modes (making use of conj symmetry)
			if (y > fs.n/2 || (x > fs.n/2 && (y == 0 || y == fs.n/2))) {
				int y_conj = (y == 0 ? 0 : fs.n - y);

				//Set the z = 0 mode for these x and y indices
				int z = 0;
				int pz = z;
				int idx = z + (fs.n/2 + 1) * (y + fs.n*x);
				int idx_conj = z + (fs.n/2 + 1) * (y_conj + fs.n*x_conj);

				// cout << "\tmode: " << x << " " << y << " " << z << endl;
				set_mode(fld, flddot, m_fld_eff, px, py, pz, idx, false);

				fld.mdata[idx_conj][0] = fld.mdata[idx][0];
				fld.mdata[idx_conj][1] = -1. * fld.mdata[idx][1];
				flddot.mdata[idx_conj][0] = flddot.mdata[idx][0];
				flddot.mdata[idx_conj][1] = -1. * flddot.mdata[idx][1];

				// Set the z = N / 2 mode
				z = fs.n/2;
				pz = z;
				idx = z + (fs.n/2 + 1) * (y + fs.n*x);
				idx_conj = z + (fs.n/2 + 1) * (y_conj + fs.n*x_conj);

				// cout << "\tmode: " << x << " " << y << " " << z << endl;
				set_mode(fld, flddot, m_fld_eff, px, py, pz, idx, false);

				fld.mdata[idx_conj][0] = fld.mdata[idx][0];
				fld.mdata[idx_conj][1] = -1. * fld.mdata[idx][1];
				flddot.mdata[idx_conj][0] = flddot.mdata[idx][0];
				flddot.mdata[idx_conj][1] = -1. * flddot.mdata[idx][1];
			}
			else if ((x == 0 || x == fs.n/2) && (y == 0 || y == fs.n/2)) {
				int z = 0;
				int pz = z;
				int idx = z + (fs.n/2 + 1) * (y + fs.n*x);

				if (x != 0 || y != 0) {
					// cout << "\tmode: " << x << " " << y << " " << z << endl;
					set_mode(fld, flddot, m_fld_eff, px, py, pz, idx, true);
				}
				
				z = fs.n/2;
				pz = z;
				idx = z + (fs.n/2 + 1) * (y + fs.n*x);

				// cout << "\tmode: " << x << " " << y << " " << z << endl;
				set_mode(fld, flddot, m_fld_eff, px, py, pz, idx, true);
			}
		}
	}
}

template <typename R>
void le_style_initializer<R>::initialize()
{
	R m_phi_eff, m_chi_eff;

	m_phi_eff = sqrt(
		(
			(mp.md_e_phi != 0) ?
			(
				2.0*mp.md_c_phi*mp.md_e_phi*pow<2>(mp.md_s_phi)*
				(pow<2>(mp.md_s_phi) + (2.0*mp.md_e_phi - 1.0)*pow<2>(mp.phi0))*
				pow(1.0 + pow<2>(mp.phi0/mp.md_s_phi), mp.md_e_phi)
			)/pow<2>(pow<2>(mp.md_s_phi) + pow<2>(mp.phi0)) :
			pow<2>(mp.m_phi)
		) + pow<2>(mp.g*mp.chi0) + 3.0*mp.lambda_phi*pow<2>(mp.phi0) + 5.0*mp.gamma_phi*pow<4>(mp.phi0)
	)/mp.rescale_B;
	m_chi_eff = sqrt(
		(
			(mp.md_e_chi != 0) ?
			(
				2.0*mp.md_c_chi*mp.md_e_chi*pow<2>(mp.md_s_chi)*
				(pow<2>(mp.md_s_chi) + (2.0*mp.md_e_chi - 1.0)*pow<2>(mp.chi0))*
				pow(1.0 + pow<2>(mp.chi0/mp.md_s_chi), mp.md_e_chi)
			)/pow<2>(pow<2>(mp.md_s_chi) + pow<2>(mp.chi0)) :
			pow<2>(mp.m_chi)
		) + pow<2>(mp.g*mp.phi0) + 3.0*mp.lambda_chi*pow<2>(mp.chi0) + 5.0*mp.gamma_chi*pow<4>(mp.chi0)
	)/mp.rescale_B;

	initialize_field(phi, phidot, m_phi_eff);
	initialize_field(chi, chidot, m_chi_eff);
}

template class le_style_initializer<double>;
#ifdef USE_LD
template class le_style_initializer<long double>;
#endif
