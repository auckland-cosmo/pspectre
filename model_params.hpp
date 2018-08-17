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
 * @brief The physical model parameters.
 */

#ifndef MODEL_PARAMS
#define MODEL_PARAMS

#include "pow/pow.hpp"

#include <iostream>
#include <string>
#include <set>

#include <cmath>

template <typename R>
struct rs_init
{
	rs_init(R m, R B, R s, R r, R A, const std::string &d)
		: mag(m), rescale_B(B), rescale_s(s),
		  rescale_r(r), rescale_A(A), desc(d) {}

	// Order largest first.
	bool operator < (const rs_init &rs) const
	{
		return mag > rs.mag;
	}

	R mag;
	R rescale_B;
	R rescale_s;
	R rescale_r;
	R rescale_A;
	std::string desc;
};

template <typename R>
inline bool operator < (const rs_init<R> &rs1, const rs_init<R> &rs2)
{
	return rs1.mag > rs2.mag;
}

/**
 * @brief Static model parameters.
 */

template <typename R>
struct model_params
{
	model_params()
	{
		using namespace std;
		
		// The defaults are compatible with DEFROST
		len = 10.0;
		lambda_phi = 0.0;
		lambda_chi = 0.0;
		gamma_phi = 0.0;
		gamma_chi = 0.0;
		m_phi = (1./2.e5)/sqrt(8. * M_PI);
		m_chi = 0.0;
		g = sqrt(1e4 * 8. * M_PI * pow<2>(m_phi));
		phi0 = 1.0093430384226378929425913902459/sqrt(8. * M_PI);
		chi0 = 0.;
		phidot0 = 0.;
		chidot0 = 0.;

		calculate_derived_params();
	}

	void calculate_derived_params(bool report = false)
	{
		using namespace std;
		
		dp = 2.*M_PI/len;

		// We assume that either a mass term or a lambda term dominates.
		R phi2 = m_phi/2.0 * pow<2>(phi0);
		R chi2 = m_chi/2.0 * pow<2>(chi0);
		R phi4 = lambda_phi/4.0 * pow<4>(phi0);
		R chi4 = lambda_chi/4.0 * pow<4>(chi0);
		R phi6 = gamma_phi/6.0 * pow<6>(phi0);
		R chi6 = gamma_chi/6.0 * pow<6>(chi0);

		R phi2mag = fabs(phi2);
		R chi2mag = fabs(chi2);
		R phi4mag = fabs(phi4);
		R chi4mag = fabs(chi4);
		R phi6mag = fabs(phi6);
		R chi6mag = fabs(chi6);

		std::set<rs_init<R> > rss;
		if (phi0 != 0) rss.insert(rs_init<R>(phi2mag, m_phi, 0, 1.5, 1/phi0, "phi^2"));
		if (chi0 != 0) rss.insert(rs_init<R>(chi2mag, m_chi, 0, 1.5, 1/chi0, "chi^2"));
		if (phi0 != 0) rss.insert(rs_init<R>(phi4mag, phi0 * sqrt(fabs(lambda_phi)), -1, 1, 1/phi0, "phi^4"));
		if (chi0 != 0) rss.insert(rs_init<R>(chi4mag, chi0 * sqrt(fabs(lambda_chi)), -1, 1, 1/chi0, "chi^4"));
		if (phi0 != 0) rss.insert(rs_init<R>(phi6mag, pow<2>(phi0) * sqrt(fabs(gamma_phi)), -3./2., 3./4., 1/phi0, "phi^6"));
		if (chi0 != 0) rss.insert(rs_init<R>(chi6mag, pow<2>(chi0) * sqrt(fabs(gamma_chi)), -3./2., 3./4., 1/chi0, "chi^6"));

		if (rss.size()) {
			rescale_B = rss.begin()->rescale_B;
			rescale_s = rss.begin()->rescale_s;
			rescale_r = rss.begin()->rescale_r;
			rescale_A = rss.begin()->rescale_A;
		}
		else {
			rescale_B = 1.0;
			rescale_s = 0;
			rescale_r = 0;
			rescale_A = 1.0;
		}

		if (report) std::cout << "Initially Dominant Term: " << rss.begin()->desc << std::endl;
	}

	/**
	 * Returns the value of the field potential at a point given the values of the fields at that point.
	 * The field values are sent in program units, and the potential is returned in program units.
	 * This is equation 6.5 from the LatticeEasy manual.
	 */

	R V(R phi, R chi, R a_t)
	{
		using namespace std;

		const R tophys = 1./rescale_A * pow(a_t, -rescale_r);
		const R phi_phys = tophys * phi;
		const R chi_phys = tophys * chi;
		return pow<2>(rescale_A / rescale_B) * pow(a_t, -2. * rescale_s + 2. * rescale_r) *
			(
				0.5*pow<2>(m_phi * phi_phys) +
				0.5*pow<2>(m_chi * chi_phys) +
				0.25*lambda_phi*pow<4>(phi_phys) +
				0.25*lambda_chi*pow<4>(chi_phys) +
				0.5*pow<2>(g * phi_phys * chi_phys) +
				gamma_phi*pow<6>(phi_phys)/6.0 +
				gamma_chi*pow<6>(chi_phys)/6.0
			);
	}

	/**
	 * This is where the equations of motion for the fields are actually evaluated.
	 * The first and second time derivatives of the fields are computed in accordance
	 * with the Klein-Gordon equation, which is written in program units and
	 * transformed to momentum-space. Note that the choice of program units has eliminated
	 * the first-time-derivative term from the second-time-derivative equation.
	 */

	void derivs(R phi, R chi, R phidot, R chidot,
		R chi2phi, R phi2chi, R phi3, R chi3, R phi5, R chi5, R a_t, R adot_t, R addot_t,
		R mom2, R &dphidt, R &dchidt, R &dphidotdt, R &dchidotdt)
	{
		using namespace std;
		
		dphidt = phidot;
		dchidt = chidot;

		dphidotdt = -pow(a_t, -2. * rescale_s - 2.) * mom2 * phi +
			rescale_r * ((rescale_s - rescale_r + 2) * pow<2>(adot_t/a_t) + addot_t/a_t)*phi -
			pow(a_t, -2.*rescale_s - 2. * rescale_r)/pow<2>(rescale_B)*(pow<2>(m_phi) * pow(a_t, 2. * rescale_r) * phi +
			lambda_phi/pow<2>(rescale_A) * phi3 +
			pow<2>(g/rescale_A)*chi2phi + gamma_phi/pow<4>(rescale_A) * pow(a_t, -2. * rescale_r) * phi5);

		dchidotdt = -pow(a_t, -2. * rescale_s - 2.) * mom2 * chi +
			rescale_r * ((rescale_s - rescale_r + 2) * pow<2>(adot_t/a_t) + addot_t/a_t)*chi -
			pow(a_t, -2.*rescale_s - 2. * rescale_r)/pow<2>(rescale_B)*(pow<2>(m_chi) * pow(a_t, 2. * rescale_r) * chi +
			lambda_chi/pow<2>(rescale_A) * chi3 +
			pow<2>(g/rescale_A)*phi2chi + gamma_chi/pow<4>(rescale_A) * pow(a_t, -2. * rescale_r) * chi5);
	}

	/**
	 * Returns the second time derivative of the scale factor in program units.
	 * See equation 6.26 of the LatticeEasy manual.
	 */

	R adoubledot(R a_t, R adot_t, R avg_gradient_phi, R avg_gradient_chi, R avg_V)
	{
		using namespace std;
		
		return (
			(-rescale_s - 2.)*pow<2>(adot_t)/a_t +
			8. * M_PI / pow<2>(rescale_A) * pow(a_t, -2.* rescale_s - 2. * rescale_r - 1.)*(
				1./3. * (avg_gradient_phi + avg_gradient_chi) + pow(a_t, 2*rescale_s + 2.)*avg_V
			)
		);
	}

	/** Returns the second time derivative of the scale factor in program units at a half-time-step.
	 * See equation 6.35/6.36 of the LatticeEasy manual.
	 */

	R adoubledot_staggered(R dt, R a_t, R adot_t, R avg_gradient_phi, R avg_gradient_chi, R avg_V)
	{
		using namespace std;
		
		return (
			-2. * adot_t - 2. * a_t / (dt * (rescale_s + 2.)) *
			(1. - sqrt(1. + 2. * dt * (rescale_s + 2.) * adot_t/a_t +
			pow<2>(dt) * (rescale_s + 2.) * 8. * M_PI / pow<2>(rescale_A) * pow(a_t, -2.* rescale_s - 2. * rescale_r - 2.)*(
				1./3. * (avg_gradient_phi + avg_gradient_chi) + pow(a_t, 2*rescale_s + 2.)*avg_V
			)))
		)/dt;
	}

	R gamma_phi, gamma_chi;
	R lambda_phi, lambda_chi, g;
	R m_phi, m_chi;
	R len;
	R phi0, chi0;
	R phidot0, chidot0;
	R rescale_A, rescale_B, rescale_s, rescale_r;
	R dp;
};

#endif // MODEL_PARAMS
