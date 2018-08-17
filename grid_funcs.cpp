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

#include "pow/pow.hpp"
#include "grid_funcs.hpp"

using namespace std;

template <typename R>
R grid_funcs<R>::compute_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return phi;
}

template <typename R>
R grid_funcs<R>::compute_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return chi;
}

template <typename R>
R grid_funcs<R>::compute_gpot(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return gpot;
}

template <typename R>
R grid_funcs<R>::compute_V_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return mp.V(phi, chi, ts.a);
}

template <typename R>
R grid_funcs<R>::compute_V(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_V_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_T_phi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return 0.5 * pow<2>(phidot) - mp.rescale_r * ts.adot/ts.a * phi*phidot +
			pow<2>(mp.rescale_r * ts.adot/ts.a ) * 0.5 * pow<2>(phi);
}

template <typename R>
R grid_funcs<R>::compute_T_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_T_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_T_chi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return 0.5 * pow<2>(chidot) - mp.rescale_r * ts.adot/ts.a * chi*chidot +
			pow<2>(mp.rescale_r * ts.adot/ts.a ) * 0.5 * pow<2>(chi);
}

template <typename R>
R grid_funcs<R>::compute_T_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_T_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_G_phi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return pow(ts.a, -2.*mp.rescale_s - 2.) * 0.5 * (
			pow<2>(phigradx) + pow<2>(phigrady) + pow<2>(phigradz)
		);
}

template <typename R>
R grid_funcs<R>::compute_G_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_G_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_G_chi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return pow(ts.a, -2.*mp.rescale_s - 2.) * 0.5 * (
			pow<2>(chigradx) + pow<2>(chigrady) + pow<2>(chigradz)
		);
}

template <typename R>
R grid_funcs<R>::compute_G_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_G_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_rho_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_V_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_T_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_T_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_G_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_G_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_rho(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_rho_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

template <typename R>
R grid_funcs<R>::compute_p_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return -compute_V_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_T_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) +
		compute_T_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot) -
		compute_G_phi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot)/3. -
		compute_G_chi_phys(fs, mp, ts, phi, chi, phidot, chidot,
		phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot)/3.;
}

template <typename R>
R grid_funcs<R>::compute_p(field_size &fs, model_params<R> &mp, time_state<R> &ts,
		R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
		R phigrady, R chigrady, R phigradz, R chigradz, R gpot)
{
	return compute_energy_scaling(mp, ts) *
		compute_p_phys(fs, mp, ts, phi, chi, phidot, chidot,
			phigradx, chigradx, phigrady, chigrady, phigradz, chigradz, gpot);
}

// Explicit instantiations
template class grid_funcs<double>;
#ifdef USE_LD
template class grid_funcs<long double>;
#endif
