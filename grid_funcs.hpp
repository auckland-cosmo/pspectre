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
 * @brief Grid point functions used for slice output, etc.
 */

#ifndef GRID_FUNCS_HPP
#define GRID_FUNCS_HPP

#include "pow/pow.hpp"
#include "model_params.hpp"
#include "field_size.hpp"
#include "time_state.hpp"

#include <cmath>

template <typename R>
struct grid_funcs {
	static R compute_energy_scaling(model_params<R> &mp, time_state<R> &ts)
	{
		using namespace std;
		
		return 8. * M_PI /(3. * pow<2>(mp.rescale_A) * pow(ts.a, 2.*mp.rescale_r) * pow<2>(ts.adot/ts.a));
	}

	// These functions all have the same signature for use with the slice output classes.
	static R compute_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_phidot(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_chidot(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_gpot(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_V_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_V(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_T_phi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_T_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_T_chi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_T_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_phys_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_phys_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_phys_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_phys_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_phys_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_phi_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_phys_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_G_chi_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_phi_phys_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_chi_phys_x(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_phi_phys_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_chi_phys_y(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_phi_phys_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_grad_chi_phys_z(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_rho_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_rho(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_p_phys(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);

	static R compute_p(field_size &fs, model_params<R> &mp, time_state<R> &ts,
			R phi, R chi, R phidot, R chidot, R phigradx, R chigradx,
			R phigrady, R chigrady, R phigradz, R chigradz, R gpot);
};

#endif // GRID_FUNCS_HPP
