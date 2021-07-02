#ifndef PRIVATE_GLOBALS_H
#define PRIVATE_GLOBALS_H

#include "../field.hpp"
#include "../field_size.hpp"
#include "../model_params.hpp"
#include "../time_state.hpp"

#include <iostream>
#include <fstream>

template <typename R>
struct private_globals
{
	// for the perturbations

	R hub;

	field<R> h[6], l[6], T_gw[6], dfld_sq[12];
	field<R> &phi, &chi;

	void evolve_perts(field_size &fs, model_params<R> &mp, time_state<R> &ts); //evolves the metric perturbations and derivatives
	void jtg_spectrum(field_size &fs, model_params<R> &mp, time_state<R> &ts); //take the power spectrum
	void perts_set_sf(model_params<R> &mp, time_state<R> &ts);
	void perts_info_file(std::ostream &info_file);

	void compute_dfld_sq(field_size &fs, model_params<R> &mp);
	void fft_stresstensor(field_size &fs, model_params<R> &mp, time_state<R> &ts);

	/**********************************************************/

	private_globals(field<R> &phi_, field<R> &chi_)
		: phi(phi_), chi(chi_), tstart(50.0), egw(false), tgw(2.0), nopad(false) {}

	inline void initialize(field_size &fs) {
		if (egw) {
			field_size upfs(fs.n);

			for (int i = 0; i < 6; ++i) {
				h[i].construct(upfs);
				l[i].construct(upfs);
				T_gw[i].construct(upfs);
			}

			for (int i = 0; i < 12; ++i) {
				dfld_sq[i].construct(nopad ? upfs : fs);
			}
		}
	}

	inline void set_sf_info(model_params<R> &mp, time_state<R> &ts) {
		perts_set_sf(mp, ts);
	}

	inline void evolve(field_size &fs, model_params<R> &mp, time_state<R> &ts, int counter) {
		int numstepsgw = (int)(tgw/ts.dt); //when to calc. gw spect

		if (egw && ts.t >= tstart) {
			evolve_perts(fs, mp, ts);
		}

		if (egw && ts.t >= tstart && counter % numstepsgw == 0) {
			jtg_spectrum(fs, mp, ts);
		}
	}

	inline void info_file_output(std::ostream &info_file) {
		perts_info_file(info_file);
	}

	const char **get_opt_names();

	bool process_opt(const char *param, const char *value);

	R tstart; //20.;
	bool egw; // compute the gw spectrum
	R tgw; //time to take gravitational wave spectrum
	bool nopad;

	std::ofstream ofs;
};

#endif

