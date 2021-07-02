/*
This file has the functions needed to evolve the metric perturbations.  It uses a 4th order runge-kutta integrator that ties into the LatticeEasy program.
Written by John T. Giblin, jr .. Last modified 4.22.2007
*/
// Ported to PSpectRe by Hal Finkel: 2009-2011

#include "private_globals.hpp"
#include "../model.hpp"
#include "../pow/pow.hpp"

#include <iostream>
#include <iomanip>

#include <string>
#include <cmath>
#include <cstdlib>
#include <cstring> // for memcpy

using namespace std;

#define c_re(v) ((v)[0])
#define c_im(v) ((v)[1])

template <typename R>
const char **private_globals<R>::get_opt_names()
{
	static const char * opt_names[] = {
                "gw", "tgw", "tstart", "nopad", 0
	};

	return opt_names;
};

template <typename R>
bool private_globals<R>::process_opt(const char *param, const char *value)
{
	if (!strcmp(param, "gw")) {
		egw = true;
		return true;
	}
	else if (!strcmp(param, "tgw") && value) {
		tgw = atof(value);
		return true;
	}
	else if (!strcmp(param, "tstart") && value) {
		tstart = atof(value);
		return true;
	}
	else if (!strcmp(param, "nopad")) {
		nopad = true;
		return true;
	}

	return false;
}

template <typename R>
void private_globals<R>::perts_set_sf(model_params<R> &mp, time_state<R> &ts)
{
      //jtg: added for coeffients in evolve_perts(R d);
      hub = ts.adot/ts.a;
}

template <typename R>
void private_globals<R>::perts_info_file(ostream &info_file) {
	info_file << "hij enabled: " << (egw ? string("yes") : string("no")) << endl;
	info_file << "hij start time: " << tstart << endl;
	info_file << "hij evolution interval: " << tgw << endl;
	info_file << "hij never pad grid: " << (nopad ? string("yes") : string("no")) << endl;
}

template <typename R>
void gwderivs(model_params<R> &, time_state<R> &, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, R, int, int, int, R*, R*, R, R*);

// jtg_stresstensor calculates the T_ab for any point (i,j,k). It's not really the SET, but appropriate combinations of the SET to get S_{ij}--the source term of the metric perturbations
template <typename R>
static inline void jtg_stresstensor(field_size &fs, field<R> *dfld_sq, int aa, int bb, int i, int j, int k, typename fft_dft_c2r_3d_plan<R>::complex_t &ret) {
	int vector_index;
	int pz, px, py;
	const int nflds = 2;

	px = (i<=fs.n/2 ? i : i-fs.n);
	py = (j<=fs.n/2 ? j : j-fs.n);
	pz = k;

	vector_index = k + (fs.n/2+1)*(j +fs.n*i);

	ret[0] = ret[1] = 0.0;

	for (int m = 0; m < 12; ++m) {
		dfld_sq[m].switch_state(momentum);
	}

	// The result here must be in physical units except for an additional factor of A^2 a^{2r}.
    
	if (aa==1) {
		if (bb==1) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += (1./3.)*(2*dfld_sq[3*m+0].mdata[vector_index][c] - dfld_sq[3*m+1].mdata[vector_index][c]
						- dfld_sq[3*m+2].mdata[vector_index][c]);
				}
			}	  
		}
		else if (bb==2) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += dfld_sq[3*m+6].mdata[vector_index][c];
				}
			}
		}
		else if (bb==3) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += dfld_sq[3*m+7].mdata[vector_index][c];
				}
			}
		}
		else 
			printf("calling wrong SET");
	}
	else if (aa==2) {
		if (bb==2) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += (1./3.)*(-dfld_sq[3*m+0].mdata[vector_index][c] + 2*dfld_sq[3*m+1].mdata[vector_index][c]
						- dfld_sq[3*m+2].mdata[vector_index][c]);
				}
			}	  
		}
		else if (bb==3) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += dfld_sq[3*m+8].mdata[vector_index][c];
				}
			}
		}
		else
			cerr << "calling wrong SET" << endl;
	}
	else if (aa==3) {
		if (bb==3) {
			for(int m=0; m < nflds; m++) {
				for (int c=0; c < 2; c++) {
					ret[c] += (1./3.)*(-dfld_sq[3*m+0].mdata[vector_index][c] - dfld_sq[3*m+1].mdata[vector_index][c]
						+ 2*dfld_sq[3*m+2].mdata[vector_index][c]);
				}
			}
		}
		else
			cerr << "calling wrong SET" << endl;
	}
	else {
		cerr << "calling wrong SET" << endl;
	}

	/* printf("%d %d %d %" RFP "g %" RFP "g\n", i, j, k, ret[0], ret[1]); */
}

template <typename R>
void private_globals<R>::compute_dfld_sq(field_size &fs, model_params<R> &mp) {
	phi.switch_state(momentum);
	chi.switch_state(momentum);

	for (int i = 0; i < 6; ++i) {
		dfld_sq[i].switch_state(uninitialized);
		dfld_sq[i].switch_state(momentum);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for (int i=0; i<fs.n; i++) {              // stores them in the correct
		int px = (i<=fs.n/2 ? i : i-fs.n);
		for (int j=0; j<fs.n; j++) {         // arrays
			int py = (j<=fs.n/2 ? j : j-fs.n);
			for (int k=0; k<fs.n/2+1; k++) {
				int pz = k;

				int vector_index = k + (fs.n/2+1)*(j +fs.n*i);

				dfld_sq[0].mdata[vector_index][0] = -mp.dp*px*phi.mdata[vector_index][1];
				dfld_sq[0].mdata[vector_index][1] =  mp.dp*px*phi.mdata[vector_index][0];

				dfld_sq[1].mdata[vector_index][0] = -mp.dp*py*phi.mdata[vector_index][1];
				dfld_sq[1].mdata[vector_index][1] =  mp.dp*py*phi.mdata[vector_index][0];

				dfld_sq[2].mdata[vector_index][0] = -mp.dp*pz*phi.mdata[vector_index][1];
				dfld_sq[2].mdata[vector_index][1] =  mp.dp*pz*phi.mdata[vector_index][0];

				dfld_sq[3].mdata[vector_index][0] = -mp.dp*px*chi.mdata[vector_index][1];
				dfld_sq[3].mdata[vector_index][1] =  mp.dp*px*chi.mdata[vector_index][0];

				dfld_sq[4].mdata[vector_index][0] = -mp.dp*py*chi.mdata[vector_index][1];
				dfld_sq[4].mdata[vector_index][1] =  mp.dp*py*chi.mdata[vector_index][0];

				dfld_sq[5].mdata[vector_index][0] = -mp.dp*pz*chi.mdata[vector_index][1];
				dfld_sq[5].mdata[vector_index][1] =  mp.dp*pz*chi.mdata[vector_index][0];
			}
		}
	}

#if 0
	{ /*!!*/
		for (int i = 0; i < 6; ++i) {
			dfld_sq[i].switch_state(position);
		}

		phi.switch_state(position);
		chi.switch_state(position);

		for (int k = 0; k < fs.n; ++k)
		for (int j = 0; j < fs.n; ++j)
		for (int i = 0; i < fs.n; ++i)
		for (int fld = 0; fld < 2; ++fld)
		for (int n = 0; n < 3; ++n) {
			int idx = k + phi.ldl*(j + fs.n*i);
			cout << "deriv " << fld << " " << i << " " << j << " " << k << " " << n << ": " << dfld_sq[3*fld+n].data[idx] <<
				" (" << (fld ? chi.data[idx] : phi.data[idx]) << ")" << endl;
			if (n == 2) {
				int kp = (k - 1 + fs.n) % fs.n;
				int kn = (k + 1 + fs.n) % fs.n;

				int kp2 = (k - 2 + fs.n) % fs.n;
				int kn2 = (k + 2 + fs.n) % fs.n;

				int kp3 = (k - 3 + fs.n) % fs.n;
				int kn3 = (k + 3 + fs.n) % fs.n;

				int kp4 = (k - 4 + fs.n) % fs.n;
				int kn4 = (k + 4 + fs.n) % fs.n;

				int idxp = kp + phi.ldl*(j + fs.n*i);
				int idxn = kn + phi.ldl*(j + fs.n*i);

				int idxp2 = kp2 + phi.ldl*(j + fs.n*i);
				int idxn2 = kn2 + phi.ldl*(j + fs.n*i);

				int idxp3 = kp3 + phi.ldl*(j + fs.n*i);
				int idxn3 = kn3 + phi.ldl*(j + fs.n*i);

				int idxp4 = kp4 + phi.ldl*(j + fs.n*i);
				int idxn4 = kn4 + phi.ldl*(j + fs.n*i);

				field<R> &f = fld ? chi : phi;				
				R dx = mp.len/fs.n;

				cout << "\t2nd odr: " << (f.data[idxn] - f.data[idxp])/(2*dx) << " (" << f.data[idxp] << " " << f.data[idx] << " " << f.data[idxn] << " " << dx << ")" << endl;
				cout << "\t4nd odr: " << (1/12.*f.data[idxp2] - 2/3.*f.data[idxp] + 2/3.*f.data[idxn] - 1/12.*f.data[idxn2])/dx << " (" << f.data[idxp2] << " " << f.data[idxp] << " " << f.data[idx] << " " << f.data[idxn] << " " << f.data[idxn2] << " " << dx << ")" << endl;
				cout << "\t6nd odr: " << (-1/60.*f.data[idxp3] + 3/20.*f.data[idxp2] - 3/4.*f.data[idxp] + 3/4.*f.data[idxn] - 3/20.*f.data[idxn2] + 1/60.*f.data[idxn3])/dx << " (" << f.data[idxp3] << " " << f.data[idxp2] << " " << f.data[idxp] << " " << f.data[idx] << " " << f.data[idxn] << " " << f.data[idxn2] << " " << f.data[idxn3] << " " << dx << ")" << endl;
				cout << "\t8nd odr: " << (1/280.*f.data[idxp4] - 4/105.*f.data[idxp3] + 1/5.*f.data[idxp2] - 4/5.*f.data[idxp] + 4/5.*f.data[idxn] - 1/5.*f.data[idxn2] + 4/105.*f.data[idxn3] - 1/280.*f.data[idxn4])/dx << " (" << f.data[idxp4] << " " << f.data[idxp3] << " " << f.data[idxp2] << " " << f.data[idxp] << " " << f.data[idx] << " " << f.data[idxn] << " " << f.data[idxn2] << " " << f.data[idxn3] << " " << f.data[idxn4] << " " << dx << ")" << endl;
			}
		}
	}
#endif

	for (int i = 0; i < 12; ++i) {
		if (i >= 6) dfld_sq[i].switch_state(uninitialized);
		dfld_sq[i].switch_state(padded_position);
	}

	int enpf = nopad ? 1 : fs.n_pad_factor;

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for(int i=0; i<fs.n*enpf; i++) {
		for (int j=0; j<fs.n*enpf; j++) {
			for (int k=0; k<fs.n*enpf; k++) {
				int vector_index = k + dfld_sq[0].pldl * (j + fs.n*enpf * i);

				// Then comes dx*dy, dx*dz, dy*dz
				dfld_sq[6].data[vector_index] = dfld_sq[0].data[vector_index]*dfld_sq[1].data[vector_index];
				dfld_sq[7].data[vector_index] = dfld_sq[0].data[vector_index]*dfld_sq[2].data[vector_index];
				dfld_sq[8].data[vector_index] = dfld_sq[1].data[vector_index]*dfld_sq[2].data[vector_index];

				dfld_sq[9].data[vector_index] = dfld_sq[3].data[vector_index]*dfld_sq[4].data[vector_index];
				dfld_sq[10].data[vector_index] = dfld_sq[3].data[vector_index]*dfld_sq[5].data[vector_index];
				dfld_sq[11].data[vector_index] = dfld_sq[4].data[vector_index]*dfld_sq[5].data[vector_index];

				// First 6 are just the values squared.
				for (int m=0; m<6; m++) {
					dfld_sq[m].data[vector_index] *= dfld_sq[m].data[vector_index];
				}
			}
		}
	}

	for (int i = 0; i < 12; ++i) {
		dfld_sq[i].switch_state(momentum);
	}
}

//put the sourceterm in momentum space
template <typename R>
void private_globals<R>::fft_stresstensor(field_size &fs, model_params<R> &mp, time_state<R> &ts) {
	R normmm = pow(mp.len/R(fs.n),3)/(mp.rescale_A*pow(ts.a,2+2.*mp.rescale_s+mp.rescale_r));
	// R normmm = /* pow(L/R(N),3) */ 1./(rescale_A*pow(a,2+2.*rescale_s+rescale_r));
	//coefficient in front of the SET in the euqations of motion (minus the 
	//  8\p/3)

	// cout << "normmm: " << normmm << endl;

	// The conversion from Tij is: Tij_pr = A/B^2 a^{-2s+r} Tij. The Tij array
	// here is in program units except for the 1/B^2 factor. The e.o.m. has a 1/a^2
	// factor in front of the Tij term. So the Tij are in physical units except:
	// 1/B^2 Tij_pr = A/B^2 a^{-2s+r} Tij = 1/B^2 (1/A a^{-2s-r}) (A^2 a^{2r} Tij)
	// To the tij[012345] are in physical units times A^2 a^{2r} (this is the square
	// of the field value rescaling).

	compute_dfld_sq(fs, mp);

	for (int i = 0; i < 6; ++i) {
		T_gw[i].switch_state(momentum);
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif

	for(int i=0; i<fs.n; i++) {              // stores them in the correct
		int px = (i<=fs.n/2 ? i : i-fs.n);
		for (int j=0; j<fs.n; j++) {         // arrays
			int py = (j<=fs.n/2 ? j : j-fs.n);
			for (int k=0; k<fs.n/2+1; k++) {
				typename fft_dft_c2r_3d_plan<R>::complex_t tij0 = {0, 0};
				typename fft_dft_c2r_3d_plan<R>::complex_t tij1 = {0, 0};
				typename fft_dft_c2r_3d_plan<R>::complex_t tij2 = {0, 0};
				typename fft_dft_c2r_3d_plan<R>::complex_t tij3 = {0, 0};
				typename fft_dft_c2r_3d_plan<R>::complex_t tij4 = {0, 0};
				typename fft_dft_c2r_3d_plan<R>::complex_t tij5 = {0, 0};

				jtg_stresstensor<R>(fs,dfld_sq,1,1,i,j,k,tij0);
				jtg_stresstensor<R>(fs,dfld_sq,1,2,i,j,k,tij1);
				jtg_stresstensor<R>(fs,dfld_sq,1,3,i,j,k,tij2);
				jtg_stresstensor<R>(fs,dfld_sq,2,2,i,j,k,tij3);
				jtg_stresstensor<R>(fs,dfld_sq,2,3,i,j,k,tij4);
				jtg_stresstensor<R>(fs,dfld_sq,3,3,i,j,k,tij5);

				// cout << "tij0_out: " << i << " " << j << " " << k << " " << tij0[0] << " " << tij0[1] << endl;

				int pz = k;
				R dpz, dpx, dpy;
				R momentum2=pow<2>((R)px)+pow<2>((R)py)+pow<2>((R)pz);
				if (px == 0 && py == 0 && pz == 0) {
					dpz = dpx = dpy = 0;
				}
				else {
					dpz = ((R)pz)/sqrt(momentum2);
					dpx = ((R)px)/sqrt(momentum2);
					dpy = ((R)py)/sqrt(momentum2);
				}

				int vector_index = k + (fs.n/2+1)*(j + fs.n*i);

				// dp[xyz] are unitless. The T_gw must be in program
				// units except for a missing factor of 1/B^2.

#ifdef HIJ_SRC_DEBUG
				T_gw[0].mdata[vector_index][0] = c_re(tij0);
				T_gw[0].mdata[vector_index][1] = c_im(tij0);
				T_gw[1].mdata[vector_index][0] = c_re(tij1);
				T_gw[1].mdata[vector_index][1] = c_im(tij1);
				T_gw[2].mdata[vector_index][0] = c_re(tij2);
				T_gw[2].mdata[vector_index][1] = c_im(tij2);
				T_gw[3].mdata[vector_index][0] = c_re(tij3);
				T_gw[3].mdata[vector_index][1] = c_im(tij3);
				T_gw[4].mdata[vector_index][0] = c_re(tij4);
				T_gw[4].mdata[vector_index][1] = c_im(tij4);
				T_gw[5].mdata[vector_index][0] = c_re(tij5);
				T_gw[5].mdata[vector_index][1] = c_im(tij5);
#else
				//T_{11}
				T_gw[0].mdata[vector_index][0] = .5*pow<2>(1-dpx*dpx)*normmm*c_re(tij0)
					-(1-dpx*dpx)*dpx*dpy*normmm*c_re(tij1)
					-(1-dpx*dpx)*dpx*dpz*normmm*c_re(tij2)
					+pow<2>(dpx*dpy)*normmm*c_re(tij3)
					-.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij3)
					+(1+dpx*dpx)*dpy*dpz*normmm*c_re(tij4)
					+pow<2>(dpx*dpz)*normmm*c_re(tij5)
					-.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij5);
				//T_{11}(i)
				T_gw[0].mdata[vector_index][1] = .5*pow<2>(1-dpx*dpx)*normmm*c_im(tij0)
					-(1-dpx*dpx)*dpx*dpy*normmm*c_im(tij1)
					-(1-dpx*dpx)*dpx*dpz*normmm*c_im(tij2)
					+pow<2>(dpx*dpy)*normmm*c_im(tij3)
					-.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij3)
					+(1+dpx*dpx)*dpy*dpz*normmm*c_im(tij4)
					+pow<2>(dpx*dpz)*normmm*c_im(tij5)
					-.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij5);
				//T_{22}
				T_gw[3].mdata[vector_index][0] = .5*pow<2>(1-dpy*dpy)*normmm*c_re(tij3)
					-(1-dpy*dpy)*dpx*dpy*normmm*c_re(tij1)
					-(1-dpy*dpy)*dpy*dpz*normmm*c_re(tij4)
					+pow<2>(dpx*dpy)*normmm*c_re(tij0)
					-.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij0)
					+(1+dpy*dpy)*dpx*dpz*normmm*c_re(tij2)
					+pow<2>(dpy*dpz)*normmm*c_re(tij5)
					-.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij5);
				//T_{22}(i)
				T_gw[3].mdata[vector_index][1] = .5*pow<2>(1-dpy*dpy)*normmm*c_im(tij3)
					-(1-dpy*dpy)*dpx*dpy*normmm*c_im(tij1)
					-(1-dpy*dpy)*dpy*dpz*normmm*c_im(tij4)
					+pow<2>(dpx*dpy)*normmm*c_im(tij0)
					-.5*(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij0)
					-.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij5);
				//T_{33}
				T_gw[5].mdata[vector_index][0] = .5*pow<2>(1-dpz*dpz)*normmm*c_re(tij5)
					-(1-dpz*dpz)*dpx*dpz*normmm*c_re(tij2)
					-(1-dpz*dpz)*dpy*dpz*normmm*c_re(tij4)
					+pow<2>(dpx*dpz)*normmm*c_re(tij0)
					-.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij0)
					+(1+dpz*dpz)*dpx*dpy*normmm*c_re(tij1)
					+pow<2>(dpy*dpz)*normmm*c_re(tij3)
					-.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij3);
				//T_{33}(i)
				T_gw[5].mdata[vector_index][1] = .5*pow<2>(1-dpz*dpz)*normmm*c_im(tij5)
					-(1-dpz*dpz)*dpx*dpz*normmm*c_im(tij2)
					-(1-dpz*dpz)*dpy*dpz*normmm*c_im(tij4)
					+pow<2>(dpx*dpz)*normmm*c_im(tij0)
					-.5*(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij0)
					+(1+dpz*dpz)*dpx*dpy*normmm*c_im(tij1)
					+pow<2>(dpy*dpz)*normmm*c_im(tij3)
					-.5*(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij3);
				//T_{12}
				T_gw[1].mdata[vector_index][0] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*c_re(tij0)
					-.5*dpx*dpy*(1-dpy*dpy)*normmm*c_re(tij3)
					+.5*dpx*dpy*(1+dpz*dpz)*normmm*c_re(tij5)
					+(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_re(tij1)
					-dpy*dpz*(1-dpx*dpx)*normmm*c_re(tij2)
					-dpx*dpz*(1-dpy*dpy)*normmm*c_re(tij4);
				//T_{12}(i)
				T_gw[1].mdata[vector_index][1] = -.5*dpx*dpy*(1-dpx*dpx)*normmm*c_im(tij0)
					-.5*dpx*dpy*(1-dpy*dpy)*normmm*c_im(tij3)
					+.5*dpx*dpy*(1+dpz*dpz)*normmm*c_im(tij5)
					+(1-dpx*dpx)*(1-dpy*dpy)*normmm*c_im(tij1)
					-dpy*dpz*(1-dpx*dpx)*normmm*c_im(tij2)
					-dpx*dpz*(1-dpy*dpy)*normmm*c_im(tij4);
				//T_{13}
				T_gw[2].mdata[vector_index][0] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*c_re(tij0)
					+.5*dpx*dpz*(1+dpy*dpy)*normmm*c_re(tij3)
					-.5*dpx*dpz*(1-dpz*dpz)*normmm*c_re(tij5)
					+(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_re(tij2)
					-dpy*dpz*(1-dpx*dpx)*normmm*c_re(tij1)
					-dpx*dpy*(1-dpz*dpz)*normmm*c_re(tij4);
				//T_{13}
				T_gw[2].mdata[vector_index][1] = -.5*dpx*dpz*(1-dpx*dpx)*normmm*c_im(tij0)
					+.5*dpx*dpz*(1+dpy*dpy)*normmm*c_im(tij3)
					-.5*dpx*dpz*(1-dpz*dpz)*normmm*c_im(tij5)
					+(1-dpx*dpx)*(1-dpz*dpz)*normmm*c_im(tij2)
					-dpy*dpz*(1-dpx*dpx)*normmm*c_im(tij1)
					-dpx*dpy*(1-dpz*dpz)*normmm*c_im(tij4);
				//T_{23}
				T_gw[4].mdata[vector_index][0] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*c_re(tij0)
					-.5*dpy*dpz*(1-dpy*dpy)*normmm*c_re(tij3)
					-.5*dpy*dpz*(1-dpz*dpz)*normmm*c_re(tij5)
					+(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_re(tij4)
					-dpx*dpy*(1-dpz*dpz)*normmm*c_re(tij2)
					-dpx*dpz*(1-dpy*dpy)*normmm*c_re(tij1);
				//T_{23}(i)
				T_gw[4].mdata[vector_index][1] = +.5*dpy*dpz*(1+dpx*dpx)*normmm*c_im(tij0)
					-.5*dpy*dpz*(1-dpy*dpy)*normmm*c_im(tij3)
					-.5*dpy*dpz*(1-dpz*dpz)*normmm*c_im(tij5)
					+(1-dpy*dpy)*(1-dpz*dpz)*normmm*c_im(tij4)
					-dpx*dpy*(1-dpz*dpz)*normmm*c_im(tij2)
					-dpx*dpz*(1-dpy*dpy)*normmm*c_im(tij1);
#endif // HIJ_SRC_DEBUG
			}
		}
	}    

	//subtract zero mode (just in case)
	for(int i=0; i<6; i++){
		T_gw[i].mdata[0][0] = T_gw[i].mdata[0][1] = 0;
	}
   
	//    FILTER(kmin_gw, kmax_gw, T_gw);   
	// allows for a filter on the set
}

//  this evolves the metric perturbations using 4th order RK
template <typename R>
void private_globals<R>::evolve_perts(field_size &fs, model_params<R> &mp, time_state<R> &ts) {

	fft_stresstensor(fs, mp, ts);   //calculate the source terms in momentum space

	for (int i = 0; i < 6; ++i) {
		T_gw[i].switch_state(momentum);
		h[i].switch_state(momentum);
		l[i].switch_state(momentum);
	}
 
//	printf("called fft_stresstensor \n");

#if 0
	{ /*!!*/
		R avg[12];
		for (int i = 0; i < 12; ++i) avg[i] = 0;

		for (int i = 0; i < fs.n; ++i)
		for (int j = 0; j < fs.n; ++j)
		for (int z = 0; z < 12; ++z) {
			int vector_index = fs.n/4 + (fs.n/2+1)*(j + i*fs.n);
			avg[z] += T_gw[z/2].mdata[vector_index][z % 2];
		}

		for (int i = 0; i < 12; ++i) {
			avg[i] /= fs.n*fs.n;
			cout << "<T_gw[" << i << "](k=N/4)> = " << avg[i] << endl;
		}
	}
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    
	for(int i=0;i<fs.n;i++){        // EOM's
		int ii = (i<=fs.n/2 ? i : i-fs.n);
		for(int j=0;j<fs.n;j++){
			int jj = (j<=fs.n/2 ? j : j-fs.n);
			for(int k=0;k<fs.n/2+1;k++){ //only need to do half the array
				R source[12];
				// values here are in program units except for a missing factor of 1/B^2
	
				int vector_index = k + (fs.n/2+1)*(j + i*fs.n);
	
				source[0] = T_gw[0].mdata[vector_index][0];      //take the correct parts
				source[1] = T_gw[0].mdata[vector_index][1];      //point on the lattice
				source[2] = T_gw[1].mdata[vector_index][0]; 
				source[3] = T_gw[1].mdata[vector_index][1]; 
				source[4] = T_gw[2].mdata[vector_index][0]; 
				source[5] = T_gw[2].mdata[vector_index][1]; 
				source[6] = T_gw[3].mdata[vector_index][0];    
				source[7] = T_gw[3].mdata[vector_index][1];    
				source[8] = T_gw[4].mdata[vector_index][0]; 
				source[9] = T_gw[4].mdata[vector_index][1]; 
				source[10] = T_gw[5].mdata[vector_index][0]; 
				source[11] = T_gw[5].mdata[vector_index][1]; 

				R hd1[12], hd2[12], hd3[12], hd4[12];
				R ld1[12], ld2[12], ld3[12], ld4[12];
		
				gwderivs<R>(mp, ts,
					h[0].mdata[vector_index][0],
					h[0].mdata[vector_index][1],
					h[1].mdata[vector_index][0],
					h[1].mdata[vector_index][1],
					h[2].mdata[vector_index][0],
					h[2].mdata[vector_index][1],
					h[3].mdata[vector_index][0],
					h[3].mdata[vector_index][1],
					h[4].mdata[vector_index][0],
					h[4].mdata[vector_index][1],
					h[5].mdata[vector_index][0],
					h[5].mdata[vector_index][1],
					l[0].mdata[vector_index][0],
					l[0].mdata[vector_index][1],
					l[1].mdata[vector_index][0],
					l[1].mdata[vector_index][1],
					l[2].mdata[vector_index][0],
					l[2].mdata[vector_index][1],
					l[3].mdata[vector_index][0],
					l[3].mdata[vector_index][1],
					l[4].mdata[vector_index][0],
					l[4].mdata[vector_index][1],
					l[5].mdata[vector_index][0],
					l[5].mdata[vector_index][1],
					ii, jj, k,
					hd1, ld1,
					ts.t, source);

				gwderivs<R>(mp, ts,
					h[0].mdata[vector_index][0]+ts.dt*hd1[0]/2,
					h[0].mdata[vector_index][1]+ts.dt*hd1[1]/2,
					h[1].mdata[vector_index][0]+ts.dt*hd1[2]/2,
					h[1].mdata[vector_index][1]+ts.dt*hd1[3]/2,
					h[2].mdata[vector_index][0]+ts.dt*hd1[4]/2,
					h[2].mdata[vector_index][1]+ts.dt*hd1[5]/2,
					h[3].mdata[vector_index][0]+ts.dt*hd1[6]/2,
					h[3].mdata[vector_index][1]+ts.dt*hd1[7]/2,
					h[4].mdata[vector_index][0]+ts.dt*hd1[8]/2,
					h[4].mdata[vector_index][1]+ts.dt*hd1[9]/2,
					h[5].mdata[vector_index][0]+ts.dt*hd1[10]/2,
					h[5].mdata[vector_index][1]+ts.dt*hd1[11]/2,
					l[0].mdata[vector_index][0]+ts.dt*ld1[0]/2,
					l[0].mdata[vector_index][1]+ts.dt*ld1[1]/2,
					l[1].mdata[vector_index][0]+ts.dt*ld1[2]/2,
					l[1].mdata[vector_index][1]+ts.dt*ld1[3]/2,
					l[2].mdata[vector_index][0]+ts.dt*ld1[4]/2,
					l[2].mdata[vector_index][1]+ts.dt*ld1[5]/2,
					l[3].mdata[vector_index][0]+ts.dt*ld1[6]/2,
					l[3].mdata[vector_index][1]+ts.dt*ld1[7]/2,
					l[4].mdata[vector_index][0]+ts.dt*ld1[8]/2,
					l[4].mdata[vector_index][1]+ts.dt*ld1[9]/2,
					l[5].mdata[vector_index][0]+ts.dt*ld1[10]/2,
					l[5].mdata[vector_index][1]+ts.dt*ld1[11]/2,
					ii, jj, k,
					hd2, ld2, 
					ts.t+ts.dt/2, source);

				gwderivs<R>(mp, ts,
					h[0].mdata[vector_index][0]+ts.dt*hd2[0]/2,
					h[0].mdata[vector_index][1]+ts.dt*hd2[1]/2,
					h[1].mdata[vector_index][0]+ts.dt*hd2[2]/2,
					h[1].mdata[vector_index][1]+ts.dt*hd2[3]/2,
					h[2].mdata[vector_index][0]+ts.dt*hd2[4]/2,
					h[2].mdata[vector_index][1]+ts.dt*hd2[5]/2,
					h[3].mdata[vector_index][0]+ts.dt*hd2[6]/2,
					h[3].mdata[vector_index][1]+ts.dt*hd2[7]/2,
					h[4].mdata[vector_index][0]+ts.dt*hd2[8]/2,
					h[4].mdata[vector_index][1]+ts.dt*hd2[9]/2,
					h[5].mdata[vector_index][0]+ts.dt*hd2[10]/2,
					h[5].mdata[vector_index][1]+ts.dt*hd2[11]/2,
					l[0].mdata[vector_index][0]+ts.dt*ld2[0]/2,
					l[0].mdata[vector_index][1]+ts.dt*ld2[1]/2,
					l[1].mdata[vector_index][0]+ts.dt*ld2[2]/2,
					l[1].mdata[vector_index][1]+ts.dt*ld2[3]/2,
					l[2].mdata[vector_index][0]+ts.dt*ld2[4]/2,
					l[2].mdata[vector_index][1]+ts.dt*ld2[5]/2, 
					l[3].mdata[vector_index][0]+ts.dt*ld2[6]/2,
					l[3].mdata[vector_index][1]+ts.dt*ld2[7]/2,
					l[4].mdata[vector_index][0]+ts.dt*ld2[8]/2,
					l[4].mdata[vector_index][1]+ts.dt*ld2[9]/2,
					l[5].mdata[vector_index][0]+ts.dt*ld2[10]/2,
					l[5].mdata[vector_index][1]+ts.dt*ld2[11]/2,  
					ii, jj, k,
					hd3, ld3, 
					ts.t+ts.dt/2, source);

				gwderivs<R>(mp, ts,
					h[0].mdata[vector_index][0]+ts.dt*hd3[0], 
					h[0].mdata[vector_index][1]+ts.dt*hd3[1], 
					h[1].mdata[vector_index][0]+ts.dt*hd3[2], 
					h[1].mdata[vector_index][1]+ts.dt*hd3[3], 
					h[2].mdata[vector_index][0]+ts.dt*hd3[4], 
					h[2].mdata[vector_index][1]+ts.dt*hd3[5], 
					h[3].mdata[vector_index][0]+ts.dt*hd3[6], 
					h[3].mdata[vector_index][1]+ts.dt*hd3[7], 
					h[4].mdata[vector_index][0]+ts.dt*hd3[8], 
					h[4].mdata[vector_index][1]+ts.dt*hd3[9], 
					h[5].mdata[vector_index][0]+ts.dt*hd3[10], 
					h[5].mdata[vector_index][1]+ts.dt*hd3[11], 
					l[0].mdata[vector_index][0]+ts.dt*ld3[0], 
					l[0].mdata[vector_index][1]+ts.dt*ld3[1], 
					l[1].mdata[vector_index][0]+ts.dt*ld3[2], 
					l[1].mdata[vector_index][1]+ts.dt*ld3[3], 
					l[2].mdata[vector_index][0]+ts.dt*ld3[4], 
					l[2].mdata[vector_index][1]+ts.dt*ld3[5], 
					l[3].mdata[vector_index][0]+ts.dt*ld3[6], 
					l[3].mdata[vector_index][1]+ts.dt*ld3[7], 
					l[4].mdata[vector_index][0]+ts.dt*ld3[8], 
					l[4].mdata[vector_index][1]+ts.dt*ld3[9], 
					l[5].mdata[vector_index][0]+ts.dt*ld3[10], 
					l[5].mdata[vector_index][1]+ts.dt*ld3[11], 
					ii, jj, k, 
					hd4, ld4,
					ts.t+ts.dt, source);

				h[0].mdata[vector_index][0] += ts.dt*(hd1[0]/6+hd2[0]/3+hd3[0]/3+hd4[0]/6);
				h[0].mdata[vector_index][1] += ts.dt*(hd1[1]/6+hd2[1]/3+hd3[1]/3+hd4[1]/6);
				h[1].mdata[vector_index][0] += ts.dt*(hd1[2]/6+hd2[2]/3+hd3[2]/3+hd4[2]/6);
				h[1].mdata[vector_index][1] += ts.dt*(hd1[3]/6+hd2[3]/3+hd3[3]/3+hd4[3]/6);
				h[2].mdata[vector_index][0] += ts.dt*(hd1[4]/6+hd2[4]/3+hd3[4]/3+hd4[4]/6);
				h[2].mdata[vector_index][1] += ts.dt*(hd1[5]/6+hd2[5]/3+hd3[5]/3+hd4[5]/6);
				h[3].mdata[vector_index][0] += ts.dt*(hd1[6]/6+hd2[6]/3+hd3[6]/3+hd4[6]/6);
				h[3].mdata[vector_index][1] += ts.dt*(hd1[7]/6+hd2[7]/3+hd3[7]/3+hd4[7]/6);
				h[4].mdata[vector_index][0] += ts.dt*(hd1[8]/6+hd2[8]/3+hd3[8]/3+hd4[8]/6);
				h[4].mdata[vector_index][1] += ts.dt*(hd1[9]/6+hd2[9]/3+hd3[9]/3+hd4[9]/6);
				h[5].mdata[vector_index][0] += ts.dt*(hd1[10]/6+hd2[10]/3+hd3[10]/3+hd4[10]/6);
				h[5].mdata[vector_index][1] += ts.dt*(hd1[11]/6+hd2[11]/3+hd3[11]/3+hd4[11]/6);

				l[0].mdata[vector_index][0] += ts.dt*(ld1[0]/6+ld2[0]/3+ld3[0]/3+ld4[0]/6);
				l[0].mdata[vector_index][1] += ts.dt*(ld1[1]/6+ld2[1]/3+ld3[1]/3+ld4[1]/6);
				l[1].mdata[vector_index][0] += ts.dt*(ld1[2]/6+ld2[2]/3+ld3[2]/3+ld4[2]/6);
				l[1].mdata[vector_index][1] += ts.dt*(ld1[3]/6+ld2[3]/3+ld3[3]/3+ld4[3]/6);
				l[2].mdata[vector_index][0] += ts.dt*(ld1[4]/6+ld2[4]/3+ld3[4]/3+ld4[4]/6);
				l[2].mdata[vector_index][1] += ts.dt*(ld1[5]/6+ld2[5]/3+ld3[5]/3+ld4[5]/6);
				l[3].mdata[vector_index][0] += ts.dt*(ld1[6]/6+ld2[6]/3+ld3[6]/3+ld4[6]/6);
				l[3].mdata[vector_index][1] += ts.dt*(ld1[7]/6+ld2[7]/3+ld3[7]/3+ld4[7]/6);
				l[4].mdata[vector_index][0] += ts.dt*(ld1[8]/6+ld2[8]/3+ld3[8]/3+ld4[8]/6);
				l[4].mdata[vector_index][1] += ts.dt*(ld1[9]/6+ld2[9]/3+ld3[9]/3+ld4[9]/6);
				l[5].mdata[vector_index][0] += ts.dt*(ld1[10]/6+ld2[10]/3+ld3[10]/3+ld4[10]/6);
				l[5].mdata[vector_index][1] += ts.dt*(ld1[11]/6+ld2[11]/3+ld3[11]/3+ld4[11]/6);

				//R-K routine
			}
		}
	}
    
// file for testing in case you want to output any variables--these are some i found usefull
#if 0    
    static FILE *tomout = NULL;
    
    if (!tomout) tomout = fopen("tomout.txt", "w");
        
    fprintf(tomout, "%10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   %10.10" RFP "f   \n", t, 
	    phi[0][0], 
	    phi[0][1],	  
	    h[2][1][5][6], 
	    h[3][1][5][6],	    
	    h[4][1][5][6], 
	    h[5][1][5][6],
	    h[6][1][5][6], 
	    h[7][1][5][6],	  
	    h[8][1][5][6], 
	    h[9][1][5][6],	    
	    h[10][1][5][6], 
	    h[11][1][5][6],
	    (h[0][1][5][6] + 5*h[2][1][5][6] + 6*h[4][1][5][6])/L,
	    (h[3][1][5][6] + 5*h[7][1][5][6] + 6*h[9][1][5][6])/L,
	    h[0][1][5][6] + h[6][1][5][6] + h[10][1][5][6],
	    T_gw[0][1][5][6] + 5*T_gw[2][1][5][6]
	    + 6*T_gw[4][1][5][6], 
	    T_gw[0][1][5][6] + T_gw[6][1][5][6]
	    + T_gw[10][1][5][6],
	    T_gw[0][1][5][6],
	    ccoef1, a);
	;
   
    fflush(tomout);
    // fclose (tomout);    
#endif    
}


//function for calculating the derivatives
template <typename R>
void gwderivs(model_params<R> &mp, time_state<R> &ts, R h11, 
	R hi11, 
	R h12, 
	R hi12, 
	R h13, 
	R hi13, 
	R h22, 
	R hi22, 
	R h23, 
	R hi23, 
	R h33, 
	R hi33,
	R l11, 
	R li11, 
	R l12, 
	R li12, 
	R l13, 
	R li13, 
	R l22, 
	R li22, 
	R l23, 
	R li23, 
	R l33, 
	R li33,	int ii, int jj, int k,
	R hd[12], R ld[12], 
	R time, R source_gw[12])
{
	R norm = (pow<2>(2.*M_PI)/(pow<2>(mp.len)*pow(ts.a,2.*mp.rescale_s+2)));
	R omega = (8.*M_PI);
	R kk = (pow<2>((R)ii) + pow<2>((R)jj) + pow<2>((R)k))*norm;

	// The result here is in program units except for a missing factor of 1/B^2

	ld[0] = -  kk*h11  + 2.*omega*source_gw[0];
	ld[1] = -  kk*hi11 + 2.*omega*source_gw[1];
	ld[2] = -  kk*h12  + 2.*omega*source_gw[2];
	ld[3] = -  kk*hi12 + 2.*omega*source_gw[3];
	ld[4] = -  kk*h13  + 2.*omega*source_gw[4];
	ld[5] = -  kk*hi13 + 2.*omega*source_gw[5];
	ld[6] = -  kk*h22  + 2.*omega*source_gw[6];
	ld[7] = -  kk*hi22 + 2.*omega*source_gw[7];
	ld[8] = -  kk*h23  + 2.*omega*source_gw[8];
	ld[9] = -  kk*hi23 + 2.*omega*source_gw[9];
	ld[10] = - kk*h33  + 2.*omega*source_gw[10];
	ld[11] = - kk*hi33 + 2.*omega*source_gw[11];

	hd[0]=l11;
	hd[1]=li11;
	hd[2]=l12;
	hd[3]=li12;
	hd[4]=l13;
	hd[5]=li13;
	hd[6]=l22;
	hd[7]=li22;
	hd[8]=l23;
	hd[9]=li23;
	hd[10]=l33;
	hd[11]=li33;

	/* printf("%" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g\n", hd[0], hd[1], hd[2], hd[3], hd[4], hd[5], hd[6], hd[7], hd[8], hd[9], hd[10], hd[11]); */
	/* printf("%" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g %" RFP "g\n", ld[0], ld[1], ld[2], ld[3], ld[4], ld[5], ld[6], ld[7], ld[8], ld[9], ld[10], ld[11]); */
}

template <typename R>
void private_globals<R>::jtg_spectrum(field_size &fs, model_params<R> &mp, time_state<R> &ts)
{
	R array_out[(int)(1.73205*(fs.n/2))+1];
	int numpoints_gw[(int)(1.73205*(fs.n/2))+1]; // fs.number of points in each momentum bin
	R p[(int)(1.73205*(fs.n/2))+1];
	R f2_gw[(int)(1.73205*(fs.n/2))+1]; // Values for each bin: Momentum, |F-k|^2, n_k
	int numbins_gw=(int)(sqrt(3.0)*(fs.n/2))+1; // Actual number of bins for the number of dimensions

	R dp_gw=2.*M_PI/(R)mp.len*mp.rescale_B*6.0e10/ts.a/sqrt(mp.rescale_B*pow(ts.a,mp.rescale_s)*hub); // Size of grid spacing in momentum space
	R fp2_gw;
	R norm1_gw=4.e-5/pow(100.,.333)/pow(ts.a,(R) 2.*mp.rescale_r)*M_PI/3/pow<2>(mp.rescale_A)/pow(mp.len,(R) 3.)/pow<2>(hub);

	R fp2_gwc[6];
	R f2_gwc[6][(int)(1.73205*(fs.n/2))+1]; // Values for each bin: Momentum, |F-k|^2, n_k
	R array_outc[6][(int)(1.73205*(fs.n/2))+1];

	// Calculate magnitude of momentum in each bin
	for(int i=0;i<numbins_gw;i++) {
		p[i]=dp_gw*(R)i;
		f2_gw[i]=0.0;
		numpoints_gw[i]=0;

		for (int j = 0; j < 6; ++j) {
			f2_gwc[j][i] = R(0);
		}
	}

	for (int i = 0; i < 6; ++i) {
		h[i].switch_state(momentum);
		l[i].switch_state(momentum);
	}
    
// #ifdef _OPENMP
// #pragma omp parallel for reduction(+:f2_gw,numpoints_gw) private(fp2_gw)
// #endif
   
	for(int i=0;i<fs.n;i++) {
		int px = (i<=fs.n/2 ? i : i-fs.n);
		for(int j=0;j<fs.n;j++) {
			int py = (j<=fs.n/2 ? j : j-fs.n);
			for(int k=0;k<fs.n/2+1;k++) { 
				int pz = k;
				int vector_index = k + (fs.n/2+1)* (j + fs.n*i);
				R pmagnitude_gw = sqrt(pow<2>(px)+pow<2>(py)+pow<2>(pz));

				fp2_gw = pow<2>(l[0].mdata[vector_index][0]-mp.rescale_r*hub*h[0].mdata[vector_index][0]) 
					+ pow<2>(l[0].mdata[vector_index][1]-mp.rescale_r*hub*h[0].mdata[vector_index][1])
					+ 2.*pow<2>(l[1].mdata[vector_index][0]-mp.rescale_r*hub*h[1].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[1].mdata[vector_index][1]-mp.rescale_r*hub*h[1].mdata[vector_index][1])
					+ 2.*pow<2>(l[2].mdata[vector_index][0]-mp.rescale_r*hub*h[2].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[2].mdata[vector_index][1]-mp.rescale_r*hub*h[2].mdata[vector_index][1])
					+ pow<2>(l[3].mdata[vector_index][0]-mp.rescale_r*hub*h[3].mdata[vector_index][0]) 
					+ pow<2>(l[3].mdata[vector_index][1]-mp.rescale_r*hub*h[3].mdata[vector_index][1])
					+ 2.*pow<2>(l[4].mdata[vector_index][0]-mp.rescale_r*hub*h[4].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[4].mdata[vector_index][1]-mp.rescale_r*hub*h[4].mdata[vector_index][1])
					+ pow<2>(l[5].mdata[vector_index][0]-mp.rescale_r*hub*h[5].mdata[vector_index][0]) 
					+ pow<2>(l[5].mdata[vector_index][1]-mp.rescale_r*hub*h[5].mdata[vector_index][1]);
				int fac = (k == 0 || k == fs.n/2) ? 1 : 2;
				numpoints_gw[(int)pmagnitude_gw] += fac;
				f2_gw[(int)pmagnitude_gw] += fac*fp2_gw;

				fp2_gwc[0] = pow<2>(l[0].mdata[vector_index][0]-mp.rescale_r*hub*h[0].mdata[vector_index][0]) 
					+ pow<2>(l[0].mdata[vector_index][1]-mp.rescale_r*hub*h[0].mdata[vector_index][1]);
				fp2_gwc[1] = 2.*pow<2>(l[1].mdata[vector_index][0]-mp.rescale_r*hub*h[1].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[1].mdata[vector_index][1]-mp.rescale_r*hub*h[1].mdata[vector_index][1]);
				fp2_gwc[2] = 2.*pow<2>(l[2].mdata[vector_index][0]-mp.rescale_r*hub*h[2].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[2].mdata[vector_index][1]-mp.rescale_r*hub*h[2].mdata[vector_index][1]);
				fp2_gwc[3] = pow<2>(l[3].mdata[vector_index][0]-mp.rescale_r*hub*h[3].mdata[vector_index][0]) 
					+ pow<2>(l[3].mdata[vector_index][1]-mp.rescale_r*hub*h[3].mdata[vector_index][1]);
				fp2_gwc[4] = 2.*pow<2>(l[4].mdata[vector_index][0]-mp.rescale_r*hub*h[4].mdata[vector_index][0]) 
					+ 2.*pow<2>(l[4].mdata[vector_index][1]-mp.rescale_r*hub*h[4].mdata[vector_index][1]);
				fp2_gwc[5] = pow<2>(l[5].mdata[vector_index][0]-mp.rescale_r*hub*h[5].mdata[vector_index][0]) 
					+ pow<2>(l[5].mdata[vector_index][1]-mp.rescale_r*hub*h[5].mdata[vector_index][1]);

				for (int dj = 0; dj < 6; ++dj) {
					f2_gwc[dj][(int)pmagnitude_gw] += fac*fp2_gwc[dj];
				}
			}
		}
	}

	// cout << "norm1_gw: " << norm1_gw << endl;

	for(int i=0;i<numbins_gw;i++) {
		if(numpoints_gw[i]>0) {// Converts sums to averages. (numpoints[i] should always be greater than zero.)
			R scl = norm1_gw*pow((R)i/mp.len,(R)3.)/(R)numpoints_gw[i];
			array_out[i] = scl*f2_gw[i];

			for (int j = 0; j < 6; ++j) {
				array_outc[j][i] = scl*f2_gwc[j][i];
			}
		}
		else {
			array_out[i] = 0.;

			for (int j = 0; j < 6; ++j) {
				array_outc[j][i] = 0.;
			}
		}
	}

	if (!ofs.is_open()) {
		ofs.open("gwspec.tsv");
		ofs << setprecision(10);
	}

	for(int i=0;i<numbins_gw;i++) {
		ofs << ts.t << "\t" << p[i] << "\t" << array_out[i] << "\t" << ts.adot/ts.a << "\t" << ts.a;
		for (int j = 0; j < 6; ++j) {
			ofs << "\t" << array_outc[j][i];
		}
		ofs << endl;
	}

	ofs << endl;
	ofs.flush();

	return;
}
 
// Explicit instantiations
template struct private_globals<double>;
#ifdef USE_LD
template struct private_globals<long double>;
#endif

