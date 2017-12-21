
#pragma once
#ifndef _PARTICLES_
#define _PARTICLES_


template<typename T>
void initParticles(T delta, af::array *particles_device,
	std::pair<T, T> startxy, std::pair<T, T> endxy);
template<typename T>
void ftle(af::array particles, T dt, T tmax, T delta, int dimx, int dimy, bool backward);

template<typename T>
inline void translateParticle(void func(
	af::array *, af::array *, af::array *, T *, af::array, af::array, int, af_interp_type),
	T * dt, T * tmax, af::array fieldx, af::array fieldy, af::array * coords,
	af::array xcoords, af::array ycoords, af_interp_type interp_t);

template<typename T>
void rk4_integrator(af::array *coords, af::array *fieldx, af::array *fieldy,
	T* dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t);
template<typename T>
void heuns_integrator(af::array *coords, af::array *fieldx, af::array *fieldy,
	T* dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t);


template<typename T>
void heuns_integrator_F_B(
		af::array *coords, af::array *fieldx, af::array *fieldy,
		T* dt, af::array xcoords, af::array ycoords,
		 int t_step, af_interp_type interp_t, bool backward
	);


template<typename T>
void rk4_integrator_F_B(
		af::array *coords, af::array *fieldx, af::array *fieldy,
		T* dt, af::array xcoords, af::array ycoords,
		 int t_step, af_interp_type interp_t, 
		 bool backward
	);
					
#include "part_trans_inline.h"
#endif // !_PARTICLES_