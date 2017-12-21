
#pragma once
#ifndef __INIT__
#define __INIT__
#include "initializing_inline.h"

template<typename T>
af::array theMask(af::array matrix);
/*
template<typename T>
inline void getSubgrid(int subgrid[], std::pair<T, T> start, std::pair<T, T> end,
	arma::Col<T> xcords, arma::Col<T> ycords);
*/

void setDims(int *x_dim, int *y_dim, int *t_dim);
/*
template<typename T>
void setSubgrid_Delta(std::pair<T, T> *start, std::pair<T, T> *end,
	arma::Col<T> xcords, arma::Col<T> ycords);
	*/

#endif // !__INIT__

