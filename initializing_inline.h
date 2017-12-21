#pragma once

//#include <armadillo>
#include <iostream>
using namespace std;
/*
template<typename T>
inline void getSubgrid(int subgrid[], std::pair<T, T> start, std::pair<T, T> end, arma::Col<T> xcords, arma::Col<T> ycords) {
	/*
	TO SAVE MEMORY, IF A SMALLER GRID IS CHOSEN, ONLY
	A SMALELR SUBGRID WILL BE CARRIED THROUGH THE PROGRAM.
	THIS WILL HOPEFULLY LESSEN THE BURDEN ON THE GPU LATER

	NOTE: ALSO ADDING EXTRA BUFFER AROUND THE AREA USED FOR COMPUTATIONS.
	THIS WILL MAKE THE PROGRAM MORE RELIABLE IF PARTICLES SHOULD LEAVE
	THE CHOSEN BOX.
	

	arma::Col<arma::uword> hank = (arma::find(abs(xcords - start.first) < xcords(1) - xcords(0)));
	subgrid[0] = (int)hank(0) - 4;
	hank = (arma::find(abs(xcords - end.first) < xcords(1) - xcords(0)));
	subgrid[1] = (int)hank(hank.n_elem - 1) + 4;
	hank = (arma::find(abs(ycords - start.second) < ycords(1) - ycords(0)));
	subgrid[2] = (int)hank(0) - 4;
	hank = (arma::find(abs(ycords - end.second) < ycords(1) - ycords(0)));
	subgrid[3] = (int)hank(hank.n_elem - 1) + 2;

	/*
	CANT FIND ANY FUNCTION IN ARRAYFIRE TO REPLACE THE SEARCH
	MUST BE DONE WITH FOR LOOPS I GUESS
	
}
*/


template<typename T>
inline af::array theMask(af::array matrix) {
	int sx = matrix.dims(0), sy = matrix.dims(1);
	af::array outmat = af::constant(0., sx, sy);
	int slice = matrix.dims(2);
	for (int i = 0; i < sx; i++)
	{
		for (int j = 0; j < sy; j++)
		{
			outmat(i, j) = ((matrix(i, j, 0) == 0) && (matrix(i, j, slice / 2) == 0));
		}
	}
	int * outmat_h = outmat.host<int>();
	ofstream file_mask("output/mask.txt", std::ios::out | std::ios::binary);
	if (file_mask.is_open()) {
	file_mask.write(reinterpret_cast<char *>(outmat_h), sx * sy * sizeof(int));
	}
	file_mask.close();
	return outmat;
	/*
	can print out the mask to check if the data
	has been parsed and loaded correctly.
	plot should replicate part of norwegian coast.
	*/
}

void setDims(int *x_dim, int *y_dim, int *t_dim) {
	cout << "Enter the number of grid points in the x-, y-direction: (must match name of file)" << endl;
	//cin >> x_dim >> y_dim;
	cout << "Enter the number of time-steps: " << endl;
	//cin >> t_dim;
}
/*
template<typename T>
inline void setSubgrid_Delta(std::pair<T, T> *start, std::pair<T, T> *end, arma::Col<T> xcords, arma::Col<T> ycords) {
	cout << "Place particles in a box within range(km) x: " << xcords(0) << " : " << xcords(xcords.n_elem - 1) << endl
		<< ", y: " << ycords(0) << " : " << ycords(ycords.n_elem - 1) << endl;
	cout << "format: x0 x1 y0 y1" << endl;
	//cin >> start->first >> end->first >> start->second >> end->second;
	cout << "Enter distance between particles: " << endl;
	//cin >> delta;
}
*/