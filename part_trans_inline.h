
#pragma once
#include <arrayfire.h>
#include <exception>
using namespace af;
using namespace std;

template <typename T>
af::array getMaxEigV2x2(af::array matr);
template <typename T>
void getAfType(af::dtype *type);

/**
 * FUNCTION FOR INITIATING PARTICLE POSITIONS
 * PARAMS: T, ARRAY*, PAIR, PAIR
 * 1. DISTANCE BETWEEN PARTICLES(IN KM)
 * 2. ARRAY OF PARTICLES
 * 3. PARAMETERS TO SET STARTING CORNER OF PARTICLE GRID
 * 4. PARAMETERS TO SET ENDING CORNER OF PARTICLE GRID
*/
template <typename T>
inline void initParticles(T delta, af::array *particles_device, std::pair<T, T> startxy, std::pair<T, T> endxy)
{
	cout << "Initiating starting positions..." << endl;
	af::dtype aftype;
	getAfType<T>(&aftype);

	unsigned int part_x = abs(endxy.first - startxy.first) / delta;
	unsigned int part_y = abs(endxy.second - startxy.second) / delta;
	dim_t t_dim = (*particles_device).dims(2);
	//SETTING THE STARTING POSITIONS ON THE HOST-SIDE
	T *host_arrx = (T *)malloc(part_x * part_y * sizeof(T));
	T *host_arry = (T *)malloc(part_x * part_y * sizeof(T));
	for (unsigned int j = 0; j < part_y; j++)
	{
		for (unsigned int i = 0; i < part_x; i++)
		{
			host_arrx[i + j * part_x] = (startxy.first + delta * i);
			host_arry[i + j * part_x] = (startxy.second + delta * j);
		}
	}
	//TRANSFERRING THE STARTING POSITIONS TO THE GPU
	af::array xstart(1, 1, (*particles_device).dims(0), 1, host_arrx);
	af::array ystart(1, 1, (*particles_device).dims(0), 1, host_arry);
	//SET THE STARTING POSITIONS INTO THE PARTICLE ARRAY
	try
	{
		(*particles_device)(span, 0, 0) = xstart;
		(*particles_device)(span, 1, 0) = ystart;
	}
	catch (std::exception &er)
	{
		cout << er.what() << endl;
		free(host_arrx);
		free(host_arry);
		throw "something gone wrong in initParts funk";
	}
	free(host_arrx);
	free(host_arry);
	//THE REST OF THE FUNCTION MAKES BINARY FILES \
	CONTAINING THE COORDINATE GRID FOR THE PARTICLES
	af::array xcord = iota(dim4(part_x, 1, 1, 1));
	xcord = xcord * (T)delta;
	xcord += startxy.first;
	af::array ycord = iota(dim4(part_y, 1, 1, 1));
	ycord = ycord * (T)delta;
	ycord += startxy.second;
	T *host_xcord = xcord.host<T>();
	T *host_ycord = ycord.host<T>();
	ofstream file_xcord("output/" + to_string(part_x - 2) + "_" + to_string(part_y - 2) + "grid_xcords.txt",
						std::ios::out | std::ios::binary);
	if (file_xcord.is_open())
	{
		file_xcord.write(reinterpret_cast<char *>(host_xcord), xcord.dims(0) * sizeof(T));
	}
	file_xcord.close();
	ofstream file_ycord("output/" + to_string(part_x - 2) + "_" + to_string(part_y - 2) + "grid_ycords.txt",
						std::ios::out | std::ios::binary);
	if (file_ycord.is_open())
	{
		file_ycord.write(reinterpret_cast<char *>(host_ycord), ycord.dims(0) * sizeof(T));
	}
	file_ycord.close();
	free(host_ycord);
	free(host_xcord);
}



/**
 * INTEGRATOR HEUNS METHOD
 * FOR FORWARD AND BACKWARD PROPAGATING PARTICLES
 * PARAMS: ARRAY*, ARRAY*, ARRAY*, ARRAY*, T*, ARRAY,ARRAY, INT, AF_INTERP_TYPE, BOOL
 * 1. PARTICLE ARRAY
 * 2. & 3. VELOCITY FIELD X,Y COMPONENTS
 * 4. dt
 * 5. & 6. SHORT ARRAY OF X- & Y-COORDINATES. ONLY USES THE 2 FIRST VALUES
 * 6. STRIDE USED BY THE INTEGRATOR
 * 7. INTERPOLATION TYPE TO USE
 * 8. BACKWARD PROPAGATION BOOL (1 = BACKWARD PROP.)
*/
template <typename T>
void heuns_integrator_F_B(af::array *coords, af::array *fieldx, af::array *fieldy,
						  T *dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t, bool backward)
{
	try
	{
		//THE FIRST FEW LINES CALCULATES THE RESOLUTION OF THE FIELD (YES I SHOULD HAVE USED GLOBALS MORE)
		T *xcord_host = xcoords.host<T>();
		T *ycord_host = ycoords.host<T>();

		T resolutionx = xcord_host[1] - xcord_host[0];
		T resolutiony = ycord_host[1] - ycord_host[0];

		dim_t steps = fieldx->dims(2);
		dim_t N_parts = coords->dims(0);
		N_parts = (unsigned int)N_parts;
		int N_saves = coords->dims(2);

		///////////////
		///////////////

		af::dtype datatype;
		getAfType<T>(&datatype);

		af::dim4 dimens(N_parts, 1, 1, 1); //dims used when working with interpolation objects

		//THIS ARRAY MIGHT BE A BIT FASTER \
		BUT COULD PROBABLY JUST USE dt DIRECTLY
		af::array dt_d = constant(*dt, dimens, datatype);

		//THE INTERPOLATION FUNCTION TAKES COORDINATES 0,1,2,...,X_DIM , 0,1,2,...,Y_DIM \
		WHICH MEANS THE COORDINATES HAD TO BE RENORMALIZED TO FIT THIS SPAN.
		af::array resx = constant(resolutionx, dimens, datatype);
		af::array resy = constant(resolutiony, dimens, datatype);
		af::array offsetx = constant(xcord_host[0], dimens, datatype);
		af::array offsety = constant(ycord_host[0], dimens, datatype);
		
		//ARRAYS FOR STORING COORDINATES
		af::array x(dimens, datatype), y(dimens, datatype);

		//FIXING INCREMENTATION IF THERE SHOULD BE A STRIDE \
		OR BACKWARD PROPAGATION.
		int index = (steps - 1) * backward;
		int save_iter = 1;
		int iter_correction = ((steps - 1) / N_saves) == 0 ? 0 : 1;
		int increment = 1 + backward * (-2);
		cout << "Moving " << N_parts << " particles through the field." << endl;
		x = (*coords)(span, 0, 0);
		y = (*coords)(span, 1, 0);
		do
		{
			try
			{
				//CREATING ARRAYS FOR THE K'S
				af::array kx1(dimens, datatype), kx2(dimens, datatype);
				af::array ky1(dimens, datatype), ky2(dimens, datatype);
				//APPROX2 IS AN INTERPOLATION FUNCTION \
				TAKES A 2D FIELD, X- Y-COORDS, INTERPOL OBJECT AND \
				A RETURN VALUE SHOULD X OR Y BE OUTSIDE THE FIELD.
				kx1 = approx2((*fieldx)(span, span, index, 0), (x - offsetx) / resx,
							  (y - offsety) / resy, interp_t, 0.0f);
				ky1 = approx2((*fieldy)(span, span, index, 0), (x - offsetx) / resx,
							  (y - offsety) / resy, interp_t, 0.0f);
				{
					//CREATING ARRAY XN,YN
					af::array xn(dimens, datatype), yn(dimens, datatype);
					xn = x + (T)increment * (dt_d)*kx1;
					yn = y + (T)increment * (dt_d)*ky1;

					kx2 = approx2((*fieldx)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					ky2 = approx2((*fieldy)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
				}

				x = x + (T)increment * (dt_d) * (kx1 + kx2) / (T)(2.0);
				eval(x); //CALL EVAL TO FORCE CALCULATION OF X. MIGH BE BETTER TO CALL \
				BOTH EVALUATIONS IN THE SAME FUNCTION -- EVAL(X,Y)
				y = y + (T)increment * (dt_d) * (ky1 + ky2) / (T)(2.0);
				eval(y);

				if ((index + 1) % (((steps/t_step) / N_saves) + iter_correction) == 0)
				{
					//WANT TO STORE POSITION AT THE END OF AN INTERVAL. \
					THEN RESET GRID.- COMMENT OUT RESET PART IF THE WHOLE TRAJECTORY IS WANTED
					(*coords)(span, 0, save_iter) = x;
					(*coords)(span, 1, save_iter++) = y;
					x = (*coords)(span, 0, 0);
					y = (*coords)(span, 1, 0);
				}
			}
			catch (std::exception &error)
			{
				free(xcord_host);
				free(ycord_host);
				throw error;
			}
			index += increment * t_step;
		} while (index < (steps - 1*t_step) && index > 0);
		(*coords)(span, 0, N_saves - 1) = x;
		(*coords)(span, 1, N_saves - 1) = y;
		free(xcord_host);
		free(ycord_host);
		//CODE FOR WRITING OUT THE ARRAY OF PARTICLES FOR PLOTTING
		if (!backward)
		{
			T *particle_endpoints = (*coords)(span, span, span).host<T>();
			ofstream particle_ends(
				"output/_particle_ends_N" + to_string(N_saves - 1) + ".txt",
				std::ios::out | std::ios::binary);
			if (particle_ends.is_open())
			{
				particle_ends.write(
					reinterpret_cast<char *>(particle_endpoints),
					(*coords).dims(0) * sizeof(T) * 2 * (N_saves));
				//cout << (dimx) * (dimy) * 2<<endl;
			}
			particle_ends.close();
			free(particle_endpoints);
		}
	}
	catch (af::exception &err)
	{
		throw err;
	}
}
/**
 * INTEGRATOR RK4
 * FOR FORWARD AND BACKWARD PROPAGATING PARTICLES
 * PARAMS: ARRAY*, ARRAY*, ARRAY*, ARRAY*, T*, ARRAY,ARRAY, INT, AF_INTERP_TYPE, BOOL
 * 1. PARTICLE ARRAY
 * 2. & 3. VELOCITY FIELD X,Y COMPONENTS
 * 4. dt
 * 5. & 6. SHORT ARRAY OF X- & Y-COORDINATES. ONLY USES THE 2 FIRST VALUES
 * 6. STRIDE USED BY THE INTEGRATOR
 * 7. INTERPOLATION TYPE TO USE
 * 8. BACKWARD PROPAGATION BOOL (1 = BACKWARD PROP.)
*/
template <typename T>
void rk4_integrator_F_B(af::array *coords, af::array *fieldx, af::array *fieldy,
						T *dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t, bool backward)
{
	//THE FIRST FEW LINES CALCULATES THE RESOLUTION OF THE FIELD (YES I SHOULD HAVE USED GLOBALS MORE)
	T *xcord_host = xcoords.host<T>();
	T *ycord_host = ycoords.host<T>();
	T resolutionx = xcord_host[1] - xcord_host[0];
	T resolutiony = ycord_host[1] - ycord_host[0];
	dim_t steps = fieldx->dims(2);

	dim_t N_parts = coords->dims(0);
	N_parts = (unsigned int)N_parts;
	int N_saves = coords->dims(2);
	///////////////
	///////////////
	cout << "Moving " << N_parts << " particles through the field(RK4)." << endl;

	af::dtype datatype;
	getAfType<T>(&datatype);
	af::dim4 dimens(N_parts, 1, 1, 1); 

	//THIS ARRAY MIGHT BE A BIT FASTER \
		BUT COULD PROBABLY JUST USE dt DIRECTLY
	af::array dt_d = constant(*dt, dimens, datatype);

	//THE INTERPOLATION FUNCTION TAKES COORDINATES 0,1,2,...,X_DIM , 0,1,2,...,Y_DIM \
		WHICH MEANS THE COORDINATES HAD TO BE RENORMALIZED TO FIT THIS SPAN.
	af::array resx = constant(resolutionx, dimens, datatype);
	af::array resy = constant(resolutiony, dimens, datatype);
	af::array offsetx = constant(xcord_host[0], dimens, datatype);
	af::array offsety = constant(ycord_host[0], dimens, datatype);

	//ARRAYS FOR STORING COORDINATES
	af::array x(dimens, datatype), y(dimens, datatype);

	//FIXING INCREMENTATION IF THERE SHOULD BE A STRIDE \
		OR BACKWARD PROPAGATION.
	int index = (steps - 1) * backward;
	int iter = 0;
	int save_iter = 1;
	int iter_correction = ((steps - 1) / N_saves) == 0 ? 0 : 1;

	int increment = 1 + backward * (-2);

	x = (*coords)(span, 0, 0);
	y = (*coords)(span, 1, 0);
	do
	{

		try
		{
			//CREATING ARRAYS FOR THE K'S
			af::array kx1(dimens, datatype), kx2(dimens, datatype), kx3(dimens, datatype), kx4(dimens, datatype);
			af::array ky1(dimens, datatype), ky2(dimens, datatype), ky3(dimens, datatype), ky4(dimens, datatype);
			//CREATING ARRAY XN,YN
			af::array xn(dimens, datatype), yn(dimens, datatype);
			//APPROX2 IS AN INTERPOLATION FUNCTION \
			TAKES A 2D FIELD, X- Y-COORDS, INTERPOL OBJECT AND \
			A RETURN VALUE SHOULD X OR Y BE OUTSIDE THE FIELD.
			kx1 = approx2((*fieldx)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);
			ky1 = approx2((*fieldy)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);
			{

				xn = x + (T)increment * (dt_d)*kx1 * (T)0.5;
				yn = y + (T)increment * (dt_d)*ky1 * (T)0.5;
				//TWO CASES:
				//1. IF THE NUMBER OF STEPS EQUAL THE DIMENSION OF THE AVAILABLE DATA \
				 		RESULTING IN LIN INTERPOLATION IN TIME\
				  2. NO INTERPOLATION - STRIDED USE OF DATA.
				if (t_step == 1)
				{

					kx2 = (approx2((*fieldx)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldx)(span, span, index + increment, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;

					ky2 = (approx2((*fieldy)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldy)(span, span, index + increment, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;

					xn = x + (T)increment * (dt_d)*kx2 * (T)0.5;
					yn = y + (T)increment * (dt_d)*ky2 * (T)0.5;

					kx3 = (approx2((*fieldx)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldx)(span, span, index + increment, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;
					ky3 = (approx2((*fieldy)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldy)(span, span, index + increment, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;
				}
				else
				{
					kx2 = approx2((*fieldx)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					ky2 = approx2((*fieldy)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					xn = x + (T)increment * (dt_d)*kx2 * 0.5;
					yn = y + (T)increment * (dt_d)*ky2 * 0.5;
					kx3 = approx2((*fieldx)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					ky3 = approx2((*fieldy)(span, span, index + increment, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
				}
				xn = x + (T)increment * (dt_d)*kx3;
				yn = y + (T)increment * (dt_d)*ky3;

				kx4 = approx2((*fieldx)(span, span, index + increment * t_step, 0), (xn - offsetx) / resx,
							  (yn - offsety) / resy, interp_t, 0.0f);
				ky4 = approx2((*fieldy)(span, span, index + increment * t_step, 0), (xn - offsetx) / resx,
							  (yn - offsety) / resy, interp_t, 0.0f);
			}

			x = x + (T)increment * (dt_d) * (kx1 + (T)2.0 * kx2 + 2.0 * kx3 + kx4) / (T)(6.0);
			eval(x);
			//CALL EVAL TO FORCE CALCULATION OF X. MIGH BE BETTER TO CALL \
				BOTH EVALUATIONS IN THE SAME FUNCTION -- EVAL(X,Y)
			y = y + (T)increment * (dt_d) * (ky1 + (T)2.0 * ky2 + (T)2.0 * ky3 + ky4) / (T)(6.0);
			eval(y);
		
			if ((iter + 1) % (((steps / t_step) / (N_saves)) + iter_correction) == 0)
				{
				//WANT TO STORE POSITION AT THE END OF AN INTERVAL. \
					THEN RESET GRID. - COMMENT OUT RESET PART IF THE WHOLE TRAJECTORY IS WANTED
					(*coords)(span, 0, save_iter) = x;
					(*coords)(span, 1, save_iter++) = y;
					x = (*coords)(span, 0, 0);
					y = (*coords)(span, 1, 0);
				}
		}
		catch (af::exception &error)
		{
			cout << error.what() << endl;
		}
		index += increment * t_step;
		iter++;
	} while (index < (steps - 1*t_step) && index > 0);
	free(xcord_host);
	free(ycord_host);
	(*coords)(span, 0, N_saves - 1) = x;
	(*coords)(span, 1, N_saves - 1) = y;
}
/**
 * FUNCTION FOR CALCULATING THE FTLE FIELD
 * PARAMS: ARRAY, T, T, T, INT, INT, BOOL
 * 1. PARTICLE ARRAY
 * 2. DT
 * 3. TMAX
 * 4. DELTA
 * 5. NUM PARTICLES X
 * 6. NUM PARTICLES Y
 * 7. BACKWARD PROPAGATION
*/
template <typename T>
void ftle(af::array particles, T dt, T tmax, T delta, int dimx, int dimy, bool backward)
{
	af::dtype datatype;
	getAfType<T>(&datatype);
	cout << "Calculating Lyapunov exponents..." << endl;

	int iter_max = particles.dims(2) - 1;
	int tot_parts = dimx * dimy;

	//FINDING 2 TIMES THE DISTANCE BETWEEN PARTICLES IN X,Y DIR
	T * x_host = particles(seq(0,3), 0, 0).host<T>();
	T * y_host = particles(seq(0,3*dimx,dimx), 1, 0).host<T>();
	T delta_x = x_host[2] - x_host[0];
	T delta_y = y_host[2] - y_host[0];

	//DIMENSIONS FOR THE FTLE-FIELD AND THE STRAIN-TENSOR
	dim4 dr_dim(2, 2, tot_parts);
	af::array ftle_dev(tot_parts, iter_max, datatype);

	string filename;
	if (backward)
	{
		filename = "output/" + to_string(dimx - 2) + "_" + to_string(dimy - 2) +
				   "_ftle_heat_t5_delta" + to_string(delta) +
				   "_dt" + to_string(dt) + "BACK_N" + to_string(iter_max) + ".txt";
	}
	else
	{
		filename = "output/" + to_string(dimx - 2) + "_" + to_string(dimy - 2) +
				   "_ftle_heat_t5_delta" + to_string(delta) +
				   "_dt" + to_string(dt) + "FORW_N" + to_string(iter_max) + ".txt";
	}

	try
	{	//LOOPING OVER THE INTERVALS
		for (int i = 1; i < iter_max + 1; i++)
		{
			//CALCULATING THE STRAINTENSOR
			af::array dr(dr_dim, datatype);
			{
				af::array rx(tot_parts, datatype), ry(tot_parts, datatype);

				rx = shift(particles(span, 0, i), 1, 0, 0) - shift(particles(span, 0, i), -1, 0, 0);

				dr(0, 0, span) = rx  / delta_x;
				rx = shift(particles(span, 0, i), dimx, 0, 0) - shift(particles(span, 0, i), -dimx, 0, 0);

				dr(0, 1, span) = rx / delta_y;
				ry = shift(particles(span, 1, i), 1, 0, 0) - shift(particles(span, 1, i), -1, 0, 0);

				dr(1, 0, span) = ry  / delta_x;
				ry = shift(particles(span, 1, i), dimx, 0, 0) - shift(particles(span, 1, i), -dimx, 0, 0);

				dr(1, 1, span) = ry  / delta_y;
			}
			af::array dr_sq(dr_dim, datatype);
			eval(dr);
			{
				//MAT-MULT	
				dr_sq(0, 0, span) = dr(0, 0, span) * dr(0, 0, span) + dr(1, 0, span) * dr(1, 0, span);
				dr_sq(0, 1, span) = dr(0, 0, span) * dr(0, 1, span) + dr(1, 0, span) * dr(1, 1, span);
				dr_sq(1, 0, span) = dr(0, 0, span) * dr(0, 1, span) + dr(1, 0, span) * dr(1, 1, span);
				dr_sq(1, 1, span) = dr(0, 1, span) * dr(0, 1, span) + dr(1, 1, span) * dr(1, 1, span);
			}
			eval(dr_sq);
			//FINDING EIGENVALUES AND CALCULATING FTLE-FIELD
			af::array temp = log(sqrt(getMaxEigV2x2<T>(dr_sq))) / (T)1.0; 
			eval(temp);
			//RESHAPE THE 1D ARRAY TO A 2D GRID
			ftle_dev(span, i - 1) = moddims(temp, tot_parts);
		}

		//THERE MIGHT BE SOME VALUES WHICH RESULT IN A ZERO IN THE LOG ABOVE
		//THIS LINE TAKES THEM OUT AND SETS THEM TO ZERO
		//THIS IS LIKELY CASED BY INCLUDING PARTICLES ON LAND.
		//////////////
		af::replace(ftle_dev, !((isNaN(ftle_dev) || isInf(ftle_dev))), 0.0);
		dim4 dims(dimx, dimy, ftle_dev.dims(1), 1);
		ftle_dev = moddims(ftle_dev, dims);
		//TRANSFER ONLY THE INNER PARTICLES BACK TO HOST
		T *host_ftle = ftle_dev(seq(1, dimx - 2), seq(1, dimy - 2), span).host<T>();
		cout << "Saving results to file: " << filename << endl;
		ofstream file_ftle(filename, std::ios::out | std::ios::binary);
		if (file_ftle.is_open())
		{
			file_ftle.write(reinterpret_cast<char *>(host_ftle), (dimx - 2) * (dimy - 2) * sizeof(T) * iter_max);
		}
		file_ftle.close();
		free(host_ftle);
	}
	catch (af::exception &err)
	{
		cout << err.what() << endl;
		cout << "Error in calculating the FTLE field..." << endl;
	}
}

//FUNCTION FOR FINDING THE LARGEST EIGENVALUE \
OF A 2X2 POSITIVE DEFINITE MATRIX
template <typename T>
af::array getMaxEigV2x2(af::array matr)
{
	dim4 dimens(matr.dims(2), 1, 1);
	af::array temp1(dimens);

	af::array a(dimens), b(dimens), c(dimens), d(dimens);
	a = (matr)(0, 0, span);
	b = (matr)(0, 1, span);
	c = (matr)(1, 0, span);
	d = (matr)(1, 1, span);

	temp1 = sqrt((a + d) * (a + d) - (T)4.0 * (a * d - b * c));

	temp1 = (a + d + temp1) / (T)2.0;

	eval(temp1);
	return temp1;
}
//METHOD FOR SETTING THE ARRAYFIRE DATATYPE
template <typename T>
void getAfType(af::dtype *type)
{
	if (sizeof(T) == sizeof(float))
	{
		*type = f32;
	}
	else
	{
		*type = f64;
	}
}


//////////////////////////////
//////////////////////////////
//////////////////////////////
/***
 * THE FUNCTIONS BELOW ARE ONLY USED WHEN RUNNING TESTS ON THE DOUBLE-GYRE FIELD
 * */
//////////////////////////////
//////////////////////////////
//////////////////////////////
template <typename T>
inline void translateParticle(void func(
								  af::array *, af::array *, af::array *, T *, af::array, af::array, int, af_interp_type),
							  T *dt, T *tmax, af::array fieldx, af::array fieldy, af::array *coords,
							  af::array xcoords, af::array ycoords, af_interp_type interp_t)
{
	//
	int index = 0;
	dim_t steps = coords->dims(2);
	std::string filename = "output/trajectory_dt_" + std::to_string(*dt) + "_.txt";

	dim_t N_parts = coords->dims(0);
	N_parts = (unsigned int)N_parts;
	cout << "Moving " << N_parts << " particles through the field." << endl;

	///////////////
	///////////////
	//af_interp_type interper = AF_INTERP_CUBIC_SPLINE;
	int t_steps = fieldx.dims(2) / steps; //(*dt >= 3600) ? (*dt) / 3600 : 1;

	(*func)(coords, &fieldx, &fieldy, dt, xcoords, ycoords, t_steps, interp_t);

	cout << "Finished moving particles." << endl;
}

template <typename T>
void rk4_integrator(af::array *coords, af::array *fieldx, af::array *fieldy,
					T *dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t)
{
	T *xcord_host = xcoords.host<T>();
	T *ycord_host = ycoords.host<T>();
	T resolutionx = xcord_host[1] - xcord_host[0];
	T resolutiony = ycord_host[1] - ycord_host[0];
	dim_t steps = coords->dims(2);

	dim_t N_parts = coords->dims(0);
	N_parts = (unsigned int)N_parts;
	///////////////
	///////////////

	af::dtype datatype;
	getAfType<T>(&datatype);
	af::dim4 dimens(N_parts, 1, 1, 1); //dims used when working with interpolation objects
	//af::dim4 dim_retur(1, 1, N_parts, 1); //dims used when putting things back into main array
	//af_interp_type interp_t = AF_INTERP_CUBIC_SPLINE;
	af::array dt_d = constant(*dt, dimens, datatype);

	af::array resx = constant(resolutionx, dimens, datatype);
	af::array resy = constant(resolutiony, dimens, datatype);
	af::array offsetx = constant(xcord_host[0], dimens, datatype);
	af::array offsety = constant(ycord_host[0], dimens, datatype);

	int index = 0;
	int index_cord = 0;
	do
	{
		af::array kx1(dimens, datatype), kx2(dimens, datatype), kx3(dimens, datatype), kx4(dimens, datatype);
		af::array ky1(dimens, datatype), ky2(dimens, datatype), ky3(dimens, datatype), ky4(dimens, datatype);
		af::array x(dimens, datatype), y(dimens, datatype);
		af::array x1(dimens, datatype), y1(dimens, datatype);
		eval(*coords);
		x = (*coords)(span, 0, index_cord);
		y = (*coords)(span, 1, index_cord);

		try
		{

			kx1 = approx2((*fieldx)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);
			ky1 = approx2((*fieldy)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);

			{
				af::array xn(dimens, datatype), yn(dimens, datatype);
				xn = x + (dt_d)*kx1 * 0.5;
				yn = y + (dt_d)*ky1 * 0.5;

				if (t_step == 1)
				{
					kx2 = (approx2((*fieldx)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldx)(span, span, index + 1, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;

					ky2 = (approx2((*fieldy)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldy)(span, span, index + 1, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;

					xn = x + (dt_d)*kx2 * 0.5;
					yn = y + (dt_d)*ky2 * 0.5;

					kx3 = (approx2((*fieldx)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldx)(span, span, index + 1, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;
					ky3 = (approx2((*fieldy)(span, span, index, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f) +
						   approx2((*fieldy)(span, span, index + 1, 0), (xn - offsetx) / resx,
								   (yn - offsety) / resy, interp_t, 0.0f)) /
						  2.0;
				}
				else
				{
					kx2 = approx2((*fieldx)(span, span, index + t_step / 2, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					ky2 = approx2((*fieldy)(span, span, index + t_step / 2, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					xn = x + (dt_d)*kx2 * 0.5;
					yn = y + (dt_d)*ky2 * 0.5;
					kx3 = approx2((*fieldx)(span, span, index + t_step / 2, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
					ky3 = approx2((*fieldy)(span, span, index + t_step / 2, 0), (xn - offsetx) / resx,
								  (yn - offsety) / resy, interp_t, 0.0f);
				}
				xn = x + (dt_d)*kx3;
				yn = y + (dt_d)*ky3;

				kx4 = approx2((*fieldx)(span, span, index + 1 * t_step, 0), (xn - offsetx) / resx,
							  (yn - offsety) / resy, interp_t, 0.0f);
				ky4 = approx2((*fieldy)(span, span, index + 1 * t_step, 0), (xn - offsetx) / resx,
							  (yn - offsety) / resy, interp_t, 0.0f);
			}

			x1 = x + (dt_d) * (kx1 + 2. * kx2 + 2. * kx3 + kx4) / (6.0);

			y1 = y + (dt_d) * (ky1 + 2. * ky2 + 2. * ky3 + ky4) / (6.0);
			//af_print(y)
			//(*coords)(index / t_step + 1, 0, span) = moddims(x1, dim_retur);
			//(*coords)(index / t_step + 1, 1, span) = moddims(y1, dim_retur);
			eval(x1);
			eval(y1);
			(*coords)(span, 0, index_cord + 1) = x1;
			(*coords)(span, 1, index_cord + 1) = y1;
		}
		catch (af::exception &error)
		{
			cout << error.what() << endl;
		}
		index_cord++;
		index = index + 1 * t_step;
	} while (index_cord < steps - 1);
	free(xcord_host);
	free(ycord_host);
	//free(xcord_host); free(ycord_host);
}

template <typename T>
void heuns_integrator(af::array *coords, af::array *fieldx, af::array *fieldy,
					  T *dt, af::array xcoords, af::array ycoords, int t_step, af_interp_type interp_t)
{
	T *xcord_host = xcoords.host<T>();
	T *ycord_host = ycoords.host<T>();
	T resolutionx = xcord_host[1] - xcord_host[0];
	T resolutiony = ycord_host[1] - ycord_host[0];
	dim_t steps = coords->dims(2);

	dim_t N_parts = coords->dims(0);
	N_parts = (unsigned int)N_parts;
	///////////////
	///////////////

	af::dtype datatype;
	getAfType<T>(&datatype);
	af::dim4 dimens(N_parts, 1, 1, 1);	//dims used when working with interpolation objects
	af::dim4 dim_retur(1, 1, N_parts, 1); //dims used when putting things back into main array
										  //af_interp_type interp_t = AF_INTERP_CUBIC_SPLINE;
	af::array dt_d = constant(*dt, dimens, datatype);

	af::array resx = constant(resolutionx, dimens, datatype);
	af::array resy = constant(resolutiony, dimens, datatype);
	af::array offsetx = constant(xcord_host[0], dimens, datatype);
	af::array offsety = constant(ycord_host[0], dimens, datatype);

	int index = 0;
	do
	{
		af::array kx1(dimens, datatype), kx2(dimens, datatype);
		af::array ky1(dimens, datatype), ky2(dimens, datatype);
		af::array xn(dimens, datatype), yn(dimens, datatype), x(dimens, datatype), y(dimens, datatype);
		af::array x1(dimens, datatype), y1(dimens, datatype);
		eval(*coords);
		x = (*coords)(span, 0, index / t_step);
		y = (*coords)(span, 1, index / t_step);
		//x = moddims(x, dimens);
		//y = moddims(y, dimens);
		try
		{

			kx1 = approx2((*fieldx)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);
			ky1 = approx2((*fieldy)(span, span, index, 0), (x - offsetx) / resx,
						  (y - offsety) / resy, interp_t, 0.0f);

			xn = x + (dt_d)*kx1;
			yn = y + (dt_d)*ky1;

			kx2 = approx2((*fieldx)(span, span, index + 1, 0), (xn - offsetx) / resx,
						  (yn - offsety) / resy, interp_t, 0.0f);
			ky2 = approx2((*fieldy)(span, span, index + 1, 0), (xn - offsetx) / resx,
						  (yn - offsety) / resy, interp_t, 0.0f);

			x1 = x + (dt_d) * (kx1 + kx2) / (2.0);
			y1 = y + (dt_d) * (ky1 + ky2) / (2.0);
			//af_print(y)
			eval(x1);
			eval(y1);
			(*coords)(span, 0, index / t_step + 1) = x1;
			(*coords)(span, 1, index / t_step + 1) = y1;
		}
		catch (std::exception &error)
		{
			cout << error.what() << endl;
		}
		index = index + 1 * t_step;
	} while (index < (steps - 1) * t_step);
	free(xcord_host);
	free(ycord_host);
	//free(xcord_host); free(ycord_host);
}
