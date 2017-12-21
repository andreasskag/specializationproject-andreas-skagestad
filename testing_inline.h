#pragma once

#include <arrayfire.h>
#include "part_trans.h"
using namespace af;
using namespace std;

#define PI 3.141592653

template <typename T>
inline void reference_analytic_rk4(af::array *particles, T dt, T tmax, T A, T omega, T epsilon);
template <typename T>
inline void getLength(af::array *coords, af::array *lengths);

template <typename T>
inline af::array vel_x(T A, T epsilon, T omega, af::array x, af::array y, af::array t);
template <typename T>
inline af::array vel_y(T A, T epsilon, T omega, af::array x, af::array y, af::array t);

template <typename T>
inline void vel_x(T A, T epsilon, T omega, af::array x, af::array y, af::array t, af::array *vel_out);
template <typename T>
inline void vel_y(T A, T epsilon, T omega, af::array x, af::array y, af::array t, af::array *vel_out);

template <typename T>
inline af::array funk(T epsilon, T omega, af::array x, af::array t);
template <typename T>
inline af::array dfunk(T epsilon, T omega, af::array x, af::array t);
template <typename T>
inline void fillFields(af::array *fieldx, af::array *fieldy, af::array x, af::array y, af::array t,
					   T A, T epsilon, T omega);
template <typename T>
inline void errorCalc(af::array ref, af::array approx, af::array *error, af::array ref_len);

template <typename T>
int saveToFile(const char *filename, af::array out);

template<typename T>
af::array findLine(const af::array lgError, const af::array lgTimes);

template <typename T>
inline void testing()
{
	af_dtype aftype;
	getAfType<T>(&aftype);
	T dt_ref;
	T tmax;
	T x_border, y_border, resolution;

	dt_ref = 0.01;
	tmax = 20.0;
	x_border = 2.0;
	y_border = 1.0;
	T delta = 0.005;
	std::pair<T, T> start;
	start.first = x_border / 20.0;
	start.second = y_border / 20.0;
	std::pair<T, T> end;
	end.first = x_border - x_border / 20.0, end.second = y_border - y_border / 20.0;
	/*
	TESTING WITH 3 DIFFERENT PARAMETERS FOR RESOLUTION
	4 FOR DT AND THEN 3 DIFFERENT INTERPOLATION METHODS. 
	WANT TO SAVE MAX ERROR, VARIANCE/STDEV AND POSSIBLY SCALING FACTOR
	COULD BE USEFUL TO SEE PLOT OF ERROR TO KNOW IF THERE IS ANY POINT OF
	TRYING TO FIT A LINE SOMEWHERE.
	*/
	const int t_steps = 4;
	const int x_steps = 3;
	const int interpol_steps = 3;
	T dx_arr[x_steps];
	dx_arr[0] = 0.1;
	dx_arr[1] = 0.05;
	dx_arr[2] = 0.01;
	T dt_arr[t_steps];
	dt_arr[0] = 0.2, dt_arr[1] = 0.1;
	dt_arr[2] = 0.05;
	dt_arr[3] = 0.01;
	af_interp_type interp_type_arr[interpol_steps];
	interp_type_arr[0] = AF_INTERP_BILINEAR;
	interp_type_arr[1] = AF_INTERP_BICUBIC;
	interp_type_arr[2] = AF_INTERP_BICUBIC_SPLINE;

	unsigned int part_x = abs(end.first - start.first) / delta;
	unsigned int part_y = abs(end.second - start.second) / delta;

	int particles = part_x * part_y;
	T A = 0.1, omega = 1.0, epsilon = 0.25;
	af::array particle_ref(particles, 2, tmax / dt_ref, aftype);
	af::array ref_length = constant(0, dim4(particles, particle_ref.dims(2)), aftype);

	initParticles(delta, &particle_ref, start, end);
	reference_analytic_rk4(&particle_ref, dt_ref, tmax, A, omega, epsilon);
	getLength<T>(&particle_ref, &ref_length);
	int res_size_x = t_steps * x_steps * interpol_steps + 2 * t_steps * x_steps;
	int res_size_y = 5;
	af::array results(res_size_x, res_size_y, aftype);

	for (int i = 3; i < interpol_steps; i++)
	{
		af_interp_type interp = interp_type_arr[i];
		for (int xx = 0; xx < x_steps; xx++)
		{
			try
			{
				T res = dx_arr[xx];

				for (int tt = 0; tt < t_steps; tt++)
				{
					T dt = dt_arr[tt];

					/*
					FIRST WE NEED TO MAKE OURSELF SOME DISCRETE DATA
					*/

					int x_len = x_border / res, y_len = y_border / res, n_t = tmax / dt;

					af::array x_cords = seq(0, x_border - res * 0.1, res);
					af::array y_cords = seq(0, y_border - res * 0.1, res);

					af::array particles_arr(particles, 2, n_t, aftype);
					
					initParticles(delta, &particles_arr, start, end);
					af::array vel_field_x(x_len, y_len, n_t, aftype);
					af::array vel_field_y(x_len, y_len, n_t, aftype);
					af::array times = seq(0, tmax, dt);
					fillFields<T>(&vel_field_x, &vel_field_y, x_cords, y_cords, times, A, epsilon, omega);
					translateParticle<T>(&rk4_integrator, &dt, &tmax, vel_field_x, vel_field_y, &particles_arr, x_cords(seq(0, 2)), y_cords(seq(0, 2)), interp);
					T err_dt = (dt_ref > dt) ? dt_ref : dt;
					af::array error = constant(0, particles, tmax / err_dt, aftype);

					errorCalc<T>(particle_ref, particles_arr, &error, ref_length);

					af::array error_t = iota(dim4(n_t - 1)); // * err_dt;

					af::array log_err = log(error(0, seq(1, n_t - 1)));
					if (xx == x_steps - 1 && tt == t_steps - 2 && i == interpol_steps - 1)
					{
						saveToFile<T>("RK4_error1.txt", error(0, span));
					}
					if (xx == x_steps - 1 && tt == t_steps - 3 && i == interpol_steps - 1)
					{
						saveToFile<T>("RK4_error2.txt", error(0, span));
					}
					results(tt + xx * t_steps + i * t_steps * x_steps, 0) = mean(error(0, span));
					results(tt + xx * t_steps + i * t_steps * x_steps, 1) = stdev(error(0, span));
					results(tt + xx * t_steps + i * t_steps * x_steps, 2) = max(error(0, span));

					af::array alpha_factor(2, aftype);// = (log_err(log_err.dims(1) - 1) - log_err(0)) / (((error_t(error_t.dims(0) - 1))));
					alpha_factor = findLine<T>(log(error(0, seq(n_t / 5, n_t - 2))), log(error_t(seq(n_t / 5, n_t - 2))));							  
			
					results(tt + xx * t_steps + i * t_steps * x_steps, 3) = alpha_factor(0);
					results(tt + xx * t_steps + i * t_steps * x_steps, 4) = exp(alpha_factor(1));
					eval(results);
				} //tt
			}
			catch (af::exception &err)
			{
				cout << "Error in testing - loops" << endl;
				cout << err.what() << endl;
			}
		} //xx

	} //interp
	for (int x = 2; x < x_steps; x++)
	{
		for (int tid = 2; tid < t_steps - 1; tid++)
		{

			T res = dx_arr[x];
			T dt = dt_arr[tid];
			int x_len = x_border / res, y_len = y_border / res, n_t = tmax / dt;

			af::array x_cords = iota(dim4(x_len)) * res;
			af::array y_cords = iota(dim4(y_len)) * res;

			af::array particles_double_step(particles, 2, n_t / 2, aftype);
			af::array particles_heuns(particles, 2, n_t, aftype);
			initParticles(delta, &particles_double_step, start, end);
			initParticles(delta, &particles_heuns, start, end);

			af::array vel_field_x(x_len, y_len, n_t, aftype);
			af::array vel_field_y(x_len, y_len, n_t, aftype);
			af::array times = iota(dim4(n_t)) * dt;
			fillFields<T>(&vel_field_x, &vel_field_y, x_cords, y_cords, times, A, epsilon, omega);
			dt = dt * 2.;
			translateParticle<T>(&rk4_integrator, &dt, &tmax, vel_field_x, vel_field_y, &particles_double_step, x_cords, y_cords, AF_INTERP_BICUBIC_SPLINE);
			dt = dt / 2.;
			translateParticle<T>(&heuns_integrator, &dt, &tmax, vel_field_x, vel_field_y, &particles_heuns, x_cords, y_cords, AF_INTERP_BICUBIC_SPLINE);
			T err_dt_double = (dt_ref > dt * 2) ? dt_ref : dt * 2;
			T err_dt_heuns = (dt_ref > dt) ? dt_ref : dt;
			af::array error_double = constant(0, particles, tmax / err_dt_double, aftype), error_heuns = constant(0, particles, tmax / err_dt_heuns, aftype);

			//eval(particles_heuns);
			errorCalc<T>(particle_ref, particles_double_step, &error_double, ref_length);
			eval(error_double);
			errorCalc<T>(particle_ref, particles_heuns, &error_heuns, ref_length);
			eval(error_heuns);
			af::array error_t_heun = iota(dim4(n_t - 1));

			af::array error_t_double = iota(dim4(n_t / 2 - 1));


			af::array log_times_h(log(error_t_heun(seq(1, n_t - 2)))), 
						log_times_d(log(error_t_double(seq(1, n_t / 2 - 2))));
			af::array log_err_h(moddims(log(error_heuns(0, seq(1, n_t - 2))), n_t - 2)), 
						log_err_d(moddims(log(error_double(0, seq(1, n_t / 2 - 2))), n_t / 2 - 2));

			af::array factors_h(2, aftype), factors_d(2, aftype);

			factors_d = findLine<T>(log_err_d(seq(n_t / 5, n_t / 2 - 3)), log_times_d(seq(n_t / 5, n_t / 2 - 3)));

			factors_h = findLine<T>(log_err_h(seq(n_t / 5, n_t - 3)), log_times_h(seq(n_t / 5, n_t - 3)));

			if (x == x_steps - 1 && tid == 1)
			{
				saveToFile<T>("RK4_2xdt_error.txt", error_double(0, span));
				
			}
			if (x == x_steps - 1 && tid == 2)
			{
				saveToFile<T>("heuns_error.txt", error_heuns(0, span));

				
				af::array forward_heuns(particles, 2, 6, aftype);
				initParticles(delta, &forward_heuns, start, end);
				
				heuns_integrator_F_B(&forward_heuns, &vel_field_x, &vel_field_y, &dt, x_cords(seq(0,3)), 
					y_cords(seq(0,3)), 1, AF_INTERP_BICUBIC_SPLINE, false);
				ftle<T>(forward_heuns, dt, tmax, delta, part_x, part_y, false);

				af::array backw_parts(particles, 2, 6, aftype);
				initParticles(delta, &backw_parts, start, end);

				heuns_integrator_F_B(&backw_parts, &vel_field_x, &vel_field_y, &dt, x_cords, 
					y_cords, 1, AF_INTERP_BICUBIC_SPLINE, true);
				ftle<T>(backw_parts, dt, tmax, delta, part_x, part_y, true);
				
				
			}
		
			int offset = t_steps * x_steps * interpol_steps;
			int coord = t_steps * x + tid;
			results(offset + coord, 0) = mean(error_double(0, span));
			results(offset + coord, 1) = stdev(error_double(0, span));
			results(offset + coord, 2) = max(error_double(0, span));
			results(offset + coord, 3) = factors_d(0);
			results(offset + coord, 4) = exp(factors_d(1));
			results(offset + t_steps * x_steps + coord, 0) = mean(error_heuns(0, span));
			results(offset + t_steps * x_steps + coord, 1) = stdev(error_heuns(0, span));
			results(offset + t_steps * x_steps + coord, 2) = max(error_heuns(0, span));
			results(offset + t_steps * x_steps + coord, 3) = factors_h(0);
			results(offset + t_steps * x_steps + coord, 4) = exp(factors_h(1));
			sync(0);
			
		} //tid
	}	 //x

	cout << "Mean   stdev   max   powerlaw   y_0" << endl;
	af_print(results);
	T *res_host = results.host<T>();
	ofstream file_res("error_Npart.txt", std::ios::out | std::ios::binary);
	if (file_res.is_open())
	{
		file_res.write(reinterpret_cast<char *>(res_host), res_size_x * res_size_y * sizeof(T));
	}
	file_res.close();
	free(res_host);
}

template <typename T>
inline void errorCalc(af::array ref, af::array approx, af::array *error, af::array ref_len)
{
	int timesteps = (*error).dims(1);
	int stride;
	int stride2;
	if (approx.dims(2) > ref.dims(2))
	{
		stride = approx.dims(2) / ref.dims(2);
		stride2 = 1;
	}
	else
	{
		stride2 = ref.dims(2) / approx.dims(2);
		stride = 1;
	}
	try
	{
		//af_print(*error);
		for (int i = 1; i < timesteps; i++)
		{
			(*error)(span, i) = (ref(span, 0, i * stride2) - approx(span, 0, stride * i)) * (ref(span, 0, i * stride2) - approx(span, 0, stride * i)) + (ref(span, 1, i * stride2) - approx(span, 1, stride * i)) * (ref(span, 1, i * stride2) - approx(span, 1, stride * i));
			(*error)(span, i) = sqrt((*error)(span, i)) / ref_len(span, i * stride2);
			(*error)(0, i) = mean((*error)(span, i));
		}
		(*error) = (*error) * 100;
	}
	catch (af::exception &er)
	{
		cout << "Failed calculating error" << endl
			 << er.what() << endl;
	}
}

template <typename T>
inline void reference_analytic_rk4(af::array *particles, T dt, T tmax, T A, T omega, T epsilon)
{
	af_dtype type;
	getAfType<T>(&type);

	unsigned long N_part = particles->dims(0);
	dim4 dimens(N_part, 1, 1, 1);
	dim4 dim_retur(1, 1, N_part, 1);

	af::array t = constant(0, dimens, type), tn = constant(0, dimens, type);
	//af::array dt_d = constant(dt, N_part, type);
	int steps = particles->dims(2);
	int index = 0;

	try
	{
		do
		{

			af::array kx1(dimens, type), kx2(dimens, type), kx3(dimens, type), kx4(dimens, type);
			af::array ky1(dimens, type), ky2(dimens, type), ky3(dimens, type), ky4(dimens, type);
			af::array xn(dimens, type), yn(dimens, type), x(dimens, type), y(dimens, type);
			eval((*particles));
			eval((*particles));
			x = (*particles)(span, 0, index);
			y = (*particles)(span, 1, index);
			
			vel_x<T>(A, epsilon, omega, x, y, t, &kx1);
			vel_y<T>(A, epsilon, omega, x, y, t, &ky1);
			
			xn = x + (dt)*kx1 * (T)0.5;
			yn = y + (dt)*ky1 * (T)0.5;
			tn = t + (T)0.5 * dt;

			vel_x<T>(A, epsilon, omega, xn, yn, tn, &kx2);
			vel_y<T>(A, epsilon, omega, xn, yn, tn, &ky2);
			
			xn = x + (dt)*kx2 * (T)0.5;
			yn = y + (dt)*ky2 * (T)0.5;
			
			vel_x<T>(A, epsilon, omega, xn, yn, tn, &kx3);
			vel_y<T>(A, epsilon, omega, xn, yn, tn, &ky3);

			xn = x + (dt)*kx3;

			yn = y + (dt)*ky3;
			tn = t + dt;
			eval(xn, yn, tn);
			vel_x<T>(A, epsilon, omega, xn, yn, tn, &kx4);

			vel_y<T>(A, epsilon, omega, xn, yn, tn, &ky4);
			
			af::array x1(dimens, type), y1(dimens, type);
			eval(kx1, kx2, kx3, kx4);
			x1 = x + (dt) * (kx1 + (T)2.0 * kx2 + (T)2.0 * kx3 + kx4) / ((T)6.0);
			eval(ky1, ky2, ky3, ky4);
			y1 = y + (dt) * (ky1 + (T)2.0 * ky2 + (T)2.0 * ky3 + ky4) / ((T)6.0);
			eval(x1, y1);
			
			(*particles)(span, 0, index + 1) = x1;
			(*particles)(span, 1, index + 1) = y1;
			

			t = tn;
			++index;

		} while (index < steps - 1);
	}
	catch (af::exception &err)
	{
		cout << err.what() << endl;
	}
}
template <typename T>
inline void getLength(af::array *coords, af::array *lengths)
{
	unsigned int len_c = coords->dims(0);
	unsigned int times = coords->dims(2);
	unsigned int len_l = lengths->dims(0);
	if (len_c != len_l)
		throw "Mismatch in dimensions in getLength";
	int i = 0;
	for (i = 1; i < times; i++)
	{
		(*lengths)(span, i) = (*lengths)(span, i - 1) +
							  sqrt(
								  ((*coords)(span, 0, i) - (*coords)(span, 0, i - 1)) *
									  ((*coords)(span, 0, i) - (*coords)(span, 0, i - 1)) +
								  ((*coords)(span, 1, i) - (*coords)(span, 1, i - 1)) *
									  ((*coords)(span, 1, i) - (*coords)(span, 1, i - 1)));
	}
}

template <typename T>
inline void fillFields(af::array *fieldx, af::array *fieldy, af::array x, af::array y, af::array t,
					   T A, T epsilon, T omega)
{
	int n_t = fieldx->dims(2);
	int n_y = fieldx->dims(1);

	T *t_c = t.host<T>();
	T *y_c = y.host<T>();
	af_dtype type = f32;
	dim4 dim(fieldx->dims(0), 1, 1, 1);

	try
	{
		for (int tim = 0; tim < n_t; tim++)
		{
			af::array time = constant(t_c[tim], dim, type);
			for (int yi = 0; yi < n_y; yi++)
			{
				af::array ys = constant(y_c[yi], dim, type);

				(*fieldx)(span, yi, tim) = vel_x<T>(A, epsilon, omega, x, ys, time);
				(*fieldy)(span, yi, tim) = vel_y<T>(A, epsilon, omega, x, ys, time);
			}
		}
	}
	catch (af::exception &err)
	{
		cout << "Error filling out the fields" << endl;
		cout << err.what() << endl;
		cout << "Field x dims:" << fieldx->dims(0) << " " << fieldx->dims(1) << " " << fieldx->dims(2) << " " << endl;
		cout << "Allocating " << sizeof(T) * fieldx->dims(0) * fieldx->dims(1) * fieldx->dims(2) / (1024 * 1024) << " Mbytes" << endl;

	}
	free(t_c);
	free(y_c);
}

template <typename T>
inline af::array vel_x(T A, T epsilon, T omega, af::array x, af::array y, af::array t)
{
	try
	{

		af::array retur = -(T)PI * A * sin((T)PI * funk<T>(epsilon, omega, x, t)) * cos((T)PI * y);
		return retur;
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function vel_x" << endl;
	}
}
template <typename T>
inline af::array vel_y(T A, T epsilon, T omega, af::array x, af::array y, af::array t)
{
	try
	{

		return (T)PI * A * cos((T)PI * funk<T>(epsilon, omega, x, t)) * sin((T)PI * y) * dfunk(epsilon, omega, x, t);
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function vel_y" << endl;
	}
}
template <typename T>
inline void vel_x(T A, T epsilon, T omega, af::array x, af::array y, af::array t, af::array *vel_out)
{
	try
	{

		(*vel_out) = -(T)PI * A * sin((T)PI * funk<T>(epsilon, omega, x, t)) * cos((T)PI * y);
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function vel_x" << endl;
	}
}
template <typename T>
inline void vel_y(T A, T epsilon, T omega, af::array x, af::array y, af::array t, af::array *vel_out)
{
	try
	{

		(*vel_out) = (T)PI * A * cos((T)PI * funk<T>(epsilon, omega, x, t)) * sin((T)PI * y) * dfunk(epsilon, omega, x, t);
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function vel_y" << endl;
	}
}
template <typename T>
inline af::array funk(T epsilon, T omega, af::array x, af::array t)
{
	af::array retur(x.numdims());
	try
	{

		af::array a_t = epsilon * sin(omega * t), b_t = 1 - epsilon * 2 * sin(omega * t);

		retur = a_t * x * x + b_t * x;
		return retur;
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function func" << endl;
		cout << "x_dims in: " << x.dims(0) << x.dims(1) << x.dims(2) << x.dims(3) << endl;
		cout << "dims out: " << retur.dims(0) << retur.dims(1) << retur.dims(2) << retur.dims(3) << endl;
		cout << "t_dims in: " << t.dims(0) << t.dims(1) << t.dims(2) << t.dims(3) << endl;
	}
}
template <typename T>
inline af::array dfunk(T epsilon, T omega, af::array x, af::array t)
{
	try
	{

		af::array a_t, b_t;
		a_t = epsilon * sin(omega * t);
		b_t = 1 - epsilon * 2 * sin(omega * t);
		return 2 * a_t * x + b_t;
	}
	catch (af::exception &er)
	{
		cout << er.what() << endl;
		cout << "crash in function dfunk" << endl;
	}
}

template <typename T>
int saveToFile(const char *filename, af::array out)
{
	T *out_host = out.host<T>();
	ofstream outFile(
		filename,
		std::ios::out | std::ios::binary);
	if (outFile.is_open())
	{
		outFile.write(
			reinterpret_cast<char *>(out_host),
			out.dims(0) * out.dims(1) * out.dims(2) * sizeof(T));
	}
	else
	{
		free(out_host);
		return 1;
	}
	outFile.close();
	free(out_host);
	return 0;
}

template<typename T>
af::array findLine(const af::array lgError, const af::array lgTimes) {
	af_dtype typ;
	getAfType<T>(&typ);

	af::array meany = mean(lgError);
	af::array meanx = mean(lgTimes);

	af::array ret(2, typ);
	af::array beta(1, typ);
	af::array divisor(1, typ);
	beta = 0;
	divisor = 0;
	int N = lgError.dims(0);
	for(int i = 0; i < N; i++) {
		beta += (lgTimes(i) - meanx) * (lgError(i) - meany);
		divisor += (lgTimes(i) - meanx) * (lgTimes(i) - meanx);
	}
	beta = beta / divisor;
	ret(1) = meany - beta * meanx;
	ret(0) = beta;
	eval(ret);
	return ret;
}
