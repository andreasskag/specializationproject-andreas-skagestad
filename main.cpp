#define NOMINMAX
#include <fstream>
#include <arrayfire.h>
#include <chrono>
//#include <armadillo>
#include "initializing.h"
#include "part_trans.h"
#include "testing.h"

using namespace std;
using namespace af;

// edit here to change datatype - note: if cuda, \
there could be compatability issues when switching to double
typedef float type; 


/**
 * PARAMS: int, int, int, int, float
 * 1. set to 1 if you want to run tests. zero if you want to use a dataset
 * 2. dim_x of the dataset
 * 3. dim_y of the dataset
 * 4. dim_t of the dataset
 * 5. distance between particles(in km).
 * */
int main(int argc, char** argv) {

	//SET benchmark FOR GETTING DATA AND COMPARING WITH ANALYTIC CASE
	bool benchmark = false;


	if(argc > 1) {
		benchmark = stoi(argv[1], 0, 10);

	} else {
		cout << "Function takes args (bool, int, int, int , float)"<< endl;
		cout << "First parameter specifies if you want to run tests, (1 / 0)."
		<<  endl
		<< "The three next specify grid dimensions in the file that is to be loaded (x, y, t)(integers)" << endl 
		<< "The fifth parameter is the spacing between the particles (float)" << endl << 
		"Note: If there are too many particles, there might not be enough memory on the gpu,"<<endl<< 
		"resulting in an error" <<endl;
		return 1;
		//benchmark = 1;
	}
	/////////////
	//SETTING UP ENVIORMENT
	/////////////
	af_dtype aftype;
	
	getAfType<type>(&aftype);
	int device = 0;
	int backends = af::getAvailableBackends();
	//CHECK WHICH BACKENDS ARE AVAILABLE
	bool cuda = backends & AF_BACKEND_CUDA;
	bool cpuu = backends & AF_BACKEND_CPU;
	bool opencl = backends & AF_BACKEND_OPENCL;
	try {
		printf("Available backends: \n CUDA: %d \n OpenCL: %d \n CPU: %d \n", cuda, opencl, cpuu);
		//CHANGE THESE TO FORCE THE USAGE OF A BACKEND
		//NOTE:CRASHES IF BACKEND IS UNAVAILABLE
		//af::setBackend(AF_BACKEND_CPU);
		//af::setBackend(AF_BACKEND_OPENCL);

		//SETTING WHICH DEVICE TO USE. IT IS SET TO ZERO
		af::setDevice(device);
	}
	catch (std::exception &er) {
		cout << er.what() << endl;
		return 1;
	}
	af::info();
	if (!benchmark)
	{
		//USE DATASET
		auto starttime_prog = std::chrono::system_clock::now();
		int x_dim, y_dim, t_dim;
		string fname_u, fname_v;
		type delta;
		af_interp_type interp_t = AF_INTERP_BICUBIC_SPLINE;

		if(argc > 5) {
			x_dim = stoi(argv[2], 0, 10);
			y_dim = stoi(argv[3], 0, 10);
			t_dim = stoi(argv[4], 0, 10);
			delta = stof(argv[5], 0);
		} else {
			x_dim = 991; y_dim = 340; t_dim = 91; //for debug purpose///
			delta = 0.09;
		}

		
		
		cout << "Loading field data..." << endl;
		
	
		/////////////
		//LOADING FILES
		/////////////	

		//LOADING COORDINATE ARRAYS
		type * xcoords_host = new type[x_dim * sizeof(type)];
		type * ycoords_host = new type[y_dim * sizeof(type)];

		ifstream file_x("ocean_data/x_coords_km.txt", std::ios::binary);
		if (file_x.is_open())
		{
			file_x.read(reinterpret_cast<char *>(xcoords_host), (x_dim) * sizeof(type));
		}
		file_x.close();
		ifstream file_y("ocean_data/y_coords_km.txt", std::ios::binary);
		if (file_y.is_open())
		{
			file_y.read(reinterpret_cast<char *>(ycoords_host), (y_dim) * sizeof(type));
		}
		file_y.close();

		af::array af_xcord(x_dim,1,1,1, xcoords_host), af_ycord(y_dim,1,1,1, ycoords_host);

		std::pair<type, type> startxy, endxy;

		startxy.first = xcoords_host[x_dim / 12]; startxy.second = ycoords_host[y_dim / 12];
		endxy.first = xcoords_host[x_dim - x_dim / 12]; endxy.second = ycoords_host[y_dim -y_dim / 12];
		delete[] xcoords_host, ycoords_host;
		//LOADING COORDINATE ARRAYS

		//CALCULATE NUMBER OF PARTICLES
		unsigned int part_x = abs(endxy.first - startxy.first) / delta;
		unsigned int part_y = abs(endxy.second - startxy.second) / delta;
		unsigned long int n_particles = part_x * part_y;


		//LOADING VELOCITY FIELDS ARRAYS
		string fname = "ocean_data/" + to_string(x_dim) + "_" + to_string(y_dim) + "_" + to_string(t_dim) + "_Norkyst_";
		fname_u = fname + "u.txt";
		fname_v = fname + "v.txt";
		type * field_v_host = new type[(x_dim)* (y_dim)* t_dim * sizeof(type)];
		type * field_u_host = new type[(x_dim)* (y_dim)* t_dim * sizeof(type)];
		ifstream file_u(fname_u, std::ios::binary);
		if (file_u.is_open())
		{
			file_u.read(reinterpret_cast<char *>(field_u_host), (x_dim)* (y_dim)* t_dim * sizeof(type));
		}
		file_u.close();
		ifstream file_v(fname_v, std::ios::binary);
		if (file_v.is_open())
		{
			file_v.read(reinterpret_cast<char *>(field_v_host), (x_dim)* (y_dim)* t_dim * sizeof(type));
		}
		file_v.close();

		af::array field_u_device(x_dim, y_dim, t_dim, field_u_host);
		af::array field_v_device(x_dim, y_dim, t_dim, field_v_host);
		/////////////
		//FINISHED LOADING
		/////////////


		delete[] field_u_host, field_v_host;

		auto end_load = std::chrono::system_clock::now();
		auto starttime = std::chrono::system_clock::now();

		//
		type dt = 3600 * 1; 
		type tmax = dt * t_dim; 

		//THIS PARAMETER IS A STRIDE USED IF THE INTEGRATOR SHOULD USE \
		ANOTHER TIMESTEP THAN THE DATASET - HARDCODED FOR DATASETS USING A DT OF 1H. \
		IF A STRIDE IS GOING TO BE USED, IT HAS TO BE SET IN THE dt VARIABLE ABOVE.
		int t_steps = (dt >= 3600) ? (dt) / 3600 : 1;
		int save_t = t_dim / (10);
		//save_t = 1; //SETTING THIS VARIABLE TO 1 TO NOT SAVE ANY INTERMEDIATE STEPS
		dim4 part_dims_short(n_particles, 2, save_t + 1, 1);


		af::array particles_forw(part_dims_short);
		af::array particles_back(part_dims_short);
		


		//DATA IN M/S; CALCS IN KM/S
		field_u_device = field_u_device / 1000;
		field_v_device = field_v_device / 1000;

		end_load = std::chrono::system_clock::now();
		

		starttime = std::chrono::system_clock::now();
		//FUNCTIONS FOR INITIATING THE POSITION OF THE PARTICLES \
		FIRST ARRAY IS FOR THE FORWARDPROPGATING PARTICLES \
		SECOND FOR BACKWARD-PROPAGATING
		initParticles(delta, &particles_forw, startxy, endxy);
		initParticles(delta, &particles_back, startxy, endxy);
		try {
			bool backw = true;
			//CHOICE OF INTEGRATOR IS HARDCODED \
			CHANGE THE COMMENTED SECTION TO USE ANOTHER.
			/*
			heuns_integrator_F_B(&particles_forw, &field_u_device, &field_v_device, &dt, af_xcord, 
			af_ycord, t_steps, interp_t, !backw);
			sync(device);
			heuns_integrator_F_B(&particles_back, &field_u_device, &field_v_device, &dt, af_xcord, 
			af_ycord, t_steps, interp_t, backw);
			*/
			rk4_integrator_F_B(&particles_forw, &field_u_device, &field_v_device, &dt, af_xcord(seq(0,2)), 
				af_ycord(seq(0,2)), t_steps, interp_t, !backw);
			sync(device);	
			rk4_integrator_F_B(&particles_back, &field_u_device, &field_v_device, &dt, af_xcord(seq(0,2)),  
				af_ycord(seq(0, 2)), t_steps, interp_t, backw);	
			
			//THE EVALUATIONS ON THE GPU ARE ASYNC, THEREFORE A SYNC IS NECESSARY FOR TIMING\
			AND FOR BETTER MEMORYMANAGEMENT (SO THAT IT DOESN'T ALLOCATE NEW ARRAYS BEFORE THE \
			PREVIOUS ONES HAVE ALL BEEN FREED)
			sync(device);
		}catch (af::exception &er) {
			cout << er.what() << endl;
			return 1;
		}
		
		std::chrono::duration<double> part_trans = (std::chrono::system_clock::now() - starttime);
		auto start_ftle = std::chrono::system_clock::now();

		//THESE FUNCTIONS CALCULATE THE FTLE-FIELD
		ftle<type>(particles_forw, dt, tmax, delta, part_x, part_y, false);
		ftle<type>(particles_back, dt, tmax, delta, part_x, part_y, true);


		std::chrono::duration<double> load_time = end_load - starttime_prog;
		std::chrono::duration<double> ftle = (std::chrono::system_clock::now() - start_ftle);
		std::chrono::duration<double> wctduration = (std::chrono::system_clock::now() - starttime_prog);
		std::cout << "Loading fields: " << load_time.count() << " seconds [Wall Clock] " << load_time.count() / wctduration.count() * 100 << "% " << std::endl;
		std::cout << "Translating particles: " << part_trans.count() << " seconds [Wall Clock] " << part_trans.count()/ wctduration.count() *100<< "% "<< std::endl;
		std::cout << "FTLE: " << ftle.count() << " seconds [Wall Clock]" << ftle.count() / wctduration.count() * 100 << "% " << std::endl;
		std::cout << "Whole program: " << wctduration.count() << " seconds [Wall Clock]" << std::endl;
	}
	else {
		testing<type>();
	}
	
	
	return 0;
}

