#include <iostream>
#include <cstdlib>

#include "./inc/godunov.hpp"
#include "./inc/sim_data.hpp"

// set basic properties of computation
sim_data computation_settings(double delta_x, double CFL) {
    
    // set domain from 0 to 1
	std::array<double, 2> x_range = {0.0, 1.0};
    
    // set time range from 0 to 0.5
	std::array<double, 2> t_range = {0.0, 0.5};
    
    // create instance of simulation data with basic properties
	sim_data data(x_range, t_range, delta_x, CFL);
	
    // set initial conditions
    data.initial_conditions(1.0, 0.0, {0.0, 0.5});
    data.initial_conditions(0.1, 0.0, {0.5, 1.0});
    
    // set visualization interval to 0.05
    data.t_vis = 0.05;
    
	return data;
}

int main(int argc, char** argv) {
	if(argc != 3) {
		std::cerr << "Usage: ./application delta_x CFL_number" << std::endl;
		return -1;
	}
	else {
        // set and run computation with
        // delta_x and CFL number as command line arguments
		sim_data data = computation_settings(atof(argv[1]), atof(argv[2]));
		godunov test(&data);
		test.run_computation();
	}
    return 0;
}

