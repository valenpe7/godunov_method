#ifndef _GODUNOV_
#define _GODUNOV_

#include <iostream>
#include <cmath>  

#include "sim_data.hpp"

class godunov {
    // class containing all the methods necessary for numerical computing system
    // of equations using godunov method
public:
    // constructor getting a pointer to the data structure for data access
	godunov(sim_data* dat);
    
    // launch computation until the end
	void run_computation();
    
    // default destructor
	~godunov() = default;
private:
	sim_data* dat;
    
    // update time step according to CFL condition
	void adapt_delta_t();
    
    // compute fluxes F+ and F- in all domain at one time step
	void compute_fluxes();
    
    // solve the system using godunov method at one time step
	void solve_godunov();
    
    // compute one complete time step
	void compute_step();
};

#endif // _GODUNOV_
