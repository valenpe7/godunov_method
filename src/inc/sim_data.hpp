#ifndef _DATA_
#define _DATA_

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>

class sim_data {
    // class containing all the data necessary for computing system of equations
    // using flux splitting and godunov method
	friend class godunov;
public:
    // standard constructor, filling the data with initial or defalt values 
	sim_data(std::array<double, 2> x_range, std::array<double, 2> t_range, double delta_x, double CFL);
    
    // set initial conditions at chosen interval
	sim_data& initial_conditions(double rho_init, double u_init, std::array<double, 2> interval);
	
    // write data into file for visualization purposes
    void print_results();
    
    // default destructor
	~sim_data() = default;
public:
    double t_vis;                                          // visualization time step
private:
	std::vector<std::array<double, 2>> U;                  // vector of rho, rho*u
	std::vector<std::array<double, 2>> F1;                 // Flux F+
	std::vector<std::array<double, 2>> F2;                 // Flux F-
    std::vector<double> x;                                 // discretized space coordinate
	unsigned int j_max;                                    // number of cells
    double t;                                              // current time value
    double t_end;                                          // time when the computation ends
    double delta_x;                                        // grid size
	double delta_t;                                        // time step
	double CFL;                                            // CFL number
};

#endif // _DATA_
