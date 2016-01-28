#include "./inc/godunov.hpp"

godunov::godunov(sim_data* dat)
	:dat(dat) {
}

void godunov::adapt_delta_t() {
	double u_max = 0.0;
	double local;
    // find maximum abs. value of eigenvalues
	for(std::size_t j = 0; j < dat->j_max; ++j) {
		local = std::max(std::abs(dat->U[j][1] / dat->U[j][0] + dat->U[j][0]), std::abs(dat->U[j][1] / dat->U[j][0] - dat->U[j][0]));
		if(local > u_max) {
			u_max = local;
		}
	}
	dat->delta_t = dat->CFL*dat->delta_x/u_max;
	return;
}

void godunov::compute_fluxes() {
	for(std::size_t j = 0; j < dat->j_max; ++j) {
		double m = dat->U[j][1] / pow(dat->U[j][0], 2);
		if(m < -1) {
			dat->F1[j] = {
				0.0,
				0.0
			};
			dat->F2[j] = {
				pow(dat->U[j][0], 2)*m,
				pow(dat->U[j][0], 3)*(pow(m, 2) + 1/3)
			};
		}
		if(m >= -1 && m <= 1) {
			dat->F1[j] = {
				pow(dat->U[j][0]*(m + 1), 2)/4,
				pow(dat->U[j][0]*(m + 1), 3)/6
			};
			dat->F2[j] = {
				-pow(dat->U[j][0]*(m - 1), 2)/4,
				-pow(dat->U[j][0]*(m - 1), 3)/6
			};
		}
		if(m > 1) {
			dat->F1[j] = {
				pow(dat->U[j][0], 2)*m,
				pow(dat->U[j][0], 3)*(pow(m, 2) + 1/3)
			};
			dat->F2[j] = {
				0.0,
				0.0
			};
		}
	}
	return;
}

void godunov::solve_godunov() {
	std::vector<std::array<double, 2>> U_old;
	U_old = dat->U;
	for(std::size_t j = 1; j < dat->j_max - 1; ++j) {
		dat->U[j] = {
			U_old[j][0] - (dat->delta_t / dat->delta_x)*(dat->F1[j][0] - dat->F1[j - 1][0] + dat->F2[j + 1][0] - dat->F2[j][0]),           
			U_old[j][1] - (dat->delta_t / dat->delta_x)*(dat->F1[j][1] - dat->F1[j - 1][1] + dat->F2[j + 1][1] - dat->F2[j][1])
		};
	}
    
    // boundary conditions
	dat->U[0] = {
        U_old[1][0],
        U_old[1][1]
    };     
	dat->U[dat->j_max - 1] = {
        U_old[dat->j_max - 2][0],
        U_old[dat->j_max - 2][1]
    };
    
	return;
}

void godunov::compute_step() {
	this->adapt_delta_t();
	this->compute_fluxes();
	this->solve_godunov();
    
    //advance time
	dat->t += dat->delta_t;
	return;
}

void godunov::run_computation() {
	double t_next_vis = 0.0;
	while(dat->t < dat->t_end) {
        while(dat->t < t_next_vis) {
			this->compute_step();
		}
        dat->print_results();
        t_next_vis += dat->t_vis;
	}
	return;
}