#include "./inc/sim_data.hpp"

sim_data::sim_data(std::array<double, 2> x_range, std::array<double, 2> t_range, double delta_x, double CFL)
	:delta_x(delta_x), CFL(CFL) {
	this->t = t_range[0];
    this->t_end = t_range[1];
    this->t_vis = t_end / 10;
	this->j_max = ceil((x_range[1] - x_range[0]) / delta_x);
	this->U.resize(j_max);
	this->F1.resize(j_max);
	this->F2.resize(j_max);
    this->x.resize(this->j_max);
    for(std::size_t j = 0; j < this->j_max; ++j) {
        x[j] = x_range[0] + j*delta_x;
    }
}

sim_data& sim_data::initial_conditions(double rho_init, double u_init, std::array<double, 2> interval) {
	for(std::size_t j = 0; j < this->j_max; ++j) {
		if(this->x[j] >= interval[0] && this->x[j] <= interval[1]) {
			this->U[j][0] = rho_init;
			this->U[j][1] = rho_init * u_init;
		}
	}
	return (*this);
}

void sim_data::print_results() {
	std::string output;
	std::ostringstream stream;
	stream << "solution_t=" << std::setprecision(3) << std::fixed << this->t << ".dat";
	output = stream.str();
	std::ofstream file(output, std::ios::out);
	if(file) {
		for(std::size_t j = 0; j < this->j_max; ++j) {
			file
				<< std::setprecision(5)
				<< std::scientific
				<< this->x[j] << " "
				<< this->U[j][0] << " "
				<< this->U[j][1]
				<< std::endl;
		}
		file.close();
	} else {
		std::cerr << "error: unable write to file " << output << std::endl;
		return;
	}
	std::cout << "data written into file " << output << std::endl;
	return;
}