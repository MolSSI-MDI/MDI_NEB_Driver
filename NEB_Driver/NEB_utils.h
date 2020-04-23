#include <cmath>
#include <vector>
#include <iostream>
namespace neb_utilities {
  void generate_tangent(std::vector<std::vector<double>> coords, std::vector<std::vector<double>> cell, std::vector<double> &tangent, std::vector<double> &tangent_pos, std::vector<double> &tangent_neg, std::vector<double> energy, int natoms, int iengine, double &plen, double &nlen, double &dotpath);
	void normalize_tangent(std::vector<double> &norm_tan, std::vector<double> tangent);
	void generate_spring_forces(std::vector<double> &spring_forces, double spring_const, std::vector<double> tangent_pos, std::vector<double> tangent_neg, std::vector<double> norm_tan);
	void update_forces(std::vector<double> &forces, std::vector<double> norm_tan, std::vector<double> spring_forces, double spring_const, int iengine, bool climbing_phase, int max_engine, double plen, double nlen, double dotpath);
	double local_fnorm_square(std::vector<double> forces, int length);
	double fnorm_square(std::vector<std::vector<double>> forces, int engines);
}
