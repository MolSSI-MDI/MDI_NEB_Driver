#include <cmath>
#include <vector>
#include <iostream>
namespace neb_utilities {
	void generate_tangent(std::vector<std::vector<double>> coords, std::vector<double> &tangent, std::vector<double> &tangent_pos, std::vector<double> &tangent_neg, std::vector<double> energy, int natoms, int iengine) {
		for (int i = 0; i < natoms*3; i++) {
			tangent_pos[i] = coords[iengine+1][i] - coords[iengine][i];
			tangent_neg[i] = coords[iengine][i] - coords[iengine-1][i];
		}
		double v_i_max = std::max(std::abs(energy[iengine+1] - energy[iengine]), std::abs(energy[iengine-1] - energy[iengine]));
		double v_i_min = std::min(std::abs(energy[iengine+1] - energy[iengine]), std::abs(energy[iengine-1] - energy[iengine]));

		if ( energy[iengine] < energy[iengine + 1] ) {
			if ( energy[iengine] > energy [iengine - 1]) { // V_i+1 > V_i > V_i-1
				for (int i = 0; i < natoms*3; i++) {
					tangent[i] = tangent_pos[i];
				}
			} else { // V_i+1 > V_i < V_i-1
				if ( energy[iengine+1] > energy[iengine-1] ) {
					for (int i = 0; i < natoms*3; i++) {
						tangent[i] = tangent_pos[i]*v_i_max + tangent_neg[i]*v_i_min;
					}
				} else {
					for (int i = 0; i < natoms*3; i++) {
						tangent[i] = tangent_pos[i]*v_i_min + tangent_neg[i]*v_i_max;
					}
				}
			}
		} else { //energy[iengine] > energy[iengine + 1]
			if ( energy[iengine] < energy[iengine - 1] ) { // V_i+1 < V_i < V_i-1
				for (int i = 0; i < natoms*3; i++) {
					tangent[i] = tangent_neg[i];
				}
			} else { // V_i+1 < V_i > V_i-1
				if ( energy[iengine+1] > energy[iengine-1] ) {
					for (int i = 0; i < natoms*3; i++) {
						tangent[i] = tangent_pos[i]*v_i_max + tangent_neg[i]*v_i_min;
					}
				} else {
					for (int i = 0; i < natoms*3; i++) {
						tangent[i] = tangent_pos[i]*v_i_min + tangent_neg[i]*v_i_max;
					}
				}
			}
		}
	}

	void normalize_tangent(std::vector<double> &norm_tan, std::vector<double> tangent) {
			// Normalize the tangent vector
			double tan_mag = 0.0;
			for (int i = 0; i < tangent.size(); i++) {
				tan_mag += pow(tangent[i], 2);
			}
			tan_mag = sqrt(tan_mag);
			for (int i = 0; i < norm_tan.size(); i++) {
				norm_tan[i] = tangent[i] / tan_mag;
			}
	}
			
	void generate_spring_forces(std::vector<double> &spring_forces, double spring_const, std::vector<double> tangent_pos, std::vector<double> tangent_neg, std::vector<double> norm_tan) {
		double mag_tan_pos = 0.0;
		double mag_tan_neg = 0.0;
		int size = spring_forces.size();
		for (int i = 0; i < size; i++) {
			mag_tan_pos += pow(tangent_pos[i], 2);
		}
		mag_tan_pos = sqrt(mag_tan_pos);
		for (int i = 0; i < size; i++) {
			mag_tan_neg += pow(tangent_neg[i], 2);
		}
		mag_tan_neg = sqrt(mag_tan_neg);
		for (int i = 0; i < size; i++) {
			spring_forces[i] = spring_const * (mag_tan_pos - mag_tan_neg) * norm_tan[i];
		}
	}

	void update_forces(std::vector<double> &forces, std::vector<double> norm_tan, std::vector<double> spring_forces, int iengine, bool climbing_phase, int max_engine) {
		int size = spring_forces.size();
		if ((!climbing_phase) || (iengine != max_engine)) {
			// Calculate the True Force
			double dot_prod = 0;
			for (int i = 0; i < size; i++) {
				dot_prod += forces[i] * norm_tan[i];
			}
			for (int i = 0; i < size; i++) {
				forces[i] = spring_forces[i] -	(-forces[i] - dot_prod);
			}
		} else {
			// Calculate the True Force for a Climbing Saddlepoint Node
			double dot_prod = 0;
			for (int i = 0; i < size; i++) {
				dot_prod += forces[i] * norm_tan[i];
			}
			for (int i = 0; i < size; i++) {
				forces[i] = spring_forces[i] -	(-forces[i] - dot_prod);
			}
		}
	}
	
	double local_fnorm_square(std::vector<double> forces, int length) {
		double local_norm2_sqr = 0.0;
		for (int i = 0; i < length; i++) {
			local_norm2_sqr += forces[i]*forces[i];
		}
		return local_norm2_sqr;
	}
	
	double fnorm_square(std::vector<std::vector<double>> forces, int engines) {
		double norm2_sqr = 0.0;
		for ( int i = 0; i < engines; i++) {
			norm2_sqr += local_fnorm_square(forces[i], forces[i].size());
			norm2_sqr = std::max(norm2_sqr, local_fnorm_square(forces[i], forces[i].size()));
		}
		return sqrt(norm2_sqr);
	}
}
