#include <cmath>
#include <vector>
#include <iostream>
namespace neb_utilities {
  void generate_tangent(std::vector<std::vector<double>> coords, std::vector<std::vector<double>> cell, std::vector<double> &tangent, std::vector<double> &tangent_pos, std::vector<double> &tangent_neg, std::vector<double> energy, int natoms, int iengine, double &plen, double &nlen, double &dotpath) {
    for (int iatom = 0; iatom < natoms; iatom++) {
      tangent_pos[iatom*3 + 0] = coords[iengine+1][iatom*3 + 0] - coords[iengine][iatom*3 + 0];
      tangent_pos[iatom*3 + 1] = coords[iengine+1][iatom*3 + 1] - coords[iengine][iatom*3 + 1];
      tangent_pos[iatom*3 + 2] = coords[iengine+1][iatom*3 + 2] - coords[iengine][iatom*3 + 2];

      for (int icell = 0; icell < 3; icell++) {
	// get the dot product of tangent_pos with this cell vector
	double dot = cell[iengine][icell*3 + 0] * tangent_pos[iatom*3 + 0] + 
	  cell[iengine][icell*3 + 1] * tangent_pos[iatom*3 + 1] + 
	  cell[iengine][icell*3 + 2] * tangent_pos[iatom*3 + 2];
	double cell_dot = cell[iengine][icell*3 + 0] * cell[iengine][icell*3 + 0] +
	  cell[iengine][icell*3 + 1] * cell[iengine][icell*3 + 1] +
	  cell[iengine][icell*3 + 2] * cell[iengine][icell*3 + 2];
	if ( dot / cell_dot > 0.5 ) {
	  tangent_pos[iatom*3 + 0] -= cell[iengine][icell*3 + 0];
	  tangent_pos[iatom*3 + 1] -= cell[iengine][icell*3 + 1];
	  tangent_pos[iatom*3 + 2] -= cell[iengine][icell*3 + 2];
	}
	else if ( dot / cell_dot < -0.5 ) {
	  tangent_pos[iatom*3 + 0] += cell[iengine][icell*3 + 0];
	  tangent_pos[iatom*3 + 1] += cell[iengine][icell*3 + 1];
	  tangent_pos[iatom*3 + 2] += cell[iengine][icell*3 + 2];
	}
      }

      tangent_neg[iatom*3 + 0] = coords[iengine][iatom*3 + 0] - coords[iengine-1][iatom*3 + 0];
      tangent_neg[iatom*3 + 1] = coords[iengine][iatom*3 + 1] - coords[iengine-1][iatom*3 + 1];
      tangent_neg[iatom*3 + 2] = coords[iengine][iatom*3 + 2] - coords[iengine-1][iatom*3 + 2];

      for (int icell = 0; icell < 3; icell++) {
	// get the dot product of tangent_neg with this cell vector
	double dot = cell[iengine][icell*3 + 0] * tangent_neg[iatom*3 + 0] + 
	  cell[iengine][icell*3 + 1] * tangent_neg[iatom*3 + 1] + 
	  cell[iengine][icell*3 + 2] * tangent_neg[iatom*3 + 2];
	double cell_dot = cell[iengine][icell*3 + 0] * cell[iengine][icell*3 + 0] +
	  cell[iengine][icell*3 + 1] * cell[iengine][icell*3 + 1] +
	  cell[iengine][icell*3 + 2] * cell[iengine][icell*3 + 2];
	if ( dot / cell_dot > 0.5 ) {
	  tangent_neg[iatom*3 + 0] -= cell[iengine][icell*3 + 0];
	  tangent_neg[iatom*3 + 1] -= cell[iengine][icell*3 + 1];
	  tangent_neg[iatom*3 + 2] -= cell[iengine][icell*3 + 2];
	}
	else if ( dot / cell_dot < -0.5 ) {
	  tangent_neg[iatom*3 + 0] += cell[iengine][icell*3 + 0];
	  tangent_neg[iatom*3 + 1] += cell[iengine][icell*3 + 1];
	  tangent_neg[iatom*3 + 2] += cell[iengine][icell*3 + 2];
	}
      }


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

		// calculate plen and nlen
		plen = 0.0;
		nlen = 0.0;
		dotpath = 0.0;
		for (int i = 0; i < natoms*3; i++) {
		  plen += tangent_neg[i] * tangent_neg[i];
		  nlen += tangent_pos[i] * tangent_pos[i];
		  dotpath += tangent_neg[i] * tangent_pos[i];
		}
		plen = sqrt(plen);
		nlen = sqrt(nlen);
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
		  spring_forces[i] = spring_const * (tangent_pos[i] - tangent_neg[i]);
		}
	}

  void update_forces(std::vector<double> &forces, std::vector<double> norm_tan, std::vector<double> spring_forces, double spring_const, int iengine, bool climbing_phase, int max_engine, double plen, double nlen, double dotpath) {
	  int size = spring_forces.size();

	  // calculate dot_prod
	  double dot = 0.0;
	  double dot_prod = 0.0;
	  for (int i = 0; i < size; i++) {
	    dot += forces[i] * norm_tan[i];
	    dot_prod += spring_forces[i] * norm_tan[i];
	  }

	  // calculate the angular contribution
	  double angular_contr = 1.0;
	  double pi = 3.14159265359;
	  angular_contr = 0.5 * (1.0 + std::cos(pi * (dotpath/(plen*nlen)) ));

	  // calculate the prefactor
	  double prefactor = 1.0;
	  if ( climbing_phase && iengine == max_engine ) {
	    // this is the Climbing Saddlepoint Node
	    prefactor = -2.0*dot;
	  }
	  else {
	    prefactor = -dot + spring_const*(nlen-plen);
	  }

	  std::cout << std::endl;
	  std::cout << "IMAGE: " << iengine << std::endl;
	  std::cout << "dotSpringTanget: " << dot_prod << std::endl;
	  std::cout << "tangent: " << norm_tan[411*3+0] << " " << norm_tan[411*3+1] << " " << norm_tan[411*3+2] << std::endl;
	  std::cout << "springF: " << spring_forces[411*3+0] << " " << spring_forces[411*3+1] << " " << spring_forces[411*3+2] << std::endl;
	  std::cout << "AngularContr: " << angular_contr << std::endl;
	  std::cout << "prefactor: " << prefactor/0.00000167580395 << std::endl;
	  std::cout << "kspring: " << spring_const << std::endl;
	  std::cout << "dot: " << dot/0.00000167580395 << std::endl;
	  std::cout << "plen: " << plen << std::endl;
	  std::cout << "nlen: " << nlen << std::endl;
	  for (int i = 0; i < size; i++) {
	    forces[i] += prefactor * norm_tan[i] + angular_contr*( spring_forces[i] - dot_prod * norm_tan[i] );
	  }
	  std::cout << std::endl;
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
