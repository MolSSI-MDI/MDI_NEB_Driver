#include <iostream>
#include <fstream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include <cmath>
extern "C" {
#include "mdi.h"
}
#include <iomanip>
#include <vector>

using namespace std;

// Method to connect to an engine and store it in the mm_comms array.
void connect_to_engine(vector<MDI_Comm> &mm_comms, int engines) {
    for (int iengine = 0; iengine < engines; iengine++) {
	MDI_Comm comm = MDI_Accept_Communicator();

	// Determine the name of this engine
	char* engine_name = new char[MDI_NAME_LENGTH];
	MDI_Send_Command("<NAME", comm);
	MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
	//  Check to see which engine is connecting to the driver. This driver assumes engine names are passed as 'MM1', 'MM2', 'MM3', ... 'MMn'.    
	if ( ( engine_name[0] == 'M' ) && (engine_name[1] == 'M')) {
		if ( mm_comms[std::atoi(&engine_name[2])-1] != MDI_NULL_COMM ) {
		      throw runtime_error("Engine trying to be overritten.");
		}
		// Store the reference to the engine in an ordered list so it can be referenced later.
		mm_comms[std::atoi(&engine_name[2])-1] = comm;
	} 
	else {
	      throw runtime_error("Unrecognized engine name.");
	}
	// De-allocate space for the engine_name.
    	delete[] engine_name;
    }
}

// Method to close all the MDI engines withing the mm_comms array.
void close_engines(vector<MDI_Comm> &mm_comms, int engines) {
	for (int i = 0; i < engines; i ++) {
		MDI_Send_Command("EXIT", mm_comms[i]);
	}
}

void normalize_tangent(double * norm_tan, double * tangent, int size) {
	// Normalize the tangent vector
	double tan_mag;
	for (int i = 0; i < size; i++) {
		tan_mag += tangent[i];
	}
	tan_mag = sqrt(tan_mag);
	for (int i = 0; i < size; i++) {
		norm_tan[i] = tangent[i] / tan_mag;
	}
}

double local_fnorm_square(vector<double> forces, int length) {
	double local_norm2_sqr = 0.0;

	for (int i = 0; i < length; i++) {
		local_norm2_sqr += forces[i]*forces[i];
	}

	return local_norm2_sqr;
}

double fnorm_square(vector<vector<double>> forces, int engines) {
	double norm2_sqr = 0.0;
	for ( int i = 0; i < engines; i++) {
//		norm2_sqr += local_fnorm_square(forces[i], forces[i].size());
		norm2_sqr = std::max(norm2_sqr, local_fnorm_square(forces[i], forces[i].size()));
	}

	return sqrt(norm2_sqr);
}

int main(int argc, char **argv) {

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  // Read through all the command line options
  int iarg = 1;
  int engines = 0;
  bool initialized_mdi = false;
  double spring_const = 0.0;
  double energy_thresh = 0.0;
  double force_thresh = 0.0;
  while ( iarg < argc ) {

    if ( strcmp(argv[iarg],"-mdi") == 0 ) {

      // Ensure that the argument to the -mdi option was provided
      if ( argc-iarg < 2 ) {
	throw runtime_error("The -mdi argument was not provided.");
      }

      // Initialize the MDI Library
      world_comm = MPI_COMM_WORLD;
      int ret = MDI_Init(argv[iarg+1], &world_comm);
      if ( ret != 0 ) {
	throw runtime_error("The MDI library was not initialized correctly.");
      }
      initialized_mdi = true;
      iarg += 2;

    }
    else if ( strcmp(argv[iarg], "-engines") == 0 ) {
      engines = std::stoi(argv[iarg+1]);
      cout << "Engines: " << engines << endl;
      iarg += 2;
    }
    else if ( strcmp(argv[iarg], "-spring") == 0) {
      spring_const = std::stod(argv[iarg+1]);
      cout << "Spring Constant: " << spring_const << endl;
      iarg += 2;
    }
    else if (strcmp(argv[iarg], "-energy_threshold") == 0) {
      energy_thresh = std::stod(argv[iarg+1]);
      cout << "Energy Tolerance: " << energy_thresh << endl;
      iarg+=2;
    }
    else if (strcmp(argv[iarg], "-force_threshold") == 0) {
      force_thresh = std::stod(argv[iarg+1]);
      cout << "Force Tolerance: " << force_thresh << endl;
      iarg+=2;
    }
    else {
      throw runtime_error("Unrecognized option.");
    }

  }
  if ( not initialized_mdi ) {
    throw runtime_error("The -mdi command line option was not provided.");
  }
cout <<"Engines to connect to: " << engines << endl;
 
 
  // Connect to the engines
  vector<MDI_Comm> mm_comms(engines, MDI_NULL_COMM);

  int nengines = engines;
  connect_to_engine(mm_comms, engines);
  // Perform the simulation
  int natoms; // Create a variable to hold the number of atoms

  //Receive the number of atoms from the starting engine; in an NEB, all images are replicas and should contain the same number of atoms..
  MDI_Send_Command("<NATOMS", mm_comms[0]);
  MDI_Recv(&natoms, 1, MDI_INT, mm_comms[0]);

  // Create arrays to store the forces, coordinates, and energies from each engine (replica).
 // double forces[engines][3*natoms];
  vector<double> atoms(3*natoms, 0.0);
  vector<vector<double>> forces(engines, atoms);
  double coords[engines][3*natoms];
  double energy[engines];

  //Start a Geometry Optimization for each node.
  for (int iengine = 0; iengine < engines; iengine++) {
	MDI_Send_Command("OPTG_INIT", mm_comms[iengine]);
  }
  
  // Perform each iteration of the simulation
  bool energy_met = false;
  bool force_met = false;
  int iteration = 0;
  int max_engine = -1;
  bool climbing_phase = false;
//  while ( (!energy_met) || (!force_met) ) {
  while (iteration <=20) {
	cout << "Timestep: " << iteration << endl;
	double old_energy[engines];
//	double old_forces[engines][3*natoms];
	vector<vector<double>> old_forces(engines, atoms);
	//Perform a geometry optimization and give back the forces and coordinates from each node.
	for (int iengine = 0; iengine < engines; iengine++) {    

		//Proceed to the forces node.
		MDI_Send_Command("@FORCES", mm_comms[iengine]);
		//Request and receive the forces from the mm engine. Also store the old forces.
		
		for (int i = 0; i < 3*natoms; i++) {
			old_forces[iengine][i] = forces[iengine][i];
		}
		MDI_Send_Command("<FORCES", mm_comms[iengine]);

//		MDI_Recv(forces[iengine], 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
		MDI_Recv(&(forces[iengine][0]), 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
		
		//Proceed to the coordinates
		MDI_Send_Command("@COORDS", mm_comms[iengine]);

		//Request and receive the coordinates from the mm engine.
		MDI_Send_Command("<COORDS", mm_comms[iengine]);
		MDI_Recv(&coords[iengine], 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
	
		// Request and recieve the energy from the mm engine. Also store the old energy.
		old_energy[iengine] = energy[iengine];
		MDI_Send_Command("<ENERGY", mm_comms[iengine]);
		MDI_Recv(&energy[iengine], 1, MDI_DOUBLE, mm_comms[iengine]);
		cout << "Engine: " << iengine+1 << " Energy: " << std::setprecision(20) << energy[iengine] << endl;
//		cout << "Engine: " << iengine+1 << " Old Energy: " << std::setprecision(20) << old_energy[iengine] << endl;
	}

// neb(int engines, );

// Perform the NEB calculation 
// These exist outside the NEB method
	if (climbing_phase == false) {
		// Find the engine with the highest energy to push to the saddlepoint.
		max_engine = 0;
		for (int iengine = 0; iengine < engines; iengine++) {
			if (energy[iengine] > energy[max_engine]) {
				max_engine = iengine;
			}
		}
	}
/*
	if (iteration > 4 ) {
		climbing_phase = true;
	}
*/

// For the NEB function to perform an iteration, it needs:
// engines, coords, energy, forces

//	neb(engines, coords, energy, forces);	
	for (int iengine = 1; iengine < engines-1; iengine++) {
		// Generate the Tangent for the current image.  
		double tangent_pos[natoms*3];
		double tangent_neg[natoms*3];
		double tangent[natoms*3];

		//Generate the Tangent for the current image.
		//generate_tangent(iengine, natoms, tangent_pos, tangent_neg, tangent, coords, energy);	
	
		for (int i = 0; i < natoms*3; i++) {
			tangent_pos[i] = coords[iengine+1][i] - coords[iengine][i];
			tangent_neg[i] = coords[iengine][i] - coords[iengine-1][i];
		}
		double v_i_max = max(abs(energy[iengine+1] - energy[iengine]), abs(energy[iengine-1] - energy[iengine]));
		double v_i_min = min(abs(energy[iengine+1] - energy[iengine]), abs(energy[iengine-1] - energy[iengine]));

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
		// Generate the normalize tangent.
		double norm_tan[natoms*3];
		normalize_tangent(norm_tan, tangent, 3*natoms);



		double spring_forces[natoms*3];
		double mag_tan_pos;
		double mag_tan_neg;
	
		for (int i = 0; i < natoms*3; i++) {
			mag_tan_pos += tangent_pos[i];
		}
		mag_tan_pos = sqrt(mag_tan_pos);

		for (int i = 0; i < natoms*3; i++) {
			mag_tan_neg += tangent_neg[i];
		}
		mag_tan_neg = sqrt(mag_tan_neg);

		for (int i = 0; i < natoms*3; i++) {
			spring_forces[i] = spring_const * (mag_tan_pos - mag_tan_neg) * norm_tan[i];
		}

		if (!climbing_phase) {	
			// Calculate the True Force
			// i.e. Current force at each atom subtracting the dot product of that force with the tangent vector.
			double temp_force[natoms*3];
			double dot_prod = 0;
			for (int i = 0; i < natoms*3; i++) {
				dot_prod += forces[iengine][i]*norm_tan[i];
			}
			for (int i = 0; i < natoms*3; i++) {
				forces[iengine][i] = spring_forces[i] -	(forces[iengine][i] - dot_prod);
			}
		} else {
			if (iengine == max_engine) {
				// Push the max energy image to the saddle point.
				cout << "Max Engine: " << max_engine+1 << endl;
				double temp_force[natoms*3];
				double dot_prod = 0;
				for (int i = 0; i < natoms*3; i++) {
					dot_prod += forces[iengine][i]*norm_tan[i];
				}
				dot_prod = 2 * dot_prod;
				for (int i = 0; i < natoms*3; i++) {
					forces[iengine][i] = (1 * forces[iengine][i]) + (dot_prod * norm_tan[i]);
				}
			} else {
				// Calculate the True Force
				// i.e. Current force at each atom subtracting the dot product of that force with the tangent vector.
				double temp_force[natoms*3];
				double dot_prod = 0;
				for (int i = 0; i < natoms*3; i++) {
					dot_prod += forces[iengine][i]*norm_tan[i];
				}
				for (int i = 0; i < natoms*3; i++) {
					forces[iengine][i] = spring_forces[i] - (forces[iengine][i] - dot_prod);
				}
			}
		}  
	}    

	double max_replica_force = fnorm_square(old_forces, engines);
//	cout << "max_replica_force: " << max_replica_force << endl;

	// Send the updated forces back to the engines.
	for (int iengine = 0; iengine < engines; iengine++) { 
		MDI_Send_Command(">FORCES", mm_comms[iengine]);
		MDI_Send(&(forces[iengine]), 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
	} 
    
	//Check if we need to perform another iteration.
	if (iteration != 0) {
		if (energy_thresh == 0 ) {
			energy_met == true;
		} else {
			energy_met = true;
			for (int iengine = 0; iengine < engines; iengine++) {
				if ( abs(old_energy[iengine] - energy[iengine]) > energy_thresh) {
					energy_met = false;
//					cout << "Energy not met, perform another iteration." << endl;
				}
			}
		}
	}

    	//Check if we need to perform another iteration.
    	if (iteration != 0) {
		if (force_thresh == 0) {
			force_met = true;
		} else {
			force_met = true;
			for (int iengine = 0; iengine < engines; iengine++) {
				for (int i = 0; i < 3*natoms; i++) {
					if ( abs(old_forces[iengine][i] - forces[iengine][i]) > force_thresh) {
						force_met = false;
//						cout << "Forces not met, perform another iteration." << endl;
					}
				}
			}
		}
    	}
	iteration++;
  	
	//Provide final output.
	  ofstream output_file;
	for ( int iengine = 0; iengine < engines; iengine++) {
		std::string filename = std::to_string(iengine) + "_final_output.xyz";
		output_file.open(filename, std::ofstream::app);
		output_file << natoms << endl;
		output_file << "Final coordinates from an NEB calculation." << endl;
		for (int i = 0; i < natoms; i++) {
			output_file << "H " << coords[iengine][3*i] << " " << coords[iengine][3*i+1] << " " << coords[iengine][3*i+2] << endl;
		}
		output_file.close();
	}
  }



  // Send the "EXIT" command to each of the engines
  close_engines(mm_comms, engines);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);

  return 0;
}
