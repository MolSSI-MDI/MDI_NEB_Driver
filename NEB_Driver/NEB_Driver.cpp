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
#include "NEB_utils.h"
using namespace std;

// Method to connect to the engines and store them in the mm_comms array.
void connect_to_engine(vector<MDI_Comm> &mm_comms, int engines) {
    for (int iengine = 0; iengine < engines; iengine++) {
      MDI_Comm comm;
      MDI_Accept_Communicator(&comm);

	// Determine the name of this engine
	char* engine_name = new char[MDI_NAME_LENGTH];
	MDI_Send_Command("<NAME", comm);
	MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
	//  Check to see which engine is connecting to the driver. This driver assumes engine names are passed as 'MM1', 'MM2', 'MM3', ... 'MMn'.    
	if ( ( engine_name[0] == 'M' ) && (engine_name[1] == 'M')) {
		if ( mm_comms[std::atoi(&engine_name[2])-1] != MDI_COMM_NULL ) {
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

// Method to close all the MDI engines within the mm_comms array.
void close_engines(vector<MDI_Comm> &mm_comms, int engines) {
	for (int i = 0; i < engines; i ++) {
		MDI_Send_Command("EXIT", mm_comms[i]);
	}
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
  vector<MDI_Comm> mm_comms(engines, MDI_COMM_NULL);

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
//  double coords[engines][3*natoms] = {0};
  vector<vector<double>> coords(engines, atoms);
//  double energy[engines] = {0};
  vector<double> energy(engines, 0);
  vector<double> cell_engine(9, 0);
  vector<vector<double>> cell(engines, cell_engine);

  //Start a Geometry Optimization for each image.
  for (int iengine = 0; iengine < engines; iengine++) {
	MDI_Send_Command("@INIT_OPTG", mm_comms[iengine]);
  }
  
  // Perform each iteration of the simulation
  bool energy_met = false;
  bool force_met = false;
  int iteration = 0;
  int max_engine = -1;
  bool climbing_phase = false;
  double plen = 0.0; // distance to previous image
  double nlen = 0.0; // distance to next image
  double dotpath = 0.0;
  while ( (!energy_met) || (!force_met) ) {
  //while (iteration <=100) {
	cout << "Timestep: " << iteration << endl;
	double old_energy[engines];
	for (int iengine; iengine < engines; iengine++) {
	  old_energy[iengine] = 0.0;
	}
//	vector<vector<double>> old_forces(engines, atoms);
	//Perform a geometry optimization and give back the forces and coordinates from each node.
	for (int iengine = 0; iengine < engines; iengine++) {    
		//Proceed to the forces node.
		MDI_Send_Command("@FORCES", mm_comms[iengine]);
		//Request and receive the forces from the mm engine. Also store the old forces.
		
//		for (int i = 0; i < 3*natoms; i++) {
//			old_forces[iengine][i] = forces[iengine][i];
//		}
		MDI_Send_Command("<FORCES", mm_comms[iengine]);

		MDI_Recv(&(forces[iengine][0]), 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
		
		//Request and receive the coordinates from the mm engine.
		MDI_Send_Command("<COORDS", mm_comms[iengine]);
		MDI_Recv(&(coords[iengine][0]), 3*natoms, MDI_DOUBLE, mm_comms[iengine]);

		//Request and receive the cell dimensions from the mm engine.
		MDI_Send_Command("<CELL", mm_comms[iengine]);
		MDI_Recv(&(cell[iengine][0]), 9, MDI_DOUBLE, mm_comms[iengine]);
	
		// Request and recieve the energy from the mm engine. Also store the old energy.
		old_energy[iengine] = energy[iengine];
		MDI_Send_Command("<ENERGY", mm_comms[iengine]);
		MDI_Recv(&energy[iengine], 1, MDI_DOUBLE, mm_comms[iengine]);
		cout << "Engine: " << iengine+1 << " Energy: " << std::setprecision(20) << energy[iengine] << endl;
	}

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

	// Perform the NEB operation over all internal replicas.
	for (int iengine = 1; iengine < engines-1; iengine++) {

		// Generate the Tangent for the current image.  
		vector<double> tangent_pos(natoms*3, 0);
		vector<double> tangent_neg(natoms*3, 0);
		vector<double> tangent(natoms*3, 0);
	
	
		//Generate the Tangent for the current image.
		neb_utilities::generate_tangent(coords, cell, tangent, tangent_pos, tangent_neg, energy, natoms, iengine, plen, nlen, dotpath);
			
		// Generate the normalize tangent.
		vector<double> norm_tan(natoms*3, 0);
		neb_utilities::normalize_tangent(norm_tan, tangent);
	
		// Generate the spring forces for the current engine.
		vector<double> spring_forces(natoms*3, 0);
		neb_utilities::generate_spring_forces(spring_forces, spring_const, tangent_pos, tangent_neg, norm_tan);
	
		// Update the forces based on the spring forces.
		neb_utilities::update_forces(forces[iengine], norm_tan, spring_forces, spring_const, iengine, climbing_phase, max_engine, plen, nlen, dotpath);
	}    


//	double max_replica_force = fnorm_square(old_forces, engines);
//	cout << "max_replica_force: " << max_replica_force << endl;
	
	// Send the updated forces back to the engines.
	for (int iengine = 0; iengine < engines; iengine++) { 
		MDI_Send_Command(">FORCES", mm_comms[iengine]);
		MDI_Send(&(forces[iengine][0]), 3*natoms, MDI_DOUBLE, mm_comms[iengine]);
	} 
	
	if (iteration > 0) {
		//Check if we need to perform another iteration based on the energy.
		if (energy_thresh == 0) {
			energy_met = true;
		} else {
			energy_met = true;
			for (int iengine = 0; iengine < engines; iengine++) {
				if ( abs(old_energy[iengine] - energy[iengine]) > energy_thresh) {
					energy_met = false;
				}
			}
		}

	    	//Check if we need to perform another iteration based on the forces..
		if (force_thresh == 0) {
			force_met = true;
		} else {
			force_met = true;
			for (int iengine = 0; iengine < engines; iengine++) {
				for (int i = 0; i < 3*natoms; i++) {
					if (forces[iengine][i] > force_thresh) {
						force_met = false;
					}
				}
			}
		}
    	}
	iteration++;
  	
	//Provide final output.
	  ofstream output_file;
	for ( int iengine = 0; iengine < engines; iengine++) {
		std::string filename = std::to_string(iengine) + "_replica_output.xyz";
		output_file.open(filename, std::ofstream::app);
		output_file << natoms << endl;
		output_file << "Coordinates from iteration " << (iteration-1) << " of an NEB calculation." << endl;
		for (int i = 0; i < natoms; i++) {
			output_file << "H " << coords[iengine][3*i] << " " << coords[iengine][3*i+1] << " " << coords[iengine][3*i+2] << endl;
		}
		output_file.close();
	}
  }


  //Provide final output.
  ofstream output_file;
  std::string filename = "Final_NEB_output.xyz";
  output_file.open(filename, std::ofstream::app);
  for ( int iengine = 0; iengine < engines; iengine++) {
	output_file << natoms << endl;
	output_file << "Replica " << iengine << " final output" << endl;
	for (int i = 0; i < natoms; i++) {
		output_file << "H " << coords[iengine][3*i] << " " << coords[iengine][3*i+1] << " " << coords[iengine][3*i+2] << endl;
	}
  }
  output_file.close();

  // Send the "EXIT" command to each of the engines
  close_engines(mm_comms, engines);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);

  return 0;
}
