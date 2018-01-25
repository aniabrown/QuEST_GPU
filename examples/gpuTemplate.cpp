/* @file 
 * Basic template for using the QuEST library. In general, leave the initialisation
 * and cleanup sections as they are and edit the rotations, measurement and phase gate
 * sections.
 */

# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <math.h>
# include <unistd.h>
# include <string.h>
# include <omp.h>

# include "QuEST/qubits.h"

//! Max number of angles used to define qubit rotation
# define MaxAngles      10
//! Max number of qubits in the system
# define maxNumQubits   40
//! 1: print end qubit state to file, 0: don't print
# define REPORT_STATE 0


//--------------------------------------------------------------
//---------------------- START OF main()  ----------------------
//--------------------------------------------------------------
int main (int narg, char** varg) {

	//
	// ===== INITIALISATION
	//
	
	// INIT ENVIRONMENT: ALWAYS REQUIRED ONCE AT BEGINNING OF PROGRAM
	// These two lines will automatically set up the environment (multinode,
	// openMP only etc)  
	QuESTEnv env;
	initQuESTEnv(&env);

	// model vars
	int numQubits;
	
	// get number of qubits from command line argument
	if (narg >= 2) {
		numQubits = atoi(varg[1]);
		if (numQubits < 1 || numQubits > maxNumQubits) {
			printf(" *** error: argument %d out of range (1 -- %d)\n", numQubits,maxNumQubits);
			exit (EXIT_FAILURE);
		}
	} else {
		printf(" *** error: too few arguments, number of qubits expected\n");
		exit (EXIT_FAILURE);
	}

	// CREATE QUBIT OBJECT: REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// Before doing any operations on a set of qubits, create the MultiQubit object that will be used to 
	// represent the qubits
	MultiQubit multiQubit; 
	createMultiQubit(&multiQubit, numQubits, env);
	
	// Reporting
	if (env.rank==0) {
		printf("Demo of single qubit rotations.\n");
	}
	reportMultiQubitParams(multiQubit);
	reportQuESTEnv(env);

	// initialise the state to |0000..0>
	initStateZero(&multiQubit);


	//
	// ===== ROTATIONS
	//

	// INITIALISE QUBIT ROTATION
	// Edit these lines to change rotation angle
	double ang1,ang2,ang3;
	Complex alpha, beta;

	// define rotation angles
	double angles[MaxAngles][3] = {
		{ 1.2320,  0.4230, -0.6523},
		{ 2.1213,  0.0000,  3.6520},
		{-3.1213,  5.0230,  0.1230},
		{ 5.2341, -3.1001, -1.2340},
		{-0.1234, -0.9876,  4.1234}
	};

	// rotate
	ang1 = angles[0][0];
	ang2 = angles[0][1];
	ang3 = angles[0][2];

	alpha.real = cos(ang1) * cos(ang2);
	alpha.imag = cos(ang1) * sin(ang2);
	beta.real  = sin(ang1) * cos(ang3);
	beta.imag  = sin(ang1) * sin(ang3);

	int rotQubit;

	// DO QUBIT ROTATION
	if (env.rank==0) printf("\nPerforming qubit rotation\n");
	// Edit these lines to perform rotations as required
	//for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
	for (rotQubit=0; rotQubit<numQubits; rotQubit++) {
		// do rotation of each qubit
		rotateQubit(multiQubit,rotQubit,alpha,beta);
	}
	// END QUBIT ROTATION

	// Verification: check vector size is unchanged
        double totalProbability;
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

        // report state vector to file
	if (REPORT_STATE){
		reportState(multiQubit);
        }

	//
	// ===== perform a measurement
	//
	int measureQubit;
	double qProbability;
	measureQubit=0;
        qProbability = findProbabilityOfZero(multiQubit, measureQubit);
        if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);

	qProbability = measureInZero(multiQubit, 0);
        if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f. Then set q0 to 0\n", 0, qProbability);

        qProbability = findProbabilityOfZero(multiQubit, measureQubit);
        if (env.rank==0) printf("Probability of 0 for qubit %d = %.14f\n", measureQubit, qProbability);


	// filter out qubits
        qProbability = probOfFilterOut111(multiQubit, 0, 1, 2);
        if (env.rank==0) printf("Probability that !(q0=1 & q1=1 & q2=1) = %.14f\n", qProbability);
        totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);
        qProbability = filterOut111(multiQubit, 1, 2, 3);
        if (env.rank==0) printf("Probability that !(q1=1 & q2=1 & q3=1) = %.14f. Also filter out those states\n", qProbability);
        totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);

        qProbability = probOfFilterOut111(multiQubit, 1, 2, 3);
        if (env.rank==0) printf("Probability that !(q1=1 & q2=1 & q3=1) = %.14f\n", qProbability);
        totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);


	//
	// ===== two qubit phase gate
	//

	if (env.rank==0) printf("\nPerforming 2 qubit phase gate\n");	
	controlPhaseGate(multiQubit, 0, 2);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);
	if (env.rank==0) printf("Performing 4 qubit phase gate\n");	
	quadCPhaseGate(multiQubit, 0, 2, 3, 4);
	totalProbability = calcTotalProbability(multiQubit);
        if (env.rank==0) printf("VERIFICATION: total probability=%.14f\n", totalProbability);
	
	//
	// ======== CLEANUP
	//
	
	// free memory

	// REQUIRED ONCE PER MULTIQUBIT OBJECT	
	// When all operations on a set of qubits are completed, destroy the object
	destroyMultiQubit(multiQubit, env);


	// ALWAYS REQUIRED ONCE AT END OF PROGRAM: 
	// These two lines will perform any necessary cleanup of the environment (multinode,
	// openMP only etc)  
	closeQuESTEnv(env);

	return EXIT_SUCCESS;
}

//--------------------------------------------------------------
//----------------------- END OF main()  -----------------------
//--------------------------------------------------------------
