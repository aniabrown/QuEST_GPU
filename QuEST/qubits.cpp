// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

/** @file qubits.c
 * The core of the QuEST Library.
 */

# include <math.h>  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "precision.h"
# include "qubits.h"
# include "qubits_internal.h"

# define DEBUG 0

const char* errorCodes[] = {
    "Success",                                              // 0
    "Invalid target qubit. Note qubits are zero indexed.",  // 1
    "Invalid control qubit. Note qubits are zero indexed.", // 2 
    "Control qubit cannot equal target qubit.",             // 3
    "Invalid number of control qubits",                     // 4
    "Invalid unitary matrix.",                              // 5
    "Invalid rotation arguments.",                          // 6
    "Invalid system size. Cannot print output for systems greater than 5 qubits.", // 7
    "Can't collapse to state with zero probability.", // 8
    "Invalid number of qubits.", // 9
    "Invalid measurement outcome -- must be either 0 or 1." // 10
};

/** Print the current state vector of probability amplitudes for a set of qubits to file.
 * File format:
 * @verbatim
real, imag
realComponent1, imagComponent1
realComponent2, imagComponent2
...
realComponentN, imagComponentN
@endverbatim
 *
 * File naming convention:
 *
 * For each node that the program runs on, a file 'state_rank_[node_rank].csv' is generated. If there is
 * more than one node, ranks after the first do not include the header
 * @verbatim
real, imag
@endverbatim
 * so that files are easier to combine.
 * @param[in,out] multiQubit object representing the set of qubits
 */
void reportState(MultiQubit multiQubit){
	FILE *state;
	char filename[100];
	long long int index;
	sprintf(filename, "state_rank_%d.csv", multiQubit.chunkId);
	state = fopen(filename, "w");
	if (multiQubit.chunkId==0) fprintf(state, "real, imag\n");

	for(index=0; index<multiQubit.numAmps; index++){
		fprintf(state, "%.12f, %.12f\n", multiQubit.stateVec.real[index], multiQubit.stateVec.imag[index]);
	}
	fclose(state);
}

/** Report metainformation about a set of qubits: number of qubits, number of probability amplitudes.
 * @param[in,out] multiQubit object representing the set of qubits
 * @param[in] env object representing the execution environment (local, multinode etc)
 */
void reportMultiQubitParams(MultiQubit multiQubit){
	long long int numAmps = 1L << multiQubit.numQubits;
	long long int numAmpsPerRank = numAmps/multiQubit.numChunks;
	if (multiQubit.chunkId==0){
                printf("QUBITS:\n");
                printf("Number of qubits is %d.\n", multiQubit.numQubits);
                printf("Number of amps is %lld.\n", numAmps);
		printf("Number of amps per rank is %lld.\n", numAmpsPerRank);
    }
}

/** Rotate a single qubit a certain angle about an axis
@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] angle angle by which to rotate in radians
@param[in] unitAxis unit vector pointing along the axis about which to rotate
*/
void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector unitAxis){
    Complex alpha, beta;
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = 0;
    beta.imag = -sin(angle/2.0)*(unitAxis.x + unitAxis.y);
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void rotateX(MultiQubit multiQubit, const int rotQubit, REAL angle){
    Complex alpha, beta;
    Vector unitAxis = {1, 0, 0};
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = 0;
    beta.imag = -sin(angle/2.0)*(unitAxis.x + unitAxis.y);
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void rotateY(MultiQubit multiQubit, const int rotQubit, REAL angle){
    Complex alpha, beta;
    Vector unitAxis = {0, 1, 0};
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = 0;
    beta.imag = -sin(angle/2.0)*(unitAxis.x + unitAxis.y);
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void rotateZ(MultiQubit multiQubit, const int rotQubit, REAL angle){
    Complex alpha, beta;
    Vector unitAxis = {0, 0, 1};
    alpha.real = cos(angle/2.0);
    alpha.imag = -sin(angle/2.0)*unitAxis.z;    
    beta.real = 0;
    beta.imag = -sin(angle/2.0)*(unitAxis.x + unitAxis.y);
    compactUnitary(multiQubit, rotQubit, alpha, beta);
}

void sigmaZ(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, SIGMA_Z);
}

void sGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, S_GATE);
} 

void tGate(MultiQubit multiQubit, const int targetQubit)
{
    phaseGate(multiQubit, targetQubit, T_GATE);
}

int validateMatrixIsUnitary(ComplexMatrix2 u){

    if ( fabs(u.r0c0.real*u.r0c0.real 
        + u.r0c0.imag*u.r0c0.imag
        + u.r1c0.real*u.r1c0.real
        + u.r1c0.imag*u.r1c0.imag - 1) > REAL_EPS ) return 0;
    // check
    if ( fabs(u.r0c1.real*u.r0c1.real 
        + u.r0c1.imag*u.r0c1.imag
        + u.r1c1.real*u.r1c1.real
        + u.r1c1.imag*u.r1c1.imag - 1) > REAL_EPS ) return 0;

    if ( fabs(u.r0c0.real*u.r0c1.real 
        + u.r0c0.imag*u.r0c1.imag
        + u.r1c0.real*u.r1c1.real
        + u.r1c0.imag*u.r1c1.imag) > REAL_EPS ) return 0;

    if ( fabs(u.r0c1.real*u.r0c0.imag
        - u.r0c0.real*u.r0c1.imag
        + u.r1c1.real*u.r1c0.imag
        - u.r1c0.real*u.r1c1.imag) > REAL_EPS ) return 0;

    return 1;
}

int validateAlphaBeta(Complex alpha, Complex beta){
    if ( fabs(alpha.real*alpha.real 
        + alpha.imag*alpha.imag
        + beta.real*beta.real 
        + beta.imag*beta.imag - 1) > REAL_EPS ) return 0;
    else return 1;
}

int validateUnitVector(REAL ux, REAL uy, REAL uz){
    if ( fabs(sqrt(ux*ux + uy*uy + uz*uz) - 1) > REAL_EPS ) return 0;
    else return 1;
}

