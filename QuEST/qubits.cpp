/** @file qubits.c
 * The core of the QuEST Library.
 */

# include <math.h>  //SCB new line
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>
# include "precision.h"
# include "qubits.h"

# define DEBUG 0

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
