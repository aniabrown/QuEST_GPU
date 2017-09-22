/** @file
An implementation of the API in qubits.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "qubits.h"
# include "precision.h"
# include "qubits_internal.h"

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env)
{
	createMultiQubitCPU(multiQubit, numQubits, env);
}

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env)
{
	destroyMultiQubitCPU(multiQubit, env);
}

void initStateVec(MultiQubit *multiQubit){
	initStateVecCPU(multiQubit);
}

void initQuESTEnv(QuESTEnv *env){
        // init MPI environment
	env->rank=0;
	env->numRanks=1;
}

void syncQuESTEnv(QuESTEnv env){
	// MPI Barrier goes here in MPI version. 
} 

void closeQuESTEnv(QuESTEnv env){
	// MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQuESTEnv(QuESTEnv env){
	printf("EXECUTION ENVIRONMENT:\n");
	printf("Running locally on one node\n");
	printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
	printf("OpenMP enabled\n");
	printf("Number of threads available is %d\n", omp_get_max_threads());
# else
	printf("OpenMP disabled\n");
# endif
}

double calcTotalProbability(MultiQubit multiQubit){

  /* IJB - implemented using Kahan summation for greater accuracy at a slight floating
     point operation overhead. For more details see https://en.wikipedia.org/wiki/Kahan_summation_algorithm */
  /* Don't change the bracketing in this routine! */
  REAL pTotal=0;
  REAL y, t, c;
  long long int index;
  long long int numAmpsPerRank = multiQubit.numAmps;

  c = 0.0;
  for (index=0; index<numAmpsPerRank; index++){
    /* Perform pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index]; by Kahan */
   // pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];

    y = multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index] - c;
    t = pTotal + y;
    c = ( t - pTotal ) - y;
    pTotal = t;

    /* Perform pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index]; by Kahan */
    //pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];


    y = multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index] - c;
    t = pTotal + y;
    c = ( t - pTotal ) - y;
    pTotal = t;


  }
  return pTotal;
}

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta) 
{
	// all values required to update state vector lie in this rank
	rotateQubitLocal(multiQubit, rotQubit, alpha, beta);
}

double findProbabilityOfZero(MultiQubit multiQubit,
                const int measureQubit)
{
	double stateProb=0;
	stateProb = findProbabilityOfZeroLocal(multiQubit, measureQubit);
	return stateProb;
}

double measureInZero(MultiQubit multiQubit, const int measureQubit)
{
        double stateProb;
	stateProb = findProbabilityOfZero(multiQubit, measureQubit);
        measureInZeroLocal(multiQubit, measureQubit, stateProb);
        return stateProb;
}

double filterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        double stateProb=0;
        stateProb = probOfFilterOut111(multiQubit, idQubit1, idQubit2, idQubit3);
        filterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3, stateProb);
        return stateProb;
}

double probOfFilterOut111(MultiQubit multiQubit, const int idQubit1, const int idQubit2, const int idQubit3)
{
        double stateProb=0;
        stateProb = probOfFilterOut111Local(multiQubit, idQubit1, idQubit2, idQubit3);
        return stateProb;
}



