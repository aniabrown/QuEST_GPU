/** @file
An implementation of the API in qubits.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "qubits.h"
# include "qubits_internal.h"

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QUESTEnv env)
{
	createMultiQubitCPU(multiQubit, numQubits, env);
	cudaMalloc(&(multiQubit->deviceStateVec.real), multiQubit->numAmps*sizeof(multiQubit->deviceStateVec.real));
	cudaMalloc(&(multiQubit->deviceStateVec.imag), multiQubit->numAmps*sizeof(multiQubit->deviceStateVec.imag));

        if (!(multiQubit->deviceStateVec.real) || !(multiQubit->deviceStateVec.imag)){
                printf("Could not allocate memory on GPU!\n");
                exit (EXIT_FAILURE);
        }

}

void destroyMultiQubit(MultiQubit multiQubit, QUESTEnv env)
{
	destroyMultiQubitCPU(multiQubit, env);
	cudaFree(multiQubit.deviceStateVec.real);
	cudaFree(multiQubit.deviceStateVec.imag);
}

int GPUExists(void){
	int deviceCount, device;
	int gpuDeviceCount = 0;
	struct cudaDeviceProp properties;
	cudaError_t cudaResultCode = cudaGetDeviceCount(&deviceCount);
	if (cudaResultCode != cudaSuccess) deviceCount = 0;
	/* machines with no GPUs can still report one emulation device */
	for (device = 0; device < deviceCount; ++device) {
		cudaGetDeviceProperties(&properties, device);
		if (properties.major != 9999) { /* 9999 means emulation only */
			++gpuDeviceCount;
		}
	}
	if (gpuDeviceCount) return 1;
	else return 0;
}

void initQUESTEnv(QUESTEnv *env){
        // init MPI environment
	if (!GPUExists()){
		printf("Trying to run GPU code with no GPU available\n");
		exit(EXIT_FAILURE);
	}
	env->rank=0;
	env->numRanks=1;
}

void syncQUESTEnv(QUESTEnv env){
	cudaDeviceSynchronize();
} 

void closeQUESTEnv(QUESTEnv env){
	// MPI finalize goes here in MPI version. Call this function anyway for consistency
}

void reportQUESTEnv(QUESTEnv env){
	printf("EXECUTION ENVIRONMENT:\n");
	printf("Running locally on one node with GPU\n");
	printf("Number of ranks is %d\n", env.numRanks);
# ifdef _OPENMP
	printf("OpenMP enabled\n");
	printf("Number of threads available is %d\n", omp_get_max_threads());
# else
	printf("OpenMP disabled\n");
# endif
}

void copyStateToGPU(MultiQubit multiQubit)
{
	printf("Copying data to GPU\n");
        cudaMemcpy(multiQubit.deviceStateVec.real, multiQubit.stateVec.real, 
			multiQubit.numAmps*sizeof(multiQubit.deviceStateVec.real), cudaMemcpyHostToDevice);
        cudaMemcpy(multiQubit.deviceStateVec.imag, multiQubit.stateVec.imag, 
			multiQubit.numAmps*sizeof(multiQubit.deviceStateVec.imag), cudaMemcpyHostToDevice);
	printf("Finished copying data to GPU\n");
}

void copyStateFromGPU(MultiQubit multiQubit)
{
	cudaDeviceSynchronize();
	printf("Copying data from GPU\n");
        cudaMemcpy(multiQubit.stateVec.real, multiQubit.deviceStateVec.real, 
			multiQubit.numAmps*sizeof(multiQubit.deviceStateVec.real), cudaMemcpyDeviceToHost);
        cudaMemcpy(multiQubit.stateVec.imag, multiQubit.deviceStateVec.imag, 
			multiQubit.numAmps*sizeof(multiQubit.deviceStateVec.imag), cudaMemcpyDeviceToHost);
	printf("Finished copying data from GPU\n");
}


void initStateVec(MultiQubit *multiQubit)
{
	initStateVecCPU(multiQubit);
	copyStateToGPU(*multiQubit);
}

double calcTotalProbability(MultiQubit multiQubit){
        double pTotal=0; 
	long long int index;
	long long int numAmpsPerRank = multiQubit.numAmps;

	copyStateFromGPU(multiQubit);

        for (index=0; index<numAmpsPerRank; index++){ 
                pTotal+=multiQubit.stateVec.real[index]*multiQubit.stateVec.real[index];      
                pTotal+=multiQubit.stateVec.imag[index]*multiQubit.stateVec.imag[index];      
        } 
	return pTotal;
}


__global__ void rotateQubitKernel (MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta){
// ----- sizes
        long long int sizeBlock,                                           // size of blocks
        sizeHalfBlock;                                       // size of blocks halved
        // ----- indices
        long long int thisBlock,                                           // current block
             indexUp,indexLo;                                     // current index and corresponding index in lower half block

        // ----- temp variables
        double   stateRealUp,stateRealLo,                             // storage for previous state values
                 stateImagUp,stateImagLo;                             // (used in updates)
        // ----- temp variables
        long long int thisTask;                                   // task based approach for expose loop with small granularity
        const long long int numTasks=multiQubit.numAmps>>1;
        // (good for shared memory parallelism)


        // ---------------------------------------------------------------- //
        //            tests                                                 //
        // ---------------------------------------------------------------- //
        //assert (rotQubit >= 0 && rotQubit < multiQubit.numQubits);


        // ---------------------------------------------------------------- //
        //            dimensions                                            //
        // ---------------------------------------------------------------- //
        sizeHalfBlock = 1LL << rotQubit;                               // size of blocks halved
        sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks


        // ---------------------------------------------------------------- //
        //            rotate                                                //
        // ---------------------------------------------------------------- //

        //
        // --- task-based shared-memory parallel implementation
        //

        // Can't use multiQubit.stateVec as a private OMP var
	//! fix -- no necessary for GPU version
        double *stateVecReal = multiQubit.deviceStateVec.real;
        double *stateVecImag = multiQubit.deviceStateVec.imag;
        double alphaImag=alpha.imag, alphaReal=alpha.real;
        double betaImag=beta.imag, betaReal=beta.real;

	thisTask = blockIdx.x*blockDim.x + threadIdx.x;
	if (thisTask>=numTasks) return;

	thisBlock   = thisTask / sizeHalfBlock;
	indexUp     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
	indexLo     = indexUp + sizeHalfBlock;

	// store current state vector values in temp variables
	stateRealUp = stateVecReal[indexUp];
	stateImagUp = stateVecImag[indexUp];

	stateRealLo = stateVecReal[indexLo];
	stateImagLo = stateVecImag[indexLo];

	// state[indexUp] = alpha * state[indexUp] - conj(beta)  * state[indexLo]
	stateVecReal[indexUp] = alphaReal*stateRealUp - alphaImag*stateImagUp 
		- betaReal*stateRealLo - betaImag*stateImagLo;
	stateVecImag[indexUp] = alphaReal*stateImagUp + alphaImag*stateRealUp 
		- betaReal*stateImagLo + betaImag*stateRealLo;

	// state[indexLo] = beta  * state[indexUp] + conj(alpha) * state[indexLo]
	stateVecReal[indexLo] = betaReal*stateRealUp - betaImag*stateImagUp 
		+ alphaReal*stateRealLo + alphaImag*stateImagLo;
	stateVecImag[indexLo] = betaReal*stateImagUp + betaImag*stateRealUp 
		+ alphaReal*stateImagLo - alphaImag*stateRealLo;
}

void rotateQubit(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta) 
{
        int threadsPerCUDABlock, CUDABlocks;

        threadsPerCUDABlock = 128;
        CUDABlocks = ceil((double)(multiQubit.numAmps>>1)/threadsPerCUDABlock);
        //printf("cuda blocks: %d\n", CUDABlocks);

        rotateQubitKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, rotQubit, alpha, beta);
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



