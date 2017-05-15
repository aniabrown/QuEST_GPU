/** @file
An implementation of the API in qubits.h for a local (non-MPI) environment.
*/

# include <stdlib.h>
# include <stdio.h>
# include <omp.h>
# include "qubits.h"
# include "qubits_internal.h"

# define REDUCE_SHARED_SIZE 256

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QUESTEnv env)
{
	createMultiQubitCPU(multiQubit, numQubits, env);
	cudaMalloc(&(multiQubit->deviceStateVec.real), multiQubit->numAmps*sizeof(multiQubit->deviceStateVec.real));
	cudaMalloc(&(multiQubit->deviceStateVec.imag), multiQubit->numAmps*sizeof(multiQubit->deviceStateVec.imag));
	cudaMalloc(&(multiQubit->firstLevelReduction), ceil(multiQubit->numAmps/(double)REDUCE_SHARED_SIZE)*sizeof(double));
	cudaMalloc(&(multiQubit->secondLevelReduction), ceil(multiQubit->numAmps/(double)(REDUCE_SHARED_SIZE*REDUCE_SHARED_SIZE))*
		sizeof(double));

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

        rotateQubitKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, rotQubit, alpha, beta);
}

__device__ __host__ unsigned int log2Int( unsigned int x )
{
        unsigned int ans = 0 ;
        while( x>>=1 ) ans++;
        return ans ;
}

__device__ void reduceBlock(double *arrayIn, double *reducedArray, int length){
        int i, l, r;
        int threadMax, maxDepth;
        threadMax = length/2;
	maxDepth = log2Int(length/2);

        for (i=0; i<maxDepth+1; i++){
                if (threadIdx.x<threadMax){
                        l = threadIdx.x;
                        r = l + threadMax;
                        arrayIn[l] = arrayIn[r] + arrayIn[l];
                }
                threadMax = threadMax >> 1;
                __syncthreads(); // optimise -- use warp shuffle instead
        }

        if (threadIdx.x==0) reducedArray[blockIdx.x] = arrayIn[0];
}

__global__ void copySharedReduceBlock(double*arrayIn, double *reducedArray, int length){
	extern __shared__ double tempReductionArray[];
	int blockOffset = blockIdx.x*blockDim.x;
	tempReductionArray[threadIdx.x*2] = arrayIn[blockOffset + threadIdx.x*2];
	tempReductionArray[threadIdx.x*2+1] = arrayIn[blockOffset + threadIdx.x*2+1];
	__syncthreads();
	reduceBlock(tempReductionArray, reducedArray, length);
}

__global__ void findProbabilityOfZeroKernel(MultiQubit multiQubit,
                const int measureQubit, double *reducedArray)
{
        // ----- sizes
        long long int sizeBlock,                                           // size of blocks
        sizeHalfBlock;                                       // size of blocks halved
        // ----- indices
        long long int thisBlock,                                           // current block
             index;                                               // current index for first half block
        // ----- temp variables
        long long int thisTask;                                   // task based approach for expose loop with small granularity
        long long int numTasks=multiQubit.numAmps>>1;
        // (good for shared memory parallelism)

	extern __shared__ double tempReductionArray[];

        // ---------------------------------------------------------------- //
        //            tests                                                 //
        // ---------------------------------------------------------------- //

        // ---------------------------------------------------------------- //
        //            dimensions                                            //
        // ---------------------------------------------------------------- //
        sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
        // and then the number to skip
        sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

        // ---------------------------------------------------------------- //
        //            find probability                                      //
        // ---------------------------------------------------------------- //

        //
        // --- task-based shared-memory parallel implementation
        //

        double *stateVecReal = multiQubit.deviceStateVec.real;
        double *stateVecImag = multiQubit.deviceStateVec.imag;

	thisTask = blockIdx.x*blockDim.x + threadIdx.x;
	if (thisTask>=numTasks) return;

	thisBlock = thisTask / sizeHalfBlock;
	index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
	double realVal, imagVal;
	realVal = stateVecReal[index];
	imagVal = stateVecImag[index]; 	
	tempReductionArray[threadIdx.x] = realVal*realVal + imagVal*imagVal;
	__syncthreads();

	if (threadIdx.x<blockDim.x/2){
		reduceBlock(tempReductionArray, reducedArray, blockDim.x);
	}
}

double findProbabilityOfZero(MultiQubit multiQubit,
                const int measureQubit)
{
	int numReductionLevels = 1;
	int threadsPerCUDABlock, CUDABlocks, sharedMemSize;
        threadsPerCUDABlock = REDUCE_SHARED_SIZE;
        CUDABlocks = ceil((double)(multiQubit.numAmps>>1)/threadsPerCUDABlock);
	sharedMemSize = threadsPerCUDABlock*sizeof(double);
        findProbabilityOfZeroKernel<<<CUDABlocks, threadsPerCUDABlock, sharedMemSize>>>(multiQubit, measureQubit, 
		multiQubit.firstLevelReduction);

	cudaDeviceSynchronize();	
/*
	threadsPerCUDABlock = CUDABlocks;
	CUDABlocks = 1;
	copySharedReduceBlock<<<threadsPerCUDABlock, CUDABlocks, sharedMemSize>>>(multiQubit.firstLevelReduction, 
		multiQubit.secondLevelReduction, threadsPerCUDABlock); 
	double stateProb=0;
        cudaMemcpy(&stateProb, multiQubit.secondLevelReduction, sizeof(double), cudaMemcpyDeviceToHost);
*/
	double stateProb=0;
        cudaMemcpy(&stateProb, multiQubit.firstLevelReduction, sizeof(double), cudaMemcpyDeviceToHost);
	return stateProb;
}

__global__ void measureInZeroKernel(MultiQubit multiQubit, int measureQubit, double totalProbability)
{
        // ----- sizes
        long long int sizeBlock,                                           // size of blocks
        sizeHalfBlock;                                       // size of blocks halved
        // ----- indices
        long long int thisBlock,                                           // current block
             index;                                               // current index for first half block
        // ----- measured probability
        double   renorm;                                    // probability (returned) value
        // ----- temp variables
        long long int thisTask;                                   // task based approach for expose loop with small granularity
        // (good for shared memory parallelism)
        long long int numTasks=multiQubit.numAmps>>1;

        // ---------------------------------------------------------------- //
        //            tests                                                 //
        // ---------------------------------------------------------------- //
        // ---------------------------------------------------------------- //
        //            dimensions                                            //
        // ---------------------------------------------------------------- //
        sizeHalfBlock = 1LL << (measureQubit);                       // number of state vector elements to sum,
        // and then the number to skip
        sizeBlock     = 2LL * sizeHalfBlock;                           // size of blocks (pairs of measure and skip entries)

        // ---------------------------------------------------------------- //
        //            find probability                                      //
        // ---------------------------------------------------------------- //

        //
        // --- task-based shared-memory parallel implementation
        //
        renorm=1/sqrt(totalProbability);
        double *stateVecReal = multiQubit.deviceStateVec.real;
        double *stateVecImag = multiQubit.deviceStateVec.imag;

	thisTask = blockIdx.x*blockDim.x + threadIdx.x;
	if (thisTask>=numTasks) return;
	thisBlock = thisTask / sizeHalfBlock;
	index     = thisBlock*sizeBlock + thisTask%sizeHalfBlock;
	stateVecReal[index]=stateVecReal[index]*renorm;
	stateVecImag[index]=stateVecImag[index]*renorm;

	stateVecReal[index+sizeHalfBlock]=0;
	stateVecImag[index+sizeHalfBlock]=0;
}

double measureInZero(MultiQubit multiQubit, const int measureQubit)
{        
        double stateProb;
	stateProb = findProbabilityOfZero(multiQubit, measureQubit);

	int threadsPerCUDABlock, CUDABlocks;
        threadsPerCUDABlock = 128;
        CUDABlocks = ceil((double)(multiQubit.numAmps>>1)/threadsPerCUDABlock);
        measureInZeroKernel<<<CUDABlocks, threadsPerCUDABlock>>>(multiQubit, measureQubit, stateProb);
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



