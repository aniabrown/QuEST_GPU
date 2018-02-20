# ifndef QUBITS
# define QUBITS

# include "precision.h"

/** @file
 * The QuEST library API and objects. 
*/

/** Represents an array of complex numbers grouped into an array of real components and an array of coressponding complex components.
*/
typedef struct ComplexArray
{
	REAL *real; 
	REAL *imag;
} ComplexArray;

/** Represents one complex number.
*/
typedef struct Complex
{
	REAL real;
	REAL imag;
} Complex;

typedef struct Vector
{
    REAL x, y, z;
} Vector;

/** Represents a system of qubits.
Qubits are zero-based and the the first qubit is the rightmost
*/
typedef struct MultiQubit
{
	//! Probablilty amplitudes for the multi qubit state
	ComplexArray stateVec; 
	//! Temporary storage for a chunk of the state vector received from another process in the MPI version
	ComplexArray pairStateVec;
	//! Storage for probability amplitudes for the multi qubit state on GPU 
	ComplexArray deviceStateVec;
	//! Storage for reduction of probabilities on GPU
	REAL *firstLevelReduction, *secondLevelReduction;
	//! Number of qubits in the state
	int numQubits;
	//! Number of probability amplitudes held in stateVec by this process
	//! In the non-MPI version, this is the total number of amplitudes
	long long int numAmps;
	//! The position of the chunk of the state vector held by this process in the full state vector
	int chunkId;
	//! Number of chunks the state vector is broken up into -- the number of MPI processes used
	int numChunks;
} MultiQubit;

/** Information about the environment the program is running in.
In practice, this holds info about MPI ranks and helps to hide MPI initialization code
*/
typedef struct QuESTEnv
{
	int rank;
	int numRanks;
} QuESTEnv;


// QuEST library functions whose implementation is independent of environment (local, MPI)
void reportState(MultiQubit multiQubit);

void getEnvironmentString(QuESTEnv env, MultiQubit multiQubit, char str[200]);

void reportStateToScreen(MultiQubit multiQubit, QuESTEnv env, int reportRank);

void reportMultiQubitParams(MultiQubit multiQubit);

void controlledPhaseGate(MultiQubit multiQubit, const int idQubit1, const int idQubit2);

void controlledNot(MultiQubit multiQubit, const int controlQubit, const int targetQubit);

// QuEST library functions whose implementation depends on environment (local, MPI)

void createMultiQubit(MultiQubit *multiQubit, int numQubits, QuESTEnv env);

void destroyMultiQubit(MultiQubit multiQubit, QuESTEnv env);

void initStateZero(MultiQubit *multiQubit);

void initStatePlus(MultiQubit *multiQubit);

/** Initialize QuEST environment. If something needs to be done to set up the execution environment, such as 
 * initializing MPI when running in distributed mode, it is handled here
 * @param[in,out] env object representing the execution environment. A single instance is used for each program
 */
void initQuESTEnv(QuESTEnv *env);

/** Close QuEST environment. If something needs to be done to clean up the execution environment, such as 
 * finalizing MPI when running in distributed mode, it is handled here
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void closeQuESTEnv(QuESTEnv env);

/** Guarantees that all code up to the given point has been executed on all nodes. 
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void syncQuESTEnv(QuESTEnv env);

int syncQuESTSuccess(int successCode);

/** Report information about the QuEST environment
 * @param[in] env object representing the execution environment. A single instance is used for each program
 */
void reportQuESTEnv(QuESTEnv env);

/** Calculate the probability of being in any state by taking the norm of the entire state vector. 
 * Should be equal to 1.
 * @param[in] multiQubit object representing a set of qubits
 * @return total probability
 */
REAL calcTotalProbability(MultiQubit multiQubit);

/** Rotate a single qubit in the state vector of probability amplitudes, given the angle rotation arguments.
alphaRe = cos(angle1) * cos(angle2) \n
alphaIm = cos(angle1) * sin(angle2) \n            
betaRe  = sin(angle1) * cos(angle3) \n            
betaIm  = sin(angle1) * sin(angle3) \n           

@remarks Qubits are zero-based and the                     
the first qubit is the rightmost                  
                                                                      
@param[in,out] multiQubit object representing the set of qubits
@param[in] rotQubit qubit to rotate
@param[in] alpha rotation angle
@param[in] beta rotation angle
 */
void compactUnitary(MultiQubit multiQubit, const int rotQubit, Complex alpha, Complex beta);

void sigmaX(MultiQubit multiQubit, const int targetQubit);

void sigmaY(MultiQubit multiQubit, const int targetQubit);

void rotateAroundAxis(MultiQubit multiQubit, const int rotQubit, REAL angle, Vector unitAxis);

/** Measure the probability
of a specified qubit being in the zero state.     

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
REAL findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit);

REAL findProbabilityOfOutcome(MultiQubit multiQubit, const int measureQubit, int outcome);


REAL collapseToOutcome(MultiQubit multiQubit, const int measureQubit, int outcome);

/** Update the state vector to be consistent with measuring measureQubit=0.
Measure in Zero performs an irreversible change to the state vector: it updates the vector according
to the event that a zero have been measured on the qubit indicated by measureQubit (where 
his label starts from 0, of course). It achieves this by setting all inconsistent amplitudes to 0 and 
then renormalising based on the total probability of measuring measureQubit=0. It then returns the 
probability of making this measurement. 

@param[in,out] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
REAL measureInZero(MultiQubit multiQubit, const int measureQubit);




# endif
