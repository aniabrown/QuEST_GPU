// Distributed under MIT licence. See https://github.com/aniabrown/QuEST_GPU/blob/master/LICENCE.txt for details

# ifndef QUBITS_INTERNAL
# define QUBITS_INTERNAL

/** @file
 * Internal functions used to implement the public facing API in qubits.h. Do not call these functions
 * directly. 
 */

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);

/** Measure the probability
of a specified qubit being in the zero state.     

@param[in] multiQubit object representing the set of qubits
@param[in] measureQubit qubit to measure
@return probability of qubit measureQubit being zero
*/
REAL findProbabilityOfZero(MultiQubit multiQubit, const int measureQubit);


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
