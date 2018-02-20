# ifndef QUBITS_INTERNAL
# define QUBITS_INTERNAL

/** @file
 * Internal functions used to implement the public facing API in qubits.h. Do not call these functions
 * directly. 
 */

void phaseGate(MultiQubit multiQubit, const int targetQubit, enum phaseGateType type);


# endif
