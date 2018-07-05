#include <stdio.h>
#include "QuEST/qubits.h"

int main (int narg, char *varg[]) {

    /*
     * PREPARE QuEST environment
     * (Required only once per program)
     */

    QuESTEnv env;
    initQuESTEnv(&env);

    printf("-------------------------------------------------------\n");
    printf("Running QuEST tutorial:\n\t Basic circuit involving a system of 3 qubits.\n");
    printf("-------------------------------------------------------\n");

    /*
     * PREPARE QUBIT SYSTEM
     */

    MultiQubit qubits; 
    createMultiQubit(&qubits, 3, env);
    initStateZero(&qubits);


    /*
     * REPORT SYSTEM AND ENVIRONMENT
     */
    printf("\nThis is our environment:\n");
    reportMultiQubitParams(qubits);
    reportQuESTEnv(env);

    /*
     * APPLY CIRCUIT
     */

    controlledNot(qubits, 0, 1);

    Complex a, b;
    a.real = .5; a.imag =  .5;
    b.real = .5; b.imag = -.5;
    compactUnitary(qubits, 1, a, b);

    Vector v;
    v.x = 1; v.y = 0; v.z = 0;
    rotateAroundAxis(qubits, 2, 3.14/2, v);

    printf("\nCircuit output:\n");

    REAL prob = findProbabilityOfOutcome(qubits, 2, 1);
    printf("Probability of qubit 2 being in state 1: %lf\n", prob);


    /*
     * FREE MEMORY
     */

    destroyMultiQubit(qubits, env); 


    /*
     * CLOSE QUEST ENVIRONMET
     * (Required once at end of program)
     */
    closeQuESTEnv(env);
    return 0;
}
