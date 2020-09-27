#%%
from qutip import tensor, qeye
import math
from functools import reduce


#%%

def unitary_builder(qubit_register, circuit):
    """
    Builds a unitary transformation representing a quantum circuit
    operating on a register of qubits.
    
    Parameters
    ----------
    qubit_register : Qobj
        Qubit register that will be operated on by circuit.
    circuit : list[dict]
        List of dictionaries describing steps of the circuit

    Returns
    -------
    circuit_unitary : Qobj
        Overall unitary of all operations in the circuit.
    """    
    
    no_of_qubits = math.log(next(x for x in qubit_register.shape if x != 1), 2)
    qubit_ordering = []
    operations_in_slice = []
    operation_list = None
    for slice in circuit:
        for step in slice["operations"]:
            qubit_ordering.extend(step[1])
            operations_in_slice.extend([step[0]])
        identity_operation_count = int(no_of_qubits - len(qubit_ordering))
        operations_in_slice.extend([qeye(2)] * identity_operation_count)
        qubit_ordering.extend([x for x in range(int(no_of_qubits)) if x not in qubit_ordering])
        operation_slice = tensor(operations_in_slice).permute(qubit_ordering)
        if operation_list is None:
            operation_list = [operation_slice]
        else:
            operation_list.extend([operation_slice])
        qubit_ordering = []
        operations_in_slice = [] 
            
    circuit_unitary = reduce((lambda x, y: x * y), operation_list)
    
    return circuit_unitary
        