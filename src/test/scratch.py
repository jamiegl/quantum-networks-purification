from qutip.qip.operations import rx, cnot
from qutip import Qobj, tensor, qeye, permute
import numpy as np
import math
from qutip.qip.qubits import qubit_states
from functools import reduce


qubit_ordering = []
operations_in_slice = []
operation_list = None

circuit = [
  {
    "stage": 1,
    "operations": [
      [cnot(), [0,4]],
      [rx(np.pi), [1]]
    ]
  },
  {
    "stage": 2,
    "operations": [
      [cnot(), [1,3]]        
    ]  
  }
]     


qubit_register = qubit_states(N=5, states=[0])


no_of_qubits = math.log(next(x for x in qubit_register.shape if x != 1), 2)


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





# register_mapper = lambda x: 1 if x == register_index else 0

# marked_indices = np.vectorize(register_mapper)(np.arange(no_of_qubits))

# operation_matrix = None