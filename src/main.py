#%%
import func.channels as channels
import circuits.sample as csamp 
from qutip.qip.qubits import qubit_states


register = qubit_states(N=5, states=[0])

circuit = csamp.sample_circuit

unitary = channels.unitary_builder(register, circuit)

result_state = unitary*register