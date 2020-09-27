#%%
import numpy as np
from qutip.qip.operations import rx, cnot

# Sample representation of 2 step quantum circuit

sample_circuit = [
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
