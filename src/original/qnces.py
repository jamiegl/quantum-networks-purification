#   Definitions related to the operation of QNC/ES schemes, and the nature
#   of the error inflicted by them

#   Parameters
#   ---------------
#   A1, A2 - Coefficients of Bell basis state 00 + 11
#   B1, B2 - Coefficients of Bell basis state 01 - 10
#   C1, C2 - Coefficients of Bell basis state 01 + 10
#   D1, D2 - Coefficients of Bell basis state 00 - 11
#   p1 - Probability of a clean one qubit operation in 
#        the presence of depolarizing noise
#   p2 - Probability of a clean two qubit operation in 
#        the presence of depolarizing noise
#   pproj - Probability of projecting onto the correct basis
#           state during a single qubit measurement (1 - pproj
#           is the probability of projecting onto the
#           orthogonal state)

#   Reference
#   ---------------
#   Results taken from "Analysis of Measurement-based Quantum Network Coding
#   over Repeater Networks under Noisy Conditions"

from qutip import *
from oxiter import diagstate

def pol1(qob,gate,N,t,p):
    X = gate_expand_1toN(sigmax(),N,t)
    Y = gate_expand_1toN(sigmay(),N,t)
    Z = gate_expand_1toN(sigmaz(),N,t)
    G = gate_expand_1toN(gate,N,t)
    
    polar1 = (p*(G*qob*G)) + (((1-p)/3) * ((X*qob*X) + (Y*qob*Y) + (Z*qob*Z)))
    return polar1

def pol2(qob,gate,N,c,t,p):
    XY = gate_expand_2toN(tensor([sigmax(),sigmay()]),N,c,t)
    XZ = gate_expand_2toN(tensor([sigmax(),sigmaz()]),N,c,t)
    YX = gate_expand_2toN(tensor([sigmay(),sigmax()]),N,c,t)
    YZ = gate_expand_2toN(tensor([sigmay(),sigmaz()]),N,c,t)
    ZX = gate_expand_2toN(tensor([sigmaz(),sigmax()]),N,c,t)
    ZY = gate_expand_2toN(tensor([sigmaz(),sigmay()]),N,c,t)
    X1 = gate_expand_1toN(sigmax(),N,t)
    Y1 = gate_expand_1toN(sigmay(),N,t)
    Z1 = gate_expand_1toN(sigmaz(),N,t)
    X2 = gate_expand_1toN(sigmax(),N,c)
    Y2 = gate_expand_1toN(sigmay(),N,c)
    Z2 = gate_expand_1toN(sigmaz(),N,c)
    G = gate_expand_2toN(gate,N,c,t)

    polar2 = ((p)*(G*qob*G)) + (((1-p)/12) * ((X1*qob*X1) + (Y1*qob*Y1) \
    + (Z1*qob*Z1) + (X2*qob*X2) + (Y2*qob*Y2) + (Z2*qob*Z2) + (XY*qob*XY) \
    + (YX*qob*YX) + (ZX*qob*ZX) + (XZ*qob*XZ) + (YZ*qob*YZ) + (ZY*qob*ZY)))
    return polar2

def ESc(p1,p2,et,A,B,C,D):
    
    ket0 = Qobj([[1],[0]])
    ket1 = Qobj([[0],[1]])  
    prjk = (et*ket0.proj()) + ((1 - et)*ket1.proj())
    
    wf = diagstate(A,B,C,D)
    inpt = tensor([wf,wf,wf])

    inpt = pol2(inpt,cnot(),6,4,3,p2)
    inpt = gate_expand_1toN(prjk,6,3)*inpt*gate_expand_1toN(prjk,6,3)*2
    inpt = pol2(inpt,cnot(),6,2,1,p2)
    inpt = gate_expand_1toN(prjk,6,1)*inpt*gate_expand_1toN(prjk,6,1)*2
    inpt = pol1(inpt,snot(),6,2,p1)
    inpt = gate_expand_1toN(prjk,6,2)*inpt*gate_expand_1toN(prjk,6,2)*2
    inpt = pol1(inpt,snot(),6,4,p1)
    inpt = gate_expand_1toN(prjk,6,4)*inpt*gate_expand_1toN(prjk,6,4)*2
    
    kq = inpt.ptrace([0,5])
    q = fidelity(kq,bell_state())**2
    return q

def QNCc(p1,p2,et,A,B,C,D):
    
    ket0 = Qobj([[1],[0]])
    ket1 = Qobj([[0],[1]])  
    prjk = (et*ket0.proj()) + ((1 - et)*ket1.proj())
    
    ket0 = Qobj([[1],[0]])
    wf = diagstate(A,B,C,D)
    inpt = tensor([wf,wf,wf,wf,wf,wf,wf])
    
    inpt = pol2(inpt,cnot(),14,13,11,p2)
    inpt = pol2(inpt,cnot(),14,9,7,p2)
    
    inpt = gate_expand_1toN(prjk,14,11)*inpt*gate_expand_1toN(prjk,14,11)*2
    inpt = gate_expand_1toN(prjk,14,7)*inpt*gate_expand_1toN(prjk,14,7)*2
    
    inpt = pol2(inpt,cnot(),14,10,5,p2)
    inpt = pol2(inpt,cnot(),14,6,5,p2)
    
    inpt = gate_expand_1toN(prjk,14,5)*inpt*gate_expand_1toN(prjk,14,5)*2
    inpt = pol2(inpt,cnot(),14,4,3,p2)
    
    inpt = pol2(inpt,cnot(),14,4,1,p2)
    inpt = gate_expand_1toN(prjk,14,3)*inpt*gate_expand_1toN(prjk,14,3)*2
    
    inpt = gate_expand_1toN(prjk,14,1)*inpt*gate_expand_1toN(prjk,14,1)*2
    inpt = pol2(inpt,cnot(),14,2,12,p2)
    
    inpt = pol2(inpt,cnot(),14,0,8,p2)
    
    inpt = pol1(inpt,snot(),14,2,p1)
    inpt = pol1(inpt,snot(),14,0,p1)
    
    inpt = gate_expand_1toN(prjk,14,2)*inpt*gate_expand_1toN(prjk,14,2)*2
    inpt = gate_expand_1toN(prjk,14,0)*inpt*gate_expand_1toN(prjk,14,0)*2
    
    inpt = pol1(inpt,snot(),14,4,p1)
    inpt = gate_expand_1toN(prjk,14,4)*inpt*gate_expand_1toN(prjk,14,4)*2
    
    inpt = pol1(inpt,snot(),14,6,p1)
    inpt = pol1(inpt,snot(),14,10,p1)
    
    inpt = gate_expand_1toN(prjk,14,6)*inpt*gate_expand_1toN(prjk,14,6)*2
    inpt = gate_expand_1toN(prjk,14,10)*inpt*gate_expand_1toN(prjk,14,10)*2
    
    
    a = inpt.ptrace([12,9])
    q = ((fidelity(a,bell_state()))**2)
    return q