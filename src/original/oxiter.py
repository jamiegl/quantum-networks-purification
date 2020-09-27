#   Definitions related to the operation of the Oxford purification scheme
#   under the influence of depolarizing noise and faulty measurement.

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
#   Results taken from Dür's masters thesis:
#   "Quantum communication over long distances
#   using quantum repeaters"
            
from qutip import ket2dm, bell_state

def diagstate(A,B,C,D):
     
#   Generates a state diagonalized is the Bell basis
   
    dstate = (A * ket2dm(bell_state())) + (B * ket2dm(bell_state('11'))) \
    + (C * ket2dm(bell_state('10'))) + (D * ket2dm(bell_state('01')))
   
    if round(A + B + C + D, 4) != 1:
       raise ValueError("Sum of coefficients must be 1")
    
    return dstate

def pparr(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2):
      
#   Calculates success probability of Oxford purification 
   
    N = (4 * (A1 + (A2 * B1) + (B2 * B1) + (B1 * C2) + (A2 * C1) + (B2 * C1) \
    + (C2 * C1) + (B1 * D2) + (C1 * D2) + (A2 * D1) + (B2 * D1) + (C2 * D1) \
    + (D2 * D1) + ((p2**2) * ((A2 * A1) + (A1 * B2) + (A2 * B1) + (B2 * B1) \
    - (A1 * C2) - (B1 * C2) - (A2 * C1) - (B2 * C1) + (C2 * C1) \
    - (A1 * D2) - (B1 * D2) + (C1 * D2) - (A2 * D1) - (B2 * D1) \
    + (C2 * D1) + (D2 * D1) - (4 * A2 * A1 * pproj) \
    - (4 * A1 * B2 * pproj) - (4 * B2 * B1 * pproj) \
    + (4 * A1 * C2 * pproj) + (4 * B1 * C2 * pproj) \
    + (4 * A2 * C1 * pproj) + (4 * B2 * C1 * pproj) \
    - (4 * C2 * C1 * pproj) + (4 * A1 * D2 * pproj) \
    + (4 * B1 * D2 * pproj) - (4 * C1 * D2 * pproj) \
    + (4 * A2 * D1 * pproj) + (4 * B2 * D1 * pproj) \
    - (4 * C2 * D1 * pproj) - (4 * D2 * D1 * pproj) \
    - (4 * A2 * B1 * pproj) \
    + (4 * A2 * A1 * (pproj**2)) + (4 * A1 * B2 * (pproj**2)) \
    + (4 * A2 * B1 * (pproj**2)) + (4 * B2 * B1 * (pproj**2)) \
    - (4 * A1 * C2 * (pproj**2)) - (4 * A2 * C1 * (pproj**2)) \
    - (4 * B2 * C1 * (pproj**2)) + (4 * C2 * C1 * (pproj**2)) \
    - (4 * A1 * D2 * (pproj**2)) - (4 * B1 * D2 * (pproj**2)) \
    + (4 * C1 * D2 * (pproj**2)) - (4 * A2 * D1 * (pproj**2)) \
    - (4 * B2 * D1 * (pproj**2)) + (4 * C2 * D1 * (pproj**2)) \
    + (4 * D2 * D1 * (pproj**2)) - (4 * B1 * C2 * (pproj**2))))))/8
                 
    return N
          
 
def coeffA(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2):
        
#   Calculates weighting of 00+11 state in ensemble on iteration 
   
    Aprime = (A1 + (A2 * B1) + (B2 * B1) + (B1 * C2) + (B2 * C1) + (C2 * C1) \
    + (B1 * D2) + (C1 * D2) + (A2 * D1) + (B2 * D1) + (C2 * D1) \
    + (D2 * D1) + (A2 * C1) + ((p2**2)*((7 * A2 * A1) - (A1 * B2) - (A2 * B1) \
    + (7 * B2 * B1) - (A1 * C2) - (B1 * C2) - (A2 * C1) - (B2 * C1) \
    - (C2 * C1) - (A1 * D2) - (B1 * D2) - (C1 * D2) - (A2 * D1) - (B2 * D1) \
    - (C2 * D1) - (D2 * D1) - (16 * A2 * A1 * pproj) - (16 * B2 * B1 * pproj) \
    + (16 * A2 * C1 * pproj) + (16 * B2 * D1 * pproj) \
    + (16 * A2 * A1 * (pproj**2)) + (16 * B2 * B1 * (pproj**2))
    - (16 * A2 * C1 *(pproj**2)) - (16 * B2 * D1 * (pproj**2))))) \
    / (8 * pparr(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2))
    
    return Aprime
    

def coeffB(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2):

#   Calculates weighting of 01-10 state in ensemble on iteration 
 
    Bprime = (A1 + (A2 * B1) + (B2 * B1) + (A2 * C1) + (B2 * C1) + (C2 * C1) \
    + (B2 * D1) + (C2  * D1) + (D2 * D1) + (B1 * C2) + (B1 * D2) + (C1 * D2) \
    + (A2 * D1) + ((p2**2) * (-(A2 * A1) - (A1 * B2) - (A2 * B1) - (B2 * B1) \
    - (A1 * C2) - (B1 * C2) - (A2 * C1) - (B2 * C1) - (C2 * C1) - (A1 * D2) \
    - (B1 * D2) + (7 * C1 * D2) - (A2 * D1) - (B2 * D1) + (7 * C2 * D1) \
    - (D2 * D1) + (16 * B1 * C2 * pproj) + (16 * A1 * D2 * pproj) \
    - (16 * C1 * D2 * pproj) - (16 * C2 * D1 * pproj) \
    - (16 * B1 * C2 * (pproj**2)) - (16 * A1 * D2 * (pproj**2)) \
    + (16 * C1 * D2 * (pproj**2)) + (16 * C2 * D1 * (pproj**2))))) \
    / (8 * pparr(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2))
    
    return Bprime
    

def coeffC(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2):

#   Calculates weighting of 01+10 state in ensemble on iteration 

    Cprime = (A1 + (A2 * B1) + (B2 * B1) + (B1 * C2) + (A2 * C1) + (B2 * C1) \
    + (C2 * C1) + (B1 * D2) + (C1 * D2) + (A2 * D1) + (B2 * D1) + (C2 * D1) \
    + (D2 * D1) + ((p2**2) * (- (A2 * A1) - (A1 * B2) - (A2 * B1) - (B2 * B1) \
    - (A1 * C2) - (B1 * C2) - (A2 * C1) - (B2 * C1) + (7 * C2 * C1) \
    - (A1 * D2) - (B1 * D2) - (C1 * D2) - (A2 * D1) - (B2 * D1) - (C2 * D1) \
    + (7 * D2 * D1) + (16 * A1 * C2 * pproj) - (16 * C2 * C1 * pproj) \
    + (16 * B1 * D2 * pproj) - (16 * D2 * D1 * pproj) \
    - (16 * A1 * C2 * (pproj**2)) + (16 * C2 * C1 * (pproj**2))
    - (16 * B1 * D2 * (pproj**2)) + (16 * D2 * D1 * (pproj**2))))) \
    / (8 * pparr(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2))
    
    return Cprime

def coeffD(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2):

#   Calculates weighting of 00-11 state in ensemble on iteration 

    Dprime = (A1 + (A2 * B1) + (B2 * B1) + (B1 * C2) + (A2 * C1) + (B2 * C1) \
    + (C2 * C1) + (B1 * D2) + (C1 * D2) + (A2 * D1) + (B2 * D1) + (C2 * D1) \
    + (D2 * D1) + ((p2**2) * (-(A2 * A1) + (7 * A1 * B2) + (7 * A2 * B1) \
    - (B2 * B1) - (A1 * C2) - (B1 * C2) - (A2 *  C1) - (B2 * C1) - (C2 * C1) \
    - (A1 * D2) - (B1 * D2) - (C1 * D2) - (A2 * D1) - (B2 * D1) - (C2 * D1) \
    - (D2 * D1) - (16 * A1 * B2 * pproj) - (16 * A2 * B1 * pproj) \
    + (16 * B2 * C1 * pproj) + (16 * A2 * D1 * pproj) \
    + (16 * A1 * B2 * (pproj**2)) + (16 * A2 * B1 * (pproj**2)) \
    - (16 * B2 * C1 * (pproj**2)) - (16 * A2 * D1 * (pproj**2))))) \
    / (8 * pparr(A1,A2,B1,B2,C1,C2,D1,D2,pproj,p2)) 
    
    return Dprime
    
def connA(A,B,C,D,pproj,p1,p2):

#   Calculates weighting of 00+11 state in ensemble on connection 

    Aprimec = p1 * (p2 * (((pproj**2)*(A**2 + B**2 + C**2 + D**2)) + (((1 - pproj)**2) \
    * ((2 * A * B) + (2 * C * D))) + ((2 * pproj) * (1 - pproj) * ((A + B) \
    * (C + D)))) + ((1 - p2)/4)) + ((1 - p1)/4)
    
    return Aprimec

def connB(A,B,C,D,pproj,p1,p2):

#   Calculates weighting of 01-10 state in ensemble on connection 

    Bprimec = p1 * (p2 * (((pproj**2)*(2 * A * B) + (2 * C * D)) + (((1 - pproj)**2) \
    * ((A**2 + B**2) + (C**2 + D**2))) + ((2 * pproj) * (1 - pproj) * (A + B) \
    * (C + D))) + ((1 - p2)/4)) + ((1 - p1)/4) \
    
    return Bprimec
    

def connC(A,B,C,D,pproj,p1,p2):

#   Calculates weighting of 01+10 state in ensemble on connection 
   
    Cprimec = p1 * (p2 * (((pproj**2)*(2 * A * C) + (2 * B * D)) + (((1 - pproj)**2) \
    * ((2 * A * D) + (2 * B * C))) + (pproj * (1 - pproj) * (A**2 + B**2 \
    + C**2 + D**2 + (2 * A * B) + (2 * C * D)))) + ((1 - p2)/4)) + ((1 - p1)/4) \
    
    return Cprimec
    
def connD(A,B,C,D,pproj,p1,p2):

#   Calculates weighting of 00-11 state in ensemble on connection 

    Dprimec = p1 * (p2 * (((pproj**2)*(2 * A * D) + (2 * B * C)) + (((1 - pproj)**2) \
    * ((2 * A * C) + (2 * B * D))) + (pproj * (1 - pproj) * (A**2 + B**2 \
    + C**2 + D**2 + (2 * A * B) + (2 * C * D)))) + ((1 -p2)/4)) + ((1 -p1)/4) \
    
    return Dprimec
    

    
    
    
    
    