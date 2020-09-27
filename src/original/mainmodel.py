from qutip import * 
import time
import matplotlib.patches as mpatches

from numpy import *
from qnces import *
import matplotlib.pyplot as plt
from oxiter import * 

att = 0.17 # atteneuation coefficient of fibre
L =  linspace(10,300,30) # length of channel
NN =  linspace(1,20,20) # number of nesting levels
p1 = 0.99 # 1 qubit gate quality
p2 = 0.99 # 2 qubit gate quality
pproj = 0.99 # measurment projector quality
fmin = 0.98 # desired output fidelity
bellmetric = ket2dm(bell_state()) # Bell state to compare fidelities
fidi = empty((0,3))
z = 0

for il in L:
    minn = empty((0,3))
    for xn in NN:
        z = 1
        G = 1
        AN = 0
        AP = 0
        j = 0
        i = 0
        pfib = (10**(-att*(il/(2**xn)/10)))
        print(pfib)
        bellfibre = pol1(bellmetric,qeye(2),2,1,pfib)
        A = (fidelity(bellfibre,bellmetric))**2
        print(A)
        if A > 1:
            A = 1
        B = C = D = (1-A)/3
        while j < xn:
            At = A
            Bt = B
            Ct = C
            Dt = D
            A = connA(At,Bt,Ct,Dt,pproj,p1,p2)
            B = connB(At,Bt,Ct,Dt,pproj,p1,p2)
            C = connC(At,Bt,Ct,Dt,pproj,p1,p2)
            D = connD(At,Bt,Ct,Dt,pproj,p1,p2)
            numseg = (2 ** (xn - 1 - j))
            if A > fmin:
                G = 2
            else:
                while A < fmin:
                    At = A
                    Bt = B
                    Ct = C
                    Dt = D
                    A = coeffA(At,At,Bt,Bt,Ct,Ct,Dt,Dt,pproj,p2)
                    B = coeffB(At,At,Bt,Bt,Ct,Ct,Dt,Dt,pproj,p2)
                    C = coeffC(At,At,Bt,Bt,Ct,Ct,Dt,Dt,pproj,p2)
                    D = coeffD(At,At,Bt,Bt,Ct,Ct,Dt,Dt,pproj,p2)
                    N = pparr(At,At,Bt,Bt,Ct,Ct,Dt,Dt,pproj,p2)
                    G = G*(2/N)
                    i = i + 1
                    if i > 1000:
                        A = fmin
                        G = -1000
                        j = xn
            AN = numseg * G
            AP = AP + AN
            j = j + 1
            i = 0
            G = 1
        if AP > 0:
            minn = append (minn , [[AP , il, xn]] , axis=0)
            print(diagstate(A,B,C,D))
        
    amin = minn[0,0]
    while z < len(minn)-1:
        if minn[z,0] == min(minn[:,0]):
            NT1 = minn[z]
            z = len(minn)
        else: 
            z = z + 1
    fidi = append(fidi , [NT1] , axis=0)
#    fidi2 = fidi
#########################################

# GENERATES FIG (A) ##
plt.xlabel('Distance/km')
plt.ylabel('Bell pair expenditure')
plt.grid()
plt.title("Resource expenditure over distance (A)")
plt.plot(fidi[:,1],fidi[:,0]*8,'g',fidi2[:,1],fidi2[:,0]*7,'r')
plt.legend(['$ES$','$QNC$'], loc='upper left')
plt.savefig("4kEXD1.eps",format="eps")

plt.show()

savetxt("4kEXDes1.csv", fidi[:,0:2], delimiter=",",header = "Expenditure,Distance")
savetxt("4kEXDQNC1.csv", fidi2[:,0:2], delimiter=",",header = "Expenditure,Distance")
    

##########################################


##########################################

# GENERATES FIG (B) ##
#fig, ax1 = plt.subplots()
#ax2 = ax1.twinx()
#ax1.grid()
#plt.title("Closer look at resource jumps between 100-400km (B)")
#ax1.plot(fidi[:,1],fidi[:,2],'go',label='ES Optimal')
#ax1.plot(fidi2[:,1],fidi2[:,2],'ro',label='QNC Optimal')
#ax2.plot(fidi[:,1],fidi[:,0]*8,'g',label='ES')
#ax2.plot(fidi2[:,1],fidi2[:,0]*7,'r',label='QNC')
#ax1.set_xlabel('Distance/km')
#ax1.set_ylabel('Optimal Nesting Level', color='k')
#ax2.set_ylabel('Bell Pair Expenditure', color='k')
#ax1.legend(loc='upper left')
#ax2.legend(loc='lower right')
#plt.savefig("4000OPMNLDIS05.eps",formate="eps")
#plt.show()
#savetxt("400OPMNLDISES05p.csv", fidi, delimiter=",",header = "Expenditure,Distance,OPM Nesting Level")
#savetxt("400OPMNLDISQNC05p.csv", fidi2, delimiter=",",header = "Expenditure,Distance,OPM Nesting Level")


##########################################
