import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc

e = 10 # Tension échelon

# Variables de l'énoncé
C = 10**(-6)
R = 3
L = 0.5
#

n = 10000
tfin = 2
h = tfin/n
Y0 = np.array([0,0])

t = np.linspace(0,tfin,n)

def rlcprim(Y,t):
    return(np.array([ (e-Y[1]-R*Y[0])/L , Y[0]/C ]))

Yrk = np.zeros((n,2))
Yrk[0][:] = Y0

# Création des vecteurs intermédiaires
k1 = np.zeros(2)
k2 = np.zeros(2)
k3 = np.zeros(2)
k4 = np.zeros(2)
#

# Méthode de Runge-Kutta d'ordre 4
for i in range(1,n):
    k1 = rlcprim(Yrk[i-1],i-1)
    k2 = Yrk[i-1][:] + (h/2)*k1
    k2 = rlcprim(k2,i-1)
    k3 = Yrk[i-1][:] + (h/2)*k2
    k3 = rlcprim(k3,i-1)
    k4 = Yrk[i-1][:] + h*k3
    k4 = rlcprim(k4,i-1)

    Yrk[i][:] = Yrk[i-1][:] + (h/6)*(k1+2*k2+2*k3+k4)
#

plt.figure(1)
plt.plot(t,Yrk[:,0],label="Range-Kutta d'ordre 4")
plt.xlabel("Temps t(s)")
plt.ylabel("Intensité")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t,Yrk[:,1],label="Range-Kutta d'ordre 4")
plt.xlabel("Temps t(s)")
plt.ylabel("Tension de sortie")
plt.legend()
plt.show()




