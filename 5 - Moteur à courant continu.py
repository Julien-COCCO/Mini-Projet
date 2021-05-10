import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc

# Variables de l'énoncé
R = 5
L = 50*10**(-3)
Ke = 0.2
Kc = 0.1
Fm = 0.01
Jm = 0.05

n = 100000
tfin = 80
h = tfin/n
t = np.linspace(0,tfin,n)

u = np.zeros(n)
for i in range(n):
    if t[i]<=50 and t[i]>=10 :
        u[i] = 5

def moteurCC(Y,t):

    Yprime = np.zeros(2)
    Yprime[0] = (1/L)*(u[t]-R*Y[0]-Ke*Y[1])
    Yprime[1] = (1/Jm)*(Kc*Y[0]-Fm*Y[1])

    return(Yprime)

Y0 = np.array([0,0])

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
    k1 = moteurCC(Yrk[i-1],i-1)
    k2 = Yrk[i-1][:] + (h/2)*k1
    k2 = moteurCC(k2,i-1)
    k3 = Yrk[i-1][:] + (h/2)*k2
    k3 = moteurCC(k3,i-1)
    k4 = Yrk[i-1][:] + h*k3
    k4 = moteurCC(k4,i-1)

    Yrk[i][:] = Yrk[i-1][:] + (h/6)*(k1+2*k2+2*k3+k4)
#

plt.figure(1)
plt.plot(t,Kc*Yrk[:,0],label="Range-Kutta d'ordre 4")
plt.xlabel("Temps t(s)")
plt.ylabel("Couple moteur")
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t,Yrk[:,1],label="Range-Kutta d'ordre 4")
plt.xlabel("Temps t(s)")
plt.ylabel("Vitesse angulaire")
plt.legend()
plt.show()