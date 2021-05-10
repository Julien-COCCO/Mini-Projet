import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc

tfin = 4 # En années
n = round(365*tfin) # Intauration d'une d'une étude jour par jour
h = tfin/n
t = np.linspace(0,tfin,n)

Y0 = np.array([5,3])

## 1 et 2

a1 = 3
b1 = 1
a2 = 2
b2 = 1

y1 =  np.zeros(n)
for i in range(n):
    y1[i] = 5*np.exp(a1*t[i])

y2 = np.zeros(n)
for i in range(n):
    y2[i] = 3*np.exp(-a2*t[i])

plt.figure(1)
#plt.plot(t,y1,label="Proies sans prédateurs")
plt.plot(t,y2,label="Prédateurs sans proies")
plt.xlabel("Temps (années)")
plt.ylabel("Effectif de l'espèce")
plt.legend()
plt.show()

## 3 et +

def modelePP(Y,t):
    a1 = 3
    b1 = 1
    a2 = 2
    b2 = 1

    if Y[0] < 2:
        a1 = 0 # Il faut un effectif minimum de 2 pour que la reproduction est lieu

    if Y[0] < 1:
        Y[0] = 0 # Lorsque l'effectif descend en-dessous de 1, on le fait passer à 0

    if Y[1] < 2:
        b2 = 0 # Il faut un effectif minimum de 2 pour que la reproduction est lieu

    if Y[1] < 1:
        Y[1] = 0 # Lorsque l'effectif descend en-dessous de 1, on le fait passer à 0

    Yprime = np.zeros(2)
    Yprime[0] = a1*Y[0] - b1*Y[0]*Y[1]
    Yprime[1] = -a2*Y[1] + b2*Y[0]*Y[1]

    return(Yprime)


Ye = np.zeros((n,2))
Ye[0][:] = Y0

# Méthode d'Euler
for i in range(1,n):
    Ye[i][:] = Ye[i-1][:] + h*modelePP(Ye[i-1][:],i-1)
#

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
    k1 = modelePP(Yrk[i-1],i-1)
    k2 = Yrk[i-1][:] + (h/2)*k1
    k2 = modelePP(k2,i-1)
    k3 = Yrk[i-1][:] + (h/2)*k2
    k3 = modelePP(k3,i-1)
    k4 = Yrk[i-1][:] + h*k3
    k4 = modelePP(k4,i-1)

    Yrk[i][:] = Yrk[i-1][:] + (h/6)*(k1+2*k2+2*k3+k4)
#


plt.figure("Evolution")
plt.plot(t,Yrk[:,0],label="Evolution des proies (Runge-Kutta d'ordre 4)")
plt.plot(t,Yrk[:,1],label="Evolution des prédateurs (Runge-Kutta d'ordre 4)")
# plt.plot(t,Ye[:,0],label="Evolution des proies (Euler)")
# plt.plot(t,Ye[:,1],label="Evolution des prédateurs (Euler)")
plt.xlabel("Temps (années) (3 1 6 1)")
plt.ylabel("Effectif de l'espèce")
plt.legend()
plt.show()

plt.figure("Portrait de phase")
plt.plot(Yrk[:,0],Yrk[:,1],label="Portrait de phase (Runge-Kutta d'ordre 4)")
# plt.plot(Ye[:,0],Ye[:,1],label="Portrait de phase (Euler)")
plt.xlabel("Lièvres")
plt.ylabel("Lynxs")
plt.legend()
plt.show()









