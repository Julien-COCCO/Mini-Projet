import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as sc


n = 10000
tfin = 200
t = np.linspace(0,tfin,n)

def fusee(Y,t):
    D = 4
    a0 = 8*(10**3)
    g = 9.81
    k0 = 0.1
    u = 2*(10**3)

    Yprime = np.zeros(3)

    if Y[1] < 80 :
        Y[1] = 80
        D=0

    if Y[2] < 0 :
        Y[0] = 0
        Y[2] = 0

    Yprime[0] = (D*u)/Y[1] - g - k0*np.exp(-Y[2]/a0)*(Y[0]**2)/Y[1]
    Yprime[1] = -D
    Yprime[2] = Y[0]

    return(Yprime)

Y0 = np.array([0,400,0])

Yode = sc.odeint(fusee,Y0,t)

plt.figure("Vitesse")
plt.plot(t,Yode[:,0])
plt.xlabel("Temps t(s)")
plt.ylabel("Vitesse (m/s)")

plt.figure("Masse")
plt.plot(t,Yode[:,1])
plt.xlabel("Temps t(s)")
plt.ylabel("Masse (kg)")

plt.figure("Position")
plt.plot(t,Yode[:,2])
plt.xlabel("Temps t(s)")
plt.ylabel("Altitude (m)")

plt.show()