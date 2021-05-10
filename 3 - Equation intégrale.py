import numpy as np
import matplotlib.pyplot as plt

N = 5000 # Nombre de points

def K(x,t):
    return((1/np.pi)*1/(1+(x-t)**2))

def f(x):
    return(1)

def u(x):
    return(np.cos(np.pi*x/2))

def Mat(a,b,K,f):

    h = (b-a)/N # Création du pas h
    x = np.linspace(a,b,N) # Création d'un vecteur à N composantes entre a et b

    F = np.zeros(N) # Création du vecteur F
    for i in range(0,N):
        F[i] = f(x[i])

    A = np.zeros((N,N)) # Création du vecteur M

    for i in range(0,N): # Lignes 1 à N
        A[i][0] = K(x[i],x[0]) # Colonne 1
        for j in range(1,N-1): # Colonnes 2 à N-1
            A[i][j] = 2*K(x[i],x[j])
        A[i][N-1] = K(x[i],x[N-1]) # Colonne N

    M = np.eye(N) - (h/2)*A # M = I - (h/2)*A

    return(x,F,M)

x,F,M = Mat(0,10,K,f)

ua = np.zeros(N)
# ue = np.zeros(N)

ua = np.linalg.solve(M,F)

# for i in range(0,N):
#     ue[i] = u(x[i])

# err = np.linalg.norm(ua-ue)/np.sqrt(N)


# print("L'erreur moyenne est de",err)

plt.plot(x,ua,label="Approché")
# plt.plot(x,ue,label="Exact")
plt.legend()
plt.show()








