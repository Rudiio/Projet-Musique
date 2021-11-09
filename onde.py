#Résolution de l'équation de d'Alembert par la méthode des différences finies

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import*

#Création de la matrice tridiagonale  
def M(n):
  """int -> array(n,n)
 Renvoie matrice diagonale de dimension n*n contenant les coefficients du systeme"""

  M = zeros((n+1,n+1),float)

  for i in range(0,n):
    if(i==0):                                   #pour la première ligne de la matrice
      M[i,i], M[i,i+1]=(2-2*alpha**2),alpha**2
    elif(i==n-1):
      M[i,i-1], M[i,i]=alpha**2,(2-2*alpha**2)
    else:                                       #Pour la dernière ligne de la matrice
      M[i,i-1], M[i,i], M[i,i+1] = alpha**2,(2-2*alpha**2),alpha**2

  return M


#Initialise la matrice colonne avec les valeurs à l'état initial t=0
#Ici,on considère que U(x,0)=sin(x)
def init_U0(n):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n+1,1),float)

  for i in range(0,n+1):
    U0[i,0]=sin((MODE*pi/L)*(i*dx));
  U0[n,0]=0
  return U0


#Initialise la matrice colonne avec les valeurs à l'état t=1
#le calcul s'effectue avec la formule trouvée en discrétisant
def init_U1(n):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=1 """
  #Initialisation des tableaux
  U0=init_U0(n+1)
  U1=zeros((n+1,1),float)

  for i in range (1,n):
    U1[i]= alpha**2*(U0[i-1] - 2*U0[i]+U0[i+1])+2*U0[i]
  U1[0]=U1[n]=0.0 #Conditions initiales à modifier 
  
  return U1

#caractéristique physique 
c=340 #m.s^-1   célérité de l'onde

#Paramètres discrétisation de l'espace - maillage spatiale - indice i
L=0.58  #Longueur de la corde en m
dx=0.01 #pas entre deux points de l'espace
Nx=int(L/dx) #Nombre de points de l'espace
X=linspace(0,L,Nx+1)

#Paramètre de discrétisation du temps -maillage temporel - indice n
Duree=0.005  #Durée de la mesure de l'onde en seconde
dt=dx/c  #pas dans le temps
Nt=int(Duree/dt)  #Nombre de points dans l'espace
T=linspace(0,Duree,Nt+1)

#Conditions aux limites
u0l=0     
unl=0

#Parametre de calcul
alpha =(c*dt/dx)    
MODE=2 #Nombre de modes spatiaux

#Création de la matrice triagonale D
D=M(Nx)

#Création de la matrice U des résultats
U=zeros((Nx+1,Nt+1),float)

#Conditions initiales
U0=init_U0(Nx)
U1=init_U1(Nx)
U[:,[0]]=U0
U[:,[1]]=U1


#Matrices des conditions limites
Cl=zeros((Nx+1,1),float)
Cl[0]=u0l
Cl[Nx]=unl

#Calcul du reste des valeurs
for j in range(1,Nt):
  U[:,[j+1]]=dot(D,U[:,[j]])-U[:,[j-1]]+Cl

# Affichage de la solution Sans animation
fig = plt.figure(figsize=(7,4))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)

ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,U,cmap = 'viridis')       
plt.show()

#Graphique 2D en fonction de la position
for i in range(0,10):
  plt.plot(X,U[:,[i]])

plt.xlabel("position")
plt.ylabel("amplitude")
plt.title("Graphique 2D en fonction de la positon et à t fixé")
plt.grid("True")
plt.show()

#Graphique 2D en fonction du temps
for i in range(0,10): 
  plt.plot(T,U[i,:])

plt.xlabel("Temps")
plt.ylabel("amplitude")
plt.title("Graphique 2D en fonction du temp et à x fixé ")
plt.grid("True")

plt.show()