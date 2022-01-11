###Code d'erreur pour graphes

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from pylab import*
import time as tm
#caractéristique physique 
c=340 #m.s^-1   célérité de l'onde

#Paramètres discrétisation de l'espace - maillage spatiale - indice i
L=0.5  #Longueur de la corde en m
dx=0.0005 #pas entre deux points de l'espace
Nx=int(L/dx) #Nombre de points de l'espace
X=linspace(0,L,Nx)
dxe=dx

#Paramètre de discrétisation du temps -maillage temporel - indice n
Duree=0.001  #Durée de la mesure de l'onde en seconde
dt=dx/c  #pas dans le temps
Nt=int(Duree//dt)  #Nombre de points dans l'espace
T=linspace(0,Duree,Nt)
dte=dt

#Conditions aux limites
u0l=0     
unl=0

#Parametre de calcul
alpha =(c*dt/dx)    
MODE=3 #Nombre de modes spatiaux

###Solution exacte
Teuler=T
Xeuler=X

def fonction(x,t):
    #Applique la fonction d'onde à un entier ou un array
    
    return cos((MODE*pi*c/L)*(t*dt))*sin((MODE*pi/L)*(x*dx))
    
#Création de la matrice U des résultats
SolE=zeros((Nx,Nt),float)

x=0
t=0
while (t<Nt):
    x=0
    while(x<Nx):
        SolE[x,t]=fonction(x,t)
        x+=1
    t+=1

#Affichage de la solution exacte
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")
ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,SolE,cmap = 'cividis')       
plt.show()
## Euler explicite


#Fonctions des conditions initiales
def y0x(i):
  """int -> float
  Fonction de condition initiale y(0,x)
  onde sinusoïdale"""
  return sin((MODE*pi/L)*(i*dx))

def y0x2(i):
  """int -> float
  Fonction de condition initiale
  onde pincée triangulaire"""

  if((i*dx)<=x0):
    return 2*(i*dx)/x0
  
  else:
    return 2*(L-i*dx)/(L-x0)


def y0x3(i):
    """int -> float
  Fonction de condition initiale
    Corde pincée triangulaire"""
    
    return 2*(i*dx)*(L-i*dx)/(x0*L)

def dy0x(i):
  """int ->float
  Fonction de condition initiale d(0,x)"""
  return 0

#Création de la matrice tridiagonale  
def M(n):
  """int -> array(n,n)
 Renvoie matrice diagonale de dimension n*n contenant les coefficients du systeme"""

  M = zeros((n,n),float)

  for i in range(0,n):
    if(i==0):                                   #pour la première ligne de la matrice
      M[i,i], M[i,i+1]=(2-2*alpha**2),alpha**2
    elif(i==n-1):
      M[i,i-1], M[i,i]=alpha**2,(2-2*alpha**2)
    else:                                       #Pour la dernière ligne de la matrice
      M[i,i-1], M[i,i], M[i,i+1] = alpha**2,(2-2*alpha**2),alpha**2

  return M

###Fonctions des conditions initiales de la matrices

#Initialise la matrice colonne avec les valeurs à l'état initial t=0
#Ici,on considère que U(x,0)=sin(x)
def init_U0(n,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n,1),float)

  for i in range(0,n):
    U0[i,0]=y(i)
  #U0[n-1,0]=0
  return U0


#Initialise la matrice colonne avec les valeurs à l'état t=1 
def init_U1(n,dy,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=1 """
  #Initialisation des tableaux
  U0=init_U0(n,y)
  U1=zeros((n,1),float)

  for i in range (1,n-1):
    U1[i]= dt*dy(i) + U0[i]
  #U1[0]=U1[n-1]=0.0 #Conditions initiales à modifier 
  
  return U1

##############################################################################################


#Création de la matrice triagonale D
D=M(Nx)

#Création de la matrice U des résultats
Uexp=zeros((Nx,Nt),float)

#Conditions initiales
U0=init_U0(Nx,y0x)
U1=init_U1(Nx,dy0x,y0x)
Uexp[:,[0]]=U0
Uexp[:,[1]]=U1


#Matrices des conditions limites
Cl=zeros((Nx,1),float)
Cl[0]=u0l
Cl[Nx-1]=unl

#Calcul du reste des valeurs
for j in range(1,Nt-1):
  Uexp[:,[j+1]]=dot(D,Uexp[:,[j]])-Uexp[:,[j-1]]+Cl


#Affichage de la solution Sans animation
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")
ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,Uexp,cmap = 'viridis')       
plt.show()

#Euler Implicite

#Fonctions des conditions initiales
def y0x(i):
  """int -> float
  Fonction de condition initiale y(0,x)
  onde sinusoïdale"""
  return sin((MODE*pi/L)*(i*dx))

def y0x2(i):
  """int -> float
  Fonction de condition initiale
  onde pincée triangulaire"""

  if((i*dx)<=x0):
    return 2*(i*dx)/x0
  
  else:
    return 2*(L-i*dx)/(L-x0)


def y0x3(i):
    """int -> float
  Fonction de condition initiale
    Corde pincée triangulaire"""
    
    return 2*(i*dx)*(L-i*dx)/(x0*L)

def dy0x(i):
  """int ->float
  Fonction de condition initiale d(0,x)"""
  return 0

#Création de la matrice A inversible
def A(n):
  """int -> array(n,n)
 Renvoie matrice diagonale de dimension n*n contenant les coefficients du systeme"""

  A = zeros((n,n),float)

  for i in range(0,n):
    if(i==0):                                   #pour la première ligne de la matrice
      A[i,i], A[i,i+1]=(alpha**2)+1/2,-(alpha**2)/2
    elif(i==n-1):                               #Pour la dernière ligne de la matrice
      A[i,i-1], A[i,i]=-(alpha**2)/2,(alpha**2)+1/2
    else:                                       
      A[i,i-1], A[i,i], A[i,i+1] = -(alpha**2)/2,(alpha**2)+1/2,-(alpha**2)/2
  return A


#Création de la matrice U0 à valeur initiale
def init_U0(n,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n,1),float)

  for i in range(0,n):
    U0[i,0]=y(i)
  U0[n-1,0]=0
  return U0

#Initialise la matrice colonne avec les valeurs à l'état t=1 
def init_U1(n,dy,y):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=1 """
  #Initialisation des tableaux
  U0=init_U0(n,y)
  U1=zeros((n,1),float)

  for i in range (1,n-1):
    U1[i]= dt*dy(i) + U0[i]
  U1[0]=U1[n-1]=0.0 #Conditions initiales à modifier 
  
  return U1

##########################################################################################

#Création de la matrice U des résultats
Uimp=zeros((Nx,Nt),float)

#Conditions initiales
Uimp[:,[0]]=init_U0(Nx,y0x)
Uimp[:,[1]]=init_U1(Nx,dy0x,y0x)

# Inverse de la matrice A
InvA = linalg.inv(A(Nx))

#Calcul du reste des valeurs
for j in range(1,Nt-1):
  Sum = Uimp[:,[j]] - 1/2 * Uimp[:,[j-1]]
  Uimp[:,[j+1]] = dot(InvA, Sum)
    
# GRAPHIQUE AMPLITUDE EN FONCTION DU TEMPS, EN FONCTION DE LA POSITION
fig = plt.figure(figsize=(14,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")
ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,Uimp,cmap = 'Blues')       
plt.show()


"""Parametres de discrétisation/Résolution"""
#caractéristique physique 
c=340 #m.s^-1   célérité de l'onde

#Paramètres discrétisation de l'espace - maillage spatiale - indice i
L=0.5 #Longueur de la corde en m
dx=0.005 #pas entre deux points de l'espace
Nx=int(L/dx) #Nombre de points de l'espace
X=linspace(0,L,Nx)

#Paramètre de discrétisation du temps -maillage temporel - indice n
Duree=0.001  #Durée de la mesure de l'onde en seconde
dt=0.000001  #pas dans le temps
Nt=int(Duree/dt)  #Nombre de points dans l'espace
T=linspace(0,Duree,Nt)

#Caractéristique spec
alpha=(c/dx)**2

#Conditions aux limites
u0l=0     
unl=0

MODE=3 #Nombre de modes spatiaux

x0=L/2

#########################################################################################

"""Initialisation des conditions initiales"""
#Fonctions des conditions initiales
def y0x(i):
  """int -> float
  Fonction de condition initiale y(0,x)
  onde sinusoïdale"""
  return sin((MODE*pi/L)*(i*dx))

def y0x2(i):
  """int -> float
  Fonction de condition initiale
  onde pincée triangulaire"""

  if((i*dx)<=x0):
    return 2*(i*dx)/x0
  
  else:
    return 2*(L-i*dx)/(L-x0)


def y0x3(i):
    """int -> float
  Fonction de condition initiale
    Corde pincée triangulaire"""
    
    return 2*(i*dx)*(L-i*dx)/(x0*L)


#Initialise la matrice colonne avec les valeurs à l'état initial t=0
#Ici,on considère que U(x,0)=sin(x)
def init_U0(n):
  """int -> array(n,1)
  Retourne une matrice colonne composée des valeurs de l'onde à l'état initial t=0 """

  U0=zeros((n,1),float)

  for i in range(0,n):
    U0[i,0]=y0x(i)
  U0[n-1,0]=0
  return U0
  

#Initialise la matrice colonne avec les conditions initiales de dU/dt en t=0
def init_dU0(n):
    """int -> array(n,1)
    Retourne une matrice colonne composée des valeurs de la dérivée de la fonction d'onde à l'état initial t=0"""
    
    dU0=zeros((n,1),float)
    
    for i in range(0,n):
        dU0[i,0]=0
    
    #U0[n-1,0]=0
    return dU0


"""Définition des fonctions f et g"""
def f(n,i,U):
    """int + int + array(n,1) -> float
    Retourne la valeur de f(n,i,U)"""
    
    if(i==0):
        return alpha* (u0l - 2*U[i,0] + U[i+1,0])
    
    elif(i==n-1):
        return alpha* (U[i-1,0] - 2*U[i,0] + unl)
        
    else:
        return alpha* (U[i-1,0] - 2*U[i,0] + U[i+1,0] )
    
def g(z):
    """float -> float
    Retourne la valeur de g(n,i,U)=g(z)=z"""
    
    return z
 
 
#FONCTIONS F ET G
def K0i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*Z[i,0]
    
def L0i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*f(Nx,i,U)
    
def K1i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + 1/2*L0i(i,U,Z))
    
def L1i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l - 2*(U[i,0] +1/2*K0i(i,U,Z))  + (U[i+1,0]+1/2*K0i(i+1,U,Z)))
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+1/2*K0i(i-1,U,Z)) - 2*(U[i,0]+1/2*K0i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+1/2*K0i(i-1,U,Z)) - 2*(U[i,0]+1/2*K0i(i,U,Z)) + (U[i+1,0]+1/2*K0i(i+1,U,Z)) )
    
    return dt*T
   
def K2i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + 1/2*L1i(i,U,Z))

def L2i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l - 2*(U[i,0] +1/2*K1i(i,U,Z))  + (U[i+1,0]+1/2*K1i(i+1,U,Z)) )
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+1/2*K1i(i-1,U,Z)) - 2*(U[i,0]+1/2*K1i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+1/2*K1i(i-1,U,Z)) - 2*(U[i,0]+1/2*K1i(i,U,Z)) + (U[i+1,0]+1/2*K1i(i+1,U,Z)) )
    
    return dt*T
    
def K3i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    
    return dt*(Z[i,0] + L2i(i,U,Z) )   
 
def L3i(i,U,Z):
    """int + array(Nx,1) + array(Nx,1)-> float"""
    if(i==0):
        T=alpha* (u0l+ - 2*(U[i,0] +K1i(i,U,Z))  + (U[i+1,0]+K1i(i+1,U,Z)) )
        
    elif(i==Nx-1):
        T= alpha* ((U[i-1,0]+K2i(i-1,U,Z)) - 2*(U[i,0]+K2i(i,U,Z)) + unl)
        
    else:
        T= alpha* ((U[i-1,0]+K2i(i-1,U,Z)) - 2*(U[i,0]+K2i(i,U,Z)) + (U[i+1,0]+K2i(i+1,U,Z)) )
    
    return dt*T

#######################################################################################################

start=tm.process_time()
 
#Création de la matrice résultat
U=zeros((Nx,Nt),float)
U[:,[0]]=init_U0(Nx)

z=zeros((Nx,Nt),float)
z[:,[0]]=init_dU0(Nx)

#Paramètres d'itération
t=1
i=0

print(Nx)
print(Nt)
tm.sleep(1)

while(t<Nt):
    i=0
    
    while(i<Nx):
        #calculs de k0 et l0
        k0,l0=K0i(i,U[:,[t-1]],z[:,[t-1]]) , L0i(i,U[:,[t-1]],z[:,[t-1]])
        #print("k0=",k0,"l0=",l0)
        #calculs de k1 et l1
        k1,l1=K1i(i,U[:,[t-1]],z[:,[t-1]]), L1i(i,U[:,[t-1]],z[:,[t-1]])

        #calculs de k2 et l2
        k2,l2=K2i(i,U[:,[t-1]],z[:,[t-1]]), L2i(i,U[:,[t-1]],z[:,[t-1]])
        
        #calculs de k3 et l3
        k3,l3=K3i(i,U[:,[t-1]],z[:,[t-1]]), L3i(i,U[:,[t-1]],z[:,[t-1]])
        #print("k3=",k3,"l3=",l3)
        #Calculs de u et de z pour i et t
        U[i,t]=U[i,t-1]+ 1/6*(k0 + 2*k1 + 2*k2 + k3)
        z[i]=z[i,t-1]+ 1/6*(l0 + 2*l1 + 2*l2 + l3)
        
        i+=1
        #pdb.set_trace()
    print("t=",t)
    t+=1
     

end=tm.process_time()
print("Le temps d'exécution est de ",end-start,"secondes")

#############################################################################################

#Affichage de la solution Sans animation
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Amplitude en fonction de x et t")

ST,SX = meshgrid(T,X)
p = ax.plot_surface(SX,ST,U,cmap = 'plasma')       
plt.show()


Trk=T
Xrk=X
###Solution exacte
def fonction(x,t):
    #Applique la fonction d'onde à un entier ou un array
    
    return cos((MODE*pi*c/L)*(t*dt))*sin((MODE*pi/L)*(x*dx))
    
#Création de la matrice U des résultats
Solrk=zeros((Nx,Nt),float)

x=0
t=0
while (t<Nt):
    x=0
    while(x<Nx):
        Solrk[x,t]=fonction(x,t)
        x+=1
    t+=1
    
#Calculs d'erreurs
Eexp=np.abs(Uexp-SolE)
Eimp=np.abs(Uimp-SolE)
Erk=np.abs(U-Solrk)

#Affichage de l'Erreur
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Erreur en fonction de x et t pour Euler explicite")
ST,SX = meshgrid(Teuler,Xeuler)
p = ax.plot_surface(SX,ST,Eexp,cmap = 'viridis')       
plt.show()

#Affichage de l'Erreur
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Erreur en fonction de x et t pour Euler implicite")
ST,SX = meshgrid(Teuler,Xeuler)
p = ax.plot_surface(SX,ST,Eimp,cmap = 'Blues')       
plt.show()

#Affichage de l'Erreur
fig = plt.figure(figsize=(15,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('position', fontsize = 16)
ax.set_ylabel('temps', fontsize = 16)
ax.set_zlabel('amplitude', fontsize = 16)
ax.view_init(elev = 15, azim = 120)
plt.title("Erreur en fonction de x et t pour Runge-Kutta")
ST,SX = meshgrid(Trk,Xrk)
p = ax.plot_surface(SX,ST,Erk,cmap = 'plasma')       
plt.show()



##Moyennes des erreurs pour chaque t pour Runge-Kutta
EE=zeros(Erk.shape[1])
N=np.arange(0,Erk.shape[1]*dt-dt,dt)

for i in range(0,Erk.shape[1]):
  EE[i]=mean(Erk[:,[i]])

plt.loglog(N,EE,label="Runge-Kutta")

##Moyennes des erreurs pour chaque t pour Euler explicite
EE=zeros(Eexp.shape[1])
N=np.arange(0,Eexp.shape[1]*dte,dte)

for i in range(0,Eexp.shape[1]):
  EE[i]=mean(Eexp[:,[i]])

plt.loglog(N,EE,label="Euler explicite")

EE=zeros(Eimp.shape[1])

#Moyennes des erreurs pour chaque t pour Euler implicite
for i in range(0,Eimp.shape[1]):
  EE[i]=mean(Eimp[:,[i]])
  
plt.loglog(N,EE,label="Euler implicite")
plt.xlabel("temps")
plt.ylabel("erreur")
plt.title("Moyennes des erreurs pour chaque t pour les 3 méthodes")
plt.grid("True")
plt.legend()
plt.show()


