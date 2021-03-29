#%%

from numpy import linspace, zeros, load, linalg
from numpy import ndarray
import plotly.express as px
from math import nan, isnan


N = 30
R = 0.01

t = 180
Nt = 10000

dr = R / N
dt = t / Nt



def sq(x):
    return x*x

#%%
totalDT = []
with open("save_thermique.npy", "rb") as f:
    totalDT = load(f)


#%%

p = 0

E = 71e9
nuc = 0.21
kc = E/(3*(1-2*nuc))
muc= E/(2*(1+nuc))
etac = 1e5


alpha_c=9e-6

def k(i) -> float:
    return kc

def mu(i) -> float:
    return muc

def nu(i) ->float:
    return nuc


def eta(i) -> float:
    return etac

def alpha(i) -> float:
    return alpha_c


#%%

temps = Nt-1
DT = totalDT[temps]

M = [[0 for i in range(2*N)] for j in range(2*N)]
Rmat = [0]*(2*N)

M[0][1] = 1

M[2*N-1][2*N-2] = "-3*k(N)"
M[2*N-1][2*N-1] = "-4*mu(4)/((dr*N)**3)"
Rmat[2*N-1] = "p - 3*k(N)*alpha(N)*DT[N-1]"

M[2*N-1][2*N-2] = -3*k(N)
M[2*N-1][2*N-1] = 4*mu(4)/((dr*N)**3)
Rmat[2*N-1] = p - 3*k(N)*alpha(N)*DT[N-1]



for i in range(1, N):
    r = dr * i
    M[ 2*(i-1) + 1][ 2*(i-1) : 2*(i-1) + 4] = [-r, -1/sq(r), r, 1/sq(r)]
    #M[ 2*(i-1) + 1][ 2*(i-1) : 2*(i-1) + 4] = ["-r"+str(i), "-1/r"+str(i)+"²","r"+str(i), "1/r"+str(i)+"²"]

for i in range(1,N):
    r = dr*i
    M[2*(i-1)+2][2*(i-1) : 2*(i-1) + 4] = [-3*k(i), 4*mu(i)/(r**3), 3*k(i+1), -4*mu(i+1)/(r**3)]
    #M[2*(i-1)+2][2*(i-1) : 2*(i-1) + 4] = ["-3k"+str(i), "4mu"+str(i)+"/(r**3)", "3k"+str(i+1), "-4mu"+str(i+1)+"/(r**3)"]

    Rmat[2*(i-1)+2] = -3*k(i)*alpha(i)*DT[i-1]+3*k(i+1)*alpha(i+1)*DT[i]
    #Rmat[2*(i-1)+2] = "-3k"+str(i)+"*alpha"+str(i)+"*DT"+str(i)+"+3k"+str(i+1)+"*alpha"+str(i+1)+"*DT"+str(i)

#%%


def Tcontraintesi(sol_ab,i,r) -> ndarray:
    '''Calcule de la contrainte a partir des solution
    Attention i>0'''

    sigmai=zeros((3,3))
    ai,bi=sol_ab[2*i-2],sol_ab[2*i-1]
    sigmai[0][0]= -4*mu(i)*bi/(r**3) + 3*k(i)*ai - 3*k(i)*alpha(i)*DT[i-1]
    sigmai[1][1]=  2*mu(i)*bi/(r**3) + 3*k(i)*ai - 3*k(i)*alpha(i)*DT[i-1]
    sigmai[2][2]=  2*mu(i)*bi/(r**3) + 3*k(i)*ai - 3*k(i)*alpha(i)*DT[i-1]

    return sigmai



def Tdeplacementsi(sol_ab,i,r) -> ndarray:
    '''Calcule de la contrainte a partir des solution
    Attention i>0'''

    epsilonai=zeros((3,3))
    ai,bi=sol_ab[2*i-2],sol_ab[2*i-1]
    epsilonai[0][0]=ai - 2 * bi / (r**3)
    epsilonai[1][1]=ai + bi/(r**3)
    epsilonai[2][2]=ai + bi/(r**3)

    return epsilonai


def deplacement(sol_ab, i, r):
    ai,bi=sol_ab[2*i-2],sol_ab[2*i-1]
    return ai * r + bi/(r**2)



#%%
sol = linalg.solve(M, Rmat)

X = linspace(0, R, N)
Y1, Y2, Y3, Y4 = [], [], [], []

for i in range(1, N+1):
    r = i*dr
    contrainte = Tcontraintesi(sol, i, r)
    Y1.append(contrainte[0][0])
    Y2.append(contrainte[1][1])
    Y3.append(contrainte[2][2])
    
    Y4.append(deplacement(sol, i, r))

px.line(x = X, y = [Y1, Y2, Y3])
#px.line(x = X, y = Y4)


#%%


# X = linspace(0, R, N)
# Y1, Y2, Y3, Y4 = [], [], [], []


# totalContrainte = []

# for temps in range(0, Nt):
#     #M = M_t[temps]
#     sol = sol_t[temps]
#     c = []
#     for i in range(N):
#         r = i*dr
#         contrainte = Tcontraintesi(sol, i, r, temps)
        
#         c.append(contrainte)

#     totalContrainte.append(c)

# #px.line(x = X, y = [Y1, Y2, Y3])
# #px.line(x = X, y = Y4)