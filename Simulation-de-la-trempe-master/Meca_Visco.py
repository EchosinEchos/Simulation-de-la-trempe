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

Tref = 900


def sq(x):
    return x*x

#%%

# with open("save_thermique.npy", "rb") as f:
#     T = load(f)


# totalDT = [ [Tref - T[m][i] for i in range(len(T[0])-0)] for m in range(len(T))]
totalDT = []
X = linspace(0, R, N)
for i in range(Nt):
    totalDT.append([ -(r*r/(R*R) * i/Nt)  for r in X] )
    #totalDT.append([0]*N)

#%%

p = 0

E = 71e9
nuc = 0.21
kc = E/(3*(1-2*nuc))
muc= E/(2*(1+nuc))
etac = 1e30


alpha_c=9e-6

def k(i, t) -> float:
    return kc

def mu(i, t) -> float:
    return muc

def nu(i, t) ->float:
    return nuc


def eta(i, t) -> float:
    return etac

def alpha(i, t) -> float:
    return alpha_c


#%%

sol_t = []

thorr = [ [0] * N for i in range(Nt)]
err = [ [0] * N for i in range(Nt)]



for temps in range (1, Nt):
    
    tho = eta(N-1, temps)/mu(N-1, temps)
    tho = dt*10 if dt *10 > tho else tho

    DT = totalDT[temps]

    M = [[0 for i in range(2*N)] for j in range(2*N)]
    Rmat = [0]*(2*N)

    M[0][1] = 1


    M[2*N-1][2*N-2] = -3*k(N, temps)
    M[2*N-1][2*N-1] = 4*mu(4, temps)/((dr*N)**3)
    Rmat[2*N-1] = p - 3*k(N, temps)*alpha(N, temps)*DT[N-1]  #- 2*mu(N, temps)*err[temps-1][N-2] + thorr[temps-1][N-2]*(1-dt/tho)



    for i in range(1, N):
        r = dr * i
        M[ 2*(i-1) + 1][ 2*(i-1) : 2*(i-1) + 4] = [-r, -1/sq(r), r, 1/sq(r)]

    for i in range(1,N):

        tho_prec = tho
        tho = eta(i, temps)/mu(i, temps)

        tho = dt*10 if dt*10 > tho else tho
        tho_prec = dt*10 if dt*10 > tho_prec else tho_prec


        r = dr*i
        M[2*(i-1)+2][2*(i-1) : 2*(i-1) + 4] = [-3*k(i, temps), 4*mu(i, temps)/(r**3), 3*k(i+1, temps), -4*mu(i+1, temps)/(r**3)]

        ajout_visco = - thorr[temps-1][i-1]*(1-dt/tho) + thorr[temps-1][i-2]*(1-dt/tho_prec) + 2*mu(i, temps)*err[temps-1][i-1] - 2*mu(i, temps)*err[temps-1][i-1]

        Rmat[2*(i-1)+2] = -3*k(i, temps)*alpha(i, temps)*DT[i-1]+3*k(i+1, temps)*alpha(i+1, temps)*DT[i-1] + ajout_visco*0

    
    sol = linalg.solve(M, Rmat)
    sol_t.append(sol)


    thorr[temps] = [-4*mu(i, temps)*sol[2*i+1]/((dr*i)**3) for i in range(1,N)]
    err[temps] = [-2*sol[2*i+1]/((dr*i)**3) for i in range(1,N)]


#%%


def Tcontraintesi(sol_ab,i,r, temps) -> ndarray:
    '''Calcule de la contrainte a partir des solution
    Attention i>0'''

    sigmai=zeros((3,3))
    ai,bi=sol_ab[2*i],sol_ab[2*i+1]
    sigmai[0][0]= -4*mu(i, temps)*bi/(r**3) + 3*k(i, temps)*ai - 3*k(i, temps)*alpha(i, temps)*DT[i]
    sigmai[1][1]=  2*mu(i, temps)*bi/(r**3) + 3*k(i, temps)*ai - 3*k(i, temps)*alpha(i, temps)*DT[i]
    sigmai[2][2]=  2*mu(i, temps)*bi/(r**3) + 3*k(i, temps)*ai - 3*k(i, temps)*alpha(i, temps)*DT[i]

    return sigmai



def Tdeplacementsi(sol_ab,i,r) -> ndarray:
    '''Calcule du déplacement à partir des solution
    Attention i>0'''

    epsilonai=zeros((3,3))
    ai,bi=sol_ab[2*i],sol_ab[2*i+1]
    epsilonai[0][0]=ai - 2 * bi / (r**3)
    epsilonai[1][1]=ai + bi/(r**3)
    epsilonai[2][2]=ai + bi/(r**3)

    return epsilonai


def deplacement(sol_ab, i, r):
    ai,bi=sol_ab[2*i],sol_ab[2*i+1]
    return ai * r + bi/(r**2)



#%%
X = linspace(0, R, N)
sol = sol_t[0]
c = []
for i in range(N):
    r = i*dr
    contrainte = Tcontraintesi(sol, i, r, 0)
    
    c.append(contrainte[0][0])

px.line(x = X, y = c)
#px.line(x = X, y = T[0])


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