#%%

import plotly.express as px
from numpy import linspace, save


#Initialisation des constante
lbda = 0.96
rho = 2550
Cp = 840

N = 30
R = 0.01

t = 180
Nt = 10000

dr = R / N
dt = t / Nt

Text = list(linspace(900,300, round(2/dt))) + [300]*(Nt-round(2/dt))
Tref = 900

######################################


#Calcule des constantes
W = lbda/(rho * Cp)
T = []

#Initialisation de la temperaure
for i in range(Nt):
    T.append( ([Tref]*N)  + [Text[i]] )


######################################
def sq(x):
    return x*x


#####################################


print("########### Calcule des T #############")
for m in range(1, Nt):

    T[m][0] = W * dt * (T[m-1][1] - T[m-1][0])/sq(dr) + T[m-1][0]

    for i in range(1,N):
        T[m][i] = W * dt * ( (T[m-1][i+1] - 2*T[m-1][i] + T[m-1][i-1])/sq(dr) + (4*(T[m-1][i+1] - T[m-1][i])) / ( sq(dr) * (2*i - 1))) + T[m-1][i]

totalDT = [ [Tref - T[m][i] for i in range(len(T[0])-1)] for m in range(len(T))]



########### Enregistrement dans un fichier ################


print("########### Enregistrement dans un fichier #############")


with open("save_thermique.csv", "w") as f:

    for i in range(len(totalDT)):
        for j in range( len(totalDT[i]) - 1):
            f.write(str(totalDT[i][j]) + ";")
        f.write(str(totalDT[i][-1]))
        f.write("\n")

    f.close()


with open("save_thermique.npy", "wb") as f:
    save(f, T)



#%%
# with open("donnes.csv", "w") as f:

#     for i in range(len(T)):
#         for j in range( len(T[i]) - 1):
#             f.write(str(T[i][j]) + ";")
#         f.write(str(T[i][-1]))
#         f.write("\n")

#     f.close()


###########################################################


#%%
X = linspace(0, R, N+1)
Y = T[200][:]

px.line(x = X , y = Y)
