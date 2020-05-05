from sympy import Eq,solve,var,symbols,Function
from pandas import read_csv
from numpy import size,clip
import numpy as np
import operator
from functools import reduce  # Required in Python 3
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

s = read_csv('stoichiometry',sep=' ')

s = s.to_numpy()
m = (-s).clip(0)

nReactions = size(s,0)
nComponents = size(s,1)

var('t')
C = []
for i in range(0,nComponents):
    C.append(Function('C'+str(i))(t))

#C = symbols('C0:%d'%nComponents)

k = symbols('k0:%d'%nReactions)

RL = [k[i]*prod([C[n]**m[i,n] for n in range(0,nComponents)]) for i in range(0,nReactions)]

r = []

for j in range(0,nComponents):
    R = 0
    for i in range(0,nReactions):
        R = R + RL[i]*s[i,j]
    r.append(R)

from sympy import Eq,solve,var,symbols
from pandas import read_csv
from numpy import size,clip
import numpy as np
import operator
from functools import reduce  # Required in Python 3
def prod(iterable):
    return reduce(operator.mul, iterable, 1)

s = read_csv('stoichiometry',sep=' ')
s = s.to_numpy()
m = (-s).clip(0)

nReactions = size(s,0)
nComponents = size(s,1)

var('t')
C = []
for i in range(0,nComponents):
    C.append(Function('C'+str(i))(t))

#C = symbols('C0:%d'%nComponents)

k = read_csv('rateConstant',sep=' ')
k = k.to_numpy()[0]

RL = [k[i]*prod([C[n]**m[i,n] for n in range(0,nComponents)]) for i in range(0,nReactions)]

r = []

for j in range(0,nComponents):
    R = 0
    for i in range(0,nReactions):
        R = R + RL[i]*s[i,j]
    r.append(R)

V = Function('V')(t)

theta = read_csv('theta',sep=' ')
theta = theta.to_numpy()[0]

dVdt = solve(Eq(sum(theta)*V.diff(t),sum(r)*V),V.diff(t))[0]

model = []

for i in range(0,nComponents):
    model.append(Eq(V*C[i].diff(t)+C[i]*dVdt,r[i]*V))
    #model.append(Eq(C[i].diff(t),r[i]))
    model[i] = solve(model[i],C[i].diff(t))[0]

totalTime = 10
N = 1000
dt = totalTime/N

c = np.empty((N,nComponents))
c[:] = np.NaN

for j in range(0,nComponents):
    c[0][j] = theta[j]

for i in range(0,N-1):
    for j in range(0,nComponents):
        rateValue = model[j]
        for M in range(0,nComponents):
            rateValue = rateValue.subs(C[M],c[i,M])
        rateValue = float(rateValue)
        c[i+1][j] = c[i][j] + rateValue*dt

np.savetxt('kineticSol',c)
