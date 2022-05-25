import numpy as np
import random
import matplotlib.pyplot as plt

plt.style.use(["science", "notebook", "grid"])

##==============================================================================

def SI_simulation(tmax, zom = 0.4, a = 0.2, d = 0.1, 
                  S0 = 999, Z0 = 1, D0 = 0, da = 0, dz = 0, dd = 0):
    '''
    tmax: [Only required variable] How many days the simulation is run for
    zom: initial chance a zombiie will infect a human
    a: chance a human will succefully permanently kill a zombie
    d: chance a zombie completely decays
    S0: Initial suscceptable population
    D0: Initial dead population
    Z0: Initial zombie population
    da: rate at which a changes
    dz: rate at which zom changes
    dd: rate at which d changes
    '''
    Nt = tmax + 1
    t = np.linspace(0,tmax,Nt)
    
    (S, Z, D) = (np.zeros(Nt, dtype = 'int'), 
                 np.zeros(Nt, dtype = 'int'), 
                 np.zeros(Nt, dtype = 'int'))
    
    S[0], Z[0], D[0] = S0, Z0, D0 # Initial variables
    
    for i in range(1, Nt):
        S[i] = S[i-1]
        Z[i] = Z[i-1]
        D[i] = D[i-1]

        sp = random.randint(0, S[i-1])
        zp = random.randint(0, Z[i-1])
        ref  = min(sp,zp)
        
        if i > 5:
            '''
            Everyday after day 5, A bunch of zombies will decay completely and
            permaanently die
            '''
            for j in range(Z[i-1]):
                p = random.random()
                if p <= d:
                    Z[i] -= 1
                    D[i] += 1
        
        for j in range(ref):
            '''
            Everyday, a couple of survivors will interact with a couple of zombies
            annd they will fight.
            '''
            p = random.random()
            q = random.random()
            if p <= zom:
                S[i] -= 1
                Z[i] += 1
            if np.all([q <= a, i > 3]):
                Z[i] -= 1
                D[i] += 1
        
        '''
        Everyday after the first five days, the survivors
        get better at killing zombies and avoiding death. The zommbies also decay
        faster
        '''
        if all([i>= 5, a < 0.75, d < 0.10, zom > 0.20]): #Limits how muuch the variables change
            a += da
            zom -= dz
            d += dd
            
        
    return t, S, Z, D

##==============================================================================

tmax = 250 # How many days the simualtions are run for 
nsims = 20 # How many simulations are run

z0 = 2 
s0 = 1_000 - z0
d0 = 0

# Important: S0 + Z0 + D0 should *always* be constant

Zrate = 0.60
Attack = 0.155
Decay = 0.001

dZrate = 0.001
dAttack = 0.001
dDecay = 0.0001

##==============================================================================

#Used for testing indivdual simulations.

'''data = SI_simulation(tmax, zom = Zrate, a = Attack, d = Decay, 
                         dz = dZrate, da = dAttack, dd = dDecay)

plt.plot(t, data[1], lw = 3, c = 'b', label = 'Susceptable')
plt.plot(t, data[2], lw = 3, c = 'r', label = 'Zombies')
plt.plot(t, data[3], lw = 3, c = 'k', label = 'Dead')
plt.plot(t, data[3]+data[2]+data[1], lw = 3, c = 'g', label = 'total')'''


##==============================================================================

fig, axs = plt.subplots(2)

axs[0].legend()

if dZrate+dAttack+dDecay == 0:
    title = f'Zombie apocalypse model with {nsims} simulations and their \
average.\n\
first {tmax} days [static parameters]'
else:
    title = f'Zombie apocalypse model with {nsims} simulations and their \
average. \n\
first {tmax} days [dynamic parameters]'

fig.suptitle(title)

##==============================================================================

(Ssims, Zsims, Dsims) = (np.zeros((tmax+1, nsims)), np.zeros((tmax+1, nsims)), 
                         np.zeros((tmax+1, nsims)))

# Runs the simmuation Nsims amount of times and plots alll of them

ref1 = 0
ref2 = 1e10*tmax
best_case = (0,0,0,0)
worst_case = (0,0,0,0)

for i in range(nsims):
    data = SI_simulation(tmax, zom = Zrate, a = Attack, d = Decay, 
                         dz = dZrate, da = dAttack, dd = dDecay,
                         S0 = s0, Z0 = z0, D0 = d0)
    
    if np.sum(data[1]) > ref1:
        ref1 = np.sum(data[1])
        best_case = data
    elif np.sum(data[1]) <= ref2:
        ref2 = np.sum(data[1])
        worst_case = data

    Ssims[:,i] = data[1]
    Zsims[:,i] = data[2]
    Dsims[:,i] = data[3]

S = np.transpose(Ssims)
Z = np.transpose(Zsims)
D = np.transpose(Dsims)
t = np.linspace(0, tmax, tmax+1)

width = 0.2

for i in range(nsims):
    axs[0].plot(t, S[i], c = 'b', lw=width)
    axs[0].plot(t, Z[i], c = 'r', lw=width)
    axs[0].plot(t, D[i], c = 'k', lw=width)
    #axs[0].plot(t, (D[i]+Z[i]+S[i]), c = 'g', lw=width)

def average_out(X):
    ans = []
    for i in X:
        ans.append(np.mean(i))
    return np.array(ans)

# The averages of the simulations

S1, Z1, D1 = average_out(Ssims), average_out(Zsims), average_out(Dsims)

lwidth = 3

axs[0].plot(t, S1, lw = lwidth, c = 'b', label = 'Susceptable')
axs[0].plot(t, Z1, lw = lwidth, c = 'r', label = 'Zombies')
axs[0].plot(t, D1, lw = lwidth, c = 'k', label = 'Dead')
#axs[0].plot(t, (S1 + Z1 + D1), c = 'g', lw = lwidth)

axs[0].grid()
axs[0].legend()

lwidth = 1.5

axs[1].plot(t, best_case[1], lw = lwidth, c = 'b')
axs[1].plot(t, best_case[2], lw = lwidth, c = 'r')
axs[1].plot(t, best_case[3], lw = lwidth, c = 'k')

axs[1].plot(t, worst_case[1], lw = lwidth, c = 'b', ls = '--')
axs[1].plot(t, worst_case[2], lw = lwidth, c = 'r', ls = '--')
axs[1].plot(t, worst_case[3], lw = lwidth, c = 'k', ls = '--')

axs[1].grid()
#axs[1].legend()
axs[1].set_title('Best Case scenario (solid) vs Worst Case scenario (Dashed)', 
                 fontdict = {'fontsize': 9})

##==============================================================================

s = f'Time [Days] \n {Zrate*100}% initial chance of human getting infected on interaction. \n\
{Attack * 100}% initial chance of a human killing a zombie on interaction \n\
{Decay*100}% initial chance of a zombie decaying completely (on a given day)'

axs[1].set_xlabel(s)

axs[1].set_ylabel('Number of People')
axs[0].set_ylabel('Number of People')

axs[1].grid()
axs[0].grid()


plt.subplots_adjust(left=0.15, right=0.9, top=0.85, bottom=0.155)
plt.show()