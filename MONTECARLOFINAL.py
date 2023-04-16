#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import math as m
import numpy as np
import matplotlib.pyplot as plt
import random as r
import time


thermalization = 1000	# number of sweeps to equilibrate at each T
num_mc_steps = 10000    # number of Monte Carlo sweeps for each T

cnt_L = 0

# Tc = 2/ln(1+sqrt(2)) = 2.269185314213 ..

Temperatures = [4, 3.5, 3, 2.7, 2.5, 2.4, 2.35, 2.3, 2.29, 2.28, 2.27, 2.2, 2.25, 2.24, 2.23, 2.22, 2.21, 2.2, 2.1, 2.0, 1.8, 1.5, 1, 0.5]

# Testing purposes:
# Sizes = [5,10,20]
Sizes = [22,50,70]

M=np.zeros((len(Sizes),len(Temperatures)))      # initialize a matrix where we insert the values of the magnetization |M|
M2=np.zeros((len(Sizes),len(Temperatures)))     # initialize a matrix where we insert the values of M^2
M4=np.zeros((len(Sizes),len(Temperatures)))     # initialize a matrix where we insert the values of M^4
En=np.zeros((len(Sizes),len(Temperatures)))     # initialize a matrix where we insert the values of the energy
En2=np.zeros((len(Sizes),len(Temperatures)))    # initialize a matrix where we insert the values of the energy square

UL=np.zeros((len(Sizes),len(Temperatures)))     # initialize a matrix for Binder cumulant UL
C=np.zeros((len(Sizes),len(Temperatures)))      # initialize a matrix fpr specific heat C

# a seep

def balayage(reseau):
    for i in range(L):
        for j in range(L):
            voisins = []
            if (i % 2 == 0 and j % 2 == 0) or (i % 2 == 1 and j % 2 == 1):  # arretes
                voisins = [(i, (j + 1) % L), ((i + 1) % L, j), (i, (j - 1) % L), ((i - 1) % L, j)]
            else:  # milieux des côtés
                voisins = [((i - 1) % L, (j - 1) % L), ((i + 1) % L, (j + 1) % L), ((i - 1) % L, (j + 1) % L), ((i + 1) % L, (j - 1) % L)]
            dE = 2*reseau[i][j]*sum([reseau[x][y] for x, y in voisins])
            if dE < 0 or r.random() < m.exp(-dE/T):
                reseau[i][j] = -reseau[i][j]
    return reseau
                
# measure magnetization

def aimantation(reseau):
        return(np.sum(reseau))


# measure energy

def energie(reseau):
    E = 0
    for i in range(L):
        for j in range(L):
            if (i % 2 == 0 and j % 2 == 0) or (i % 2 == 1 and j % 2 == 1):  
                E += -reseau[i][j]*(reseau[(i+1)%L][j]+reseau[i][(j+1)%L])
    return E

# do simulations

for L in Sizes:
	reseau=[[np.sign(-1 + 2 *r.random()) for i in range(L)] for i in range(L)]   # initialize a with random spins
	cnt_T =0
	for T in Temperatures:
		print("L =", L, ", T =", T)

# Thermalize (required for each T)
		print("thermalization")
		for mc_step in range(thermalization):
			balayage(reseau)

# Start Monte Carlo simulation with calculation of M, U, and C
		print("measuring")
		moy, moy2, moy4 = 0, 0, 0
		moyE, moyE2 = 0, 0
		for mc_step in range(num_mc_steps):
			balayage(reseau)
			Mact = aimantation(reseau)/L/L  #L au carré = N
			moy += abs(Mact)
			moy2 += Mact**2
			moy4 += Mact**4
			Eact = energie(reseau)
			moyE += Eact
			moyE2 += Eact**2
    #if mc_step % 10 ==0:
				#plt.imshow(reseau)
				#plt.pause(0.1)
				#plt.clf()
		moy = moy/num_mc_steps
		moy2 = moy2/num_mc_steps
		moy4 = moy4/num_mc_steps
		moyE = moyE/num_mc_steps
		moyE2 = moyE2/num_mc_steps
		print("|M| =", moy, ", E =", moyE)
		M[cnt_L][cnt_T] = moy #pour chaque taille et chaque temp
		M2[cnt_L][cnt_T] = moy2 
		M4[cnt_L][cnt_T] = moy4 #cumulant de blander
		En[cnt_L][cnt_T] = moyE
		En2[cnt_L][cnt_T] = moyE2 #pour la chaleur specifqiue 
		print("# =======")

		cnt_T += 1
		
	cnt_L += 1
  
# magnetization per site

plt.subplot(221)
line1, = plt.plot(Temperatures, M[0], label = "L = " + str(Sizes[0]))  
line2, = plt.plot(Temperatures, M[1], label = "L = " + str(Sizes[1])) 
line3, = plt.plot(Temperatures, M[2], label = "L = " + str(Sizes[2]))
plt.legend(handles = [line1,line2,line3])
plt.xlabel('T')
plt.ylabel('<|M|>')	

# Binder cumulant

for L in range(len(Sizes)):
        for T in range(len(Temperatures)):
                UL[L][T] = 1-M4[L][T]/(3*M2[L][T]**2)
plt.subplot(222)
line1, = plt.plot(Temperatures, UL[0], label = "L = " + str(Sizes[0]))
line2, = plt.plot(Temperatures, UL[1], label = "L = " + str(Sizes[1]))
line3, = plt.plot(Temperatures, UL[2], label = "L = " + str(Sizes[2]))
plt.legend(handles = [line1,line2,line3])
plt.xlabel('T')
plt.ylabel('B')

# Specific heat per site

for L in range(len(Sizes)):
        for T in range(len(Temperatures)):
                C[L][T] = (En2[L][T] - En[L][T]**2)/(Temperatures[T]*Sizes[L])**2
plt.subplot(223)
line1, = plt.plot(Temperatures, C[0], label = "L = " + str(Sizes[0]))
line2, = plt.plot(Temperatures, C[1], label = "L = " + str(Sizes[1]))
line3, = plt.plot(Temperatures, C[2], label = "L = " + str(Sizes[2]))
plt.legend(handles = [line1,line2,line3])
plt.xlabel('T')
plt.ylabel('U_L')
plt.show()		



