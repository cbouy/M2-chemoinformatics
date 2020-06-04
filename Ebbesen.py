#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Created on 14 févr. 2017

@author: Ced
@coauthor: Phil
@source: Ebbesen, T. W. Interference of Hydrogen with the Electron Transfer to Colloidal Platinum Catalyst and Consequences for Photochemical Water Reduction. J. Phys. Chem. 88, 4131–4135 (1984).
'''

from scipy.integrate import odeint
import numpy
import matplotlib.pyplot as plt


# constantes de vitesse
k1 = 10**8
km1 = 10**2
k2 = 10**10
km2 = 10**2
k3 = 10**3
km3 = 10**16
k5 = 1.47*(10**6)
k6 = 1.4*(10**9)
k7 = 4.4*(10**9)
k8 = 10**7

# autres constantes
Iabs = 10**(-7)        # variable
phi_T = 1
phi_C = 1


# définition du système étudié 
def sys_cinetique(Z,t):
       
    # espèces intervenant dans le système
    Am = Z[0]
    A = Z[1]
    Pt = Z[2]
    Ptm = Z[3]
    Hp = Z[4]
    H2 = Z[5]
    PtH = Z[6]
    S_etoile = Z[7]
    Sp = Z[8]
    S = Z[9]   
    D = Z[10]
   
    # equadiff
    dAmdt = k6*phi_C*S_etoile*A - k7*Sp*Am - k1*Am*Pt + km1*Ptm*A
    dAdt = -k6*phi_C*S_etoile*A + k7*Sp*Am + k1*Am*Pt - km1*Ptm*A
    dPtdt = -k1*Am*Pt + km1*Ptm*A + 2*k3*(PtH**2) - 2*km3*H2*(Pt**2)
    dPtmdt = k1*Am*Pt - km1*Ptm*A - k2*Ptm*Hp + km2*PtH
    dHpdt = km2*PtH - k2*Ptm*Hp
    dH2dt = k3*(PtH**2) - km3*H2*(Pt**2)
    dPtHdt = k2*Ptm*Hp - km2*PtH - 2*k3*(PtH**2) + 2*km3*H2*(Pt**2)
    dS_etoiledt = Iabs*phi_T*S - k5*S_etoile - k6*phi_C*S_etoile*A
    dSpdt = k6*phi_C*S_etoile*A - k7*Sp*Am - k8*Sp*D
    dSdt = -Iabs*phi_T*S + k5*S_etoile + k7*Sp*Am + k8*Sp*D
    dDdt = -k8*Sp*D
    
    return [dAmdt, dAdt, dPtdt, dPtmdt, dHpdt, dH2dt, dPtHdt, dS_etoiledt, dSpdt, dSdt, dDdt]

  

# conditions initiales
Am0 = 0
A0 = 10**(-3)    # varie
Pt0 = 10**(-5)   # varie
Ptm0 = 0
Hp0 = 10**(-5)
H20 = 0
PtH0 = 0
S_etoile0 = 0
Sp0 = 0
S0 = 10**(-4)
D0 = 10**(-2)
Z0 = [Am0, A0, Pt0, Ptm0, Hp0, H20, PtH0, S_etoile0, Sp0, S0, D0]


# plage de temps calculée ( t_min, t_max, nb_de_points_calculés )
t = numpy.linspace(0, 2*10**8, 100000)  
# équivalent à t = numpy.logspace(1, 8, num=100000)  pour avoir directement en échelle log


# ODE
soln = odeint(sys_cinetique, Z0, t, rtol=1.49012*10**(-15), atol=1.49012*10**(-15), mxstep=5000000)



# graphes concentration en fonction du temps

plt.figure(1)    # utile si on veut plusieurs figures
plt.subplot(221)    # figure 1 est divisée en 4 graphes (2*2), ici le N°1
plt.plot(t,soln[:,0],color='blue', linestyle='solid', label='A-')   # (axe_x, axe_y, ...)
plt.title('Evolution de [A-]')
#plt.legend()            pas nécessaire ici car une seule courbe par graphe
#plt.xlabel("Temps en secondes")
plt.ylabel("Concentration en M")
plt.yscale('log')    # linear par défaut, autre option : log
plt.xscale('log')

plt.subplot(222)
plt.plot(t,soln[:,5], color='red', linestyle='solid', label='H2')
plt.title('Evolution de [H2]')
#plt.legend()
#plt.xlabel("Temps en secondes")
#plt.ylabel("Concentration en M")
plt.yscale('log')    # linear par défaut, autre option : log
plt.xscale('log')

plt.subplot(223)
plt.plot(t,soln[:,6],color='green', linestyle='solid', label='PtH')
plt.title('Evolution de [PtH]')
#plt.legend()
plt.xlabel("Temps en secondes")
plt.ylabel("Concentration en M")
plt.yscale('log')    # linear par défaut, autre option : log
plt.xscale('log')

plt.subplot(224)
plt.plot(t,soln[:,2],color='brown', linestyle='solid', label='Pt')
plt.title('Evolution de [Pt]')
#plt.legend()
plt.xlabel("Temps en secondes")
#plt.ylabel("Concentration en M")
plt.yscale('log')    # linear par défaut, autre option : log
plt.xscale('log')


plt.figure(2)
plt.plot(t,soln[:,0],color='blue', linestyle='solid', label='A-')   # (axe_x, axe_y, ...)
plt.plot(t,soln[:,5], color='red', linestyle='solid', label='H2')
plt.plot(t,soln[:,6],color='green', linestyle='solid', label='PtH')
#plt.plot(t,soln[:,3],color='orange', linestyle='solid', label='Pt-')
plt.plot(t,soln[:,2],color='brown', linestyle='solid', label='Pt')
plt.title('Evolution du système')
plt.legend()
plt.xlabel("Temps en secondes")
plt.ylabel("Concentration en M")
plt.yscale('log')    # linear par défaut, autre option : log
plt.xscale('log')

plt.show()
#plt.savefig('test.png')