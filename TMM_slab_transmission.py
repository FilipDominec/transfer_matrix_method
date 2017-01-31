#!/usr/bin/env python
#coding:utf8

from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm


## settings
polarisation =  "TM"        ## also known as "p-polarisation"
#polarisation =  "TE"        ## AKA "s-polarisation"
plot_type =     "amplitude"
#plot_type =     "phase"

k       = 1                         # wavenumber
beta_0   = arange(0,pi/2,.01)       # incidence angle [rad]
Z_t = 1.
Z_0 = 1.

eps1s   = arange(-2, 8,.05)         # relative permittivity
#eps1s = 2.

mu1     = 1
d1      = 4                         # thickness of the slab
## /settings


t = []
aux = []
aux2 = []
for eps1 in eps1s:
    N1      = sqrt(0j+eps1*mu1)     # index of refraction of the slab (todo: rewrite to eps,mu)
    Z1      = sqrt(0j+mu1/eps1)     # impedance of the slab

    cos_beta_1   = sqrt(1+0j - (1/N1 * sin(beta_0))**2)
    delta_1  = k*N1*d1*cos_beta_1
    if polarisation == "TE":
        gamma_0 = cos(beta_0)
        gamma_t = cos(beta_0)
        gamma_1 = cos_beta_1/Z1
    elif polarisation == "TM":
        gamma_0 = 1/cos(beta_0)
        gamma_t = 1/cos(beta_0)
        gamma_1 = 1/(cos_beta_1*Z1)

    M_1 = array([[cos(delta_1),                1j/gamma_1*sin(delta_1)   ],
               [1j*gamma_1*sin(delta_1),    cos(delta_1)             ]])
    M = M_1     # compute the transfer matrix of whole structure

    if polarisation == "TE":
        t_TE = 2*gamma_0 / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1])
        t.append(t_TE)
        aux.append(delta_1/pi)
    elif polarisation == "TM":
        t_TM = 2*gamma_t / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1]) * (Z_t / Z_0)
        t.append(t_TM)
        aux.append(delta_1/pi)

plt.figure(figsize=[8,8])

if plot_type == "amplitude":
    plt.contourf(beta_0, eps1s, abs(array(t)), 100, cmap=cm.hot) # transmission amplitude
elif plot_type == "phase":
    plt.contourf(beta_0, eps1s, angle(array(t)), 100, cmap=cm.hsv) # transmission phase
    plt.contour(beta_0, eps1s, abs(array(t)), [.998], colors='k') # transmission amplitude high
    plt.contour(beta_0, eps1s, abs(array(t)), [.002], colors='w') # transmission amplitude low

## For lossless dielectric, there are three sufficient conditions for 100% transmission:
## 1) Fabry-Perot resonance
plt.clabel(plt.contour(beta_0, eps1s, array(aux), range(1,10), colors='black'))
## 2) If mu=1.0 and epsilon=1.0:
plt.plot(beta_0, ones_like(beta_0), color='black', lw=2)
## 3) Brewster angle of incidence
if polarisation == "TM":
    plt.plot(arccos((1-eps1s)/(1+eps1s))/2, eps1s, color='black', lw=2)
    plt.plot(arcsin(((1-1/eps1s)/(1-eps1s**2))**.5), eps1s, color='red', lw=1)



## 2) Zero impedance mismatch
#plt.plot(beta_0, ones_like(beta_0), color='blue')
    
plt.xlabel('Angle of incidence [rad]')
plt.ylabel('Relative permittivity of the slab')
plt.grid()
#plt.show()
plt.savefig('%s_slab_%s_XX.png' % (polarisation, plot_type))
