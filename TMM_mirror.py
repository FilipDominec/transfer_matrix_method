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

#ks       = arange(0,4*pi,.003)                         # wavenumber
ks = [3.4, 3.5, 3.7, 3.8]

#beta_0  = arange(0,pi/2,.01)       # incidence angle [rad]
beta_0  = 0       # incidence angle [rad]
Z_t = 1.
Z_0 = 1.

eps1    = 1.46**2
mu1     = 1
d1      = .25/eps1**.5                         # thickness of the slab (quarter-wave plate)

eps2    = 6.40**2
mu2     = 1
d2      = .25/eps2**.5                         # thickness of the slab

## /settings


t = []      ## store transmission results
r = []      ## store reflection results
deltas=[]       ## For assessing the bandgap width
t_product=[]
debug = []
for k in ks:
    N1      = sqrt(0j+eps1*mu1)     # index of refraction of the slab for SiO2
    Z1      = sqrt(0j+mu1/eps1)     # impedance of the slab
    N2      = sqrt(0j+eps2*mu2)     # index of refraction for TiO2
    Z2      = sqrt(0j+mu2/eps2)     # impedance 

    cos_beta_1   = sqrt(1+0j - (1/N1 * sin(beta_0))**2)
    delta_1  = k*N1*d1*cos_beta_1
    cos_beta_2   = sqrt(1+0j - (1/N2 * sin(beta_0))**2)
    delta_2  = k*N2*d2*cos_beta_2
    deltas.append(delta_1)
    if polarisation == "TE":
        gamma_0 = cos(beta_0)
        gamma_t = cos(beta_0)
        gamma_1 = cos_beta_1/Z1
        gamma_2 = cos_beta_2/Z2
    elif polarisation == "TM":
        gamma_0 = 1/cos(beta_0)
        gamma_t = 1/cos(beta_0)
        gamma_1 = 1/(cos_beta_1*Z1)
        gamma_2 = 1/(cos_beta_2*Z2)

    ## Transfer matrices for each layer
    M_1 = array([[cos(delta_1),                1j/gamma_1*sin(delta_1)   ],
               [1j*gamma_1*sin(delta_1),    cos(delta_1)             ]])
    M_2 = array([[cos(delta_2),                1j/gamma_2*sin(delta_2)   ],
               [1j*gamma_2*sin(delta_2),    cos(delta_2)             ]])

    ## Build the transfer matrix for whole structure
    M = reduce(lambda x, y: dot(x, y), [dot(M_1, M_2)]*20)
    import scipy.linalg as la
    #print M, dot(M,M), la.det(M)
    #print
    A = dot(M_1, M_2)
    #print mpow(A,1000); print abs(la.eig(A)[0]); print A[0,0]+A[1,1]
    print la.eig(A)[0], 'Tr=', A[0,0]+A[1,1]

    ## Calculate and store reflection and transmission
    r.append(
        (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] - M[1,0] - gamma_t*M[1,1]) / \
        (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1])
        )
    if polarisation == "TE":
        t_TE = 2*gamma_0 / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1])
        t.append(t_TE)
    elif polarisation == "TM":
        t_TM = 2*gamma_t / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1]) * (Z_t / Z_0)
        t.append(t_TM)

    ## For assessing the bandgap width
    t_12 = 2*gamma_1 / (gamma_1*M_1[0,0] + gamma_1*gamma_2*M_1[0,1] + M_1[1,0] + gamma_2*M_1[1,1])
    t_21 = 2*gamma_2 / (gamma_2*M_2[0,0] + gamma_2*gamma_1*M_2[0,1] + M_2[1,0] + gamma_1*M_2[1,1])
    t_product.append(sqrt(abs((t_12 * t_21))))
    #debug.append(M[0,0]+M[1,1])
    #debug.append(abs(2-(2+gamma_1/gamma_2+gamma_2/gamma_1)*sin(delta_1)**2)) ## condition 1
    #debug.append(abs(2*cos(delta_1)**2 -(gamma_1/gamma_2+gamma_2/gamma_1)*sin(delta_1)**2)) 

    debug.append(abs(4*gamma_1*gamma_2/(gamma_1+gamma_2)**2) - sin(delta_1)**2)  ## condition 2
    #debug.append(abs((2+gamma_1/gamma_2 + gamma_2/gamma_1)*sin(delta_1)**2-2))  
plt.figure(figsize=[8,8])

#if plot_type == "amplitude":
    #plt.contourf(beta_0, ks, log(1-abs(array(t))), 100, cmap=cm.hot) # transmission amplitude
#elif plot_type == "phase":
    #plt.contourf(beta_0, ks, angle(array(t)), 100, cmap=cm.hsv) # transmission phase
    #plt.contour(beta_0, ks, abs(array(t)), [.998], colors='k') # transmission amplitude high
    #plt.contour(beta_0, ks, abs(array(t)), [.002], colors='w') # transmission amplitude low

## Plot reflection and transmission
t = array(t)
r = array(r)
plt.plot(ks, abs(t), color='blue', label='$t$')
plt.plot(ks, abs(r), color='red', label='$r$')
#plt.plot(ks, angle(array(t)), color='blue', ls='--', label='\\phi_t') ## phase
#plt.plot(ks, angle(array(r))/pi, color='red', ls='--', label='\\phi_r')

## For assessing the bandgap width
plt.plot(ks, t_product, color='green', label='$\\sqrt{t_{12} t_{21}}$')
plt.plot(ks, sin(deltas), color='green', label='$\\sin(\\delta_1) = \\sin(\\delta_2)$')

#plt.plot(ks, log(array(debug))/log(2), color='k', label='debug')
plt.plot(ks, debug, color='k', label='debug')

    
plt.xlabel('Wavenumber [rad]')
plt.ylabel('Amplitude')
plt.grid()
plt.legend()
plt.show()
#plt.savefig('%s_slab_%s_layers=100c.png' % (polarisation, plot_type))
