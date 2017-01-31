#!/usr/bin/env python
#coding:utf8

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


## ==== User settings ====
#polarisation =  "TM"        ## also known as "p-polarisation"
polarisation =  "TE"        ## AKA "s-polarisation"

plot_type =     "amplitude"
#plot_type =     "phase"


## Note: this script plots 2-D contours. Exactly two of the following quantities
## should be defined as an np.arrayG

freqs     = np.linspace(0,.6e12, 200)
beta_0s   = np.linspace(0,np.pi/2, 200)       # incidence angle [rad]
#eps1s   = arange(-2, 8,.05)        # relative permittivity
eps1    = 4.
mu1     = 1.
d1      = 350e-6               # thickness of the slab (normalized for THz frequency range, 

## impedance of outer medium (front side, back side)
(Z_t, Z_0) = 1., 1.            

                                    # as needed for interpolation)
## ==== </user settings> ====


t, r = [], []
aux = []
aux2 = []
k      = 2*np.pi*freqs / 3e8                          # wavenumber
for beta_0 in beta_0s:
    N1      = np.sqrt(0j+eps1*mu1)     # index of refraction of the slab
    Z1      = np.sqrt(0j+mu1/eps1)     # impedance of the slab

    cos_beta_1   = np.sqrt(1+0j - (1/N1 * np.sin(beta_0))**2)
    delta_1  = k*N1*d1*cos_beta_1
    if polarisation == "TE":
        gamma_0 = np.cos(beta_0)
        gamma_t = np.cos(beta_0)
        gamma_1 = cos_beta_1/Z1
    elif polarisation == "TM":
        gamma_0 = 1/np.cos(beta_0)
        gamma_t = 1/np.cos(beta_0)
        gamma_1 = 1/(cos_beta_1*Z1)

    M_1 = np.array([[np.cos(delta_1),                1j/gamma_1*np.sin(delta_1)   ],
               [1j*gamma_1*np.sin(delta_1),    np.cos(delta_1)             ]])
    M = M_1     # compute the transfer matrix of whole structure

    if polarisation == "TE":
        t_TE = 2*gamma_0 / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1])
        t.append(t_TE)
    elif polarisation == "TM":
        t_TM = 2*gamma_t / (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1]) * (Z_t / Z_0)
        t.append(t_TM)
    r_TE = 1* ((gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] - M[1,0] - gamma_t*M[1,1]) / 
            (gamma_0*M[0,0] + gamma_0*gamma_t*M[0,1] + M[1,0] + gamma_t*M[1,1]))
    r.append(r_TE)
    aux.append(delta_1/np.pi)

plt.figure(figsize=(12,10))
if plot_type == "amplitude":
    CS = plt.contourf(freqs/1e12, beta_0s, abs(np.array(r)), 100, cmap=cm.jet) # transmission amplitude
    for c in CS.collections: 
         c.set_antialiased(False)

    plt.colorbar() # draw colorbar
    plt.clabel(plt.contour(freqs/1e12, beta_0s, abs(np.array(r)), linewidths=0.5, levels=np.arange(0.,1.,.25), colors='k'))
elif plot_type == "phase":
    plt.contourf(freqs/1e12, beta_0s, angle(np.array(t)), 100, cmap=cm.hsv) # transmission phase
    plt.contour (freqs/1e12, beta_0s, abs(np.array(t)), [.998], colors='k') # transmission amplitude high
    plt.contour (freqs/1e12, beta_0s, abs(np.array(t)), [.002], colors='w') # transmission amplitude low

## For lossless dielectric, there are three sufficient conditions for 100% transmission:
## 1) Fabry-Perot resonance
#plt.clabel(plt.contour(beta_0, eps1s, np.array(aux), range(1,10), colors='black'))
## 2) If mu=1.0 and epsilon=1.0:
#plt.plot(beta_0, ones_like(beta_0), color='black', lw=2)
## 3) Brewster angle of incidence
#if polarisation == "TM":
    #plt.plot(np.arccos((1-eps1s)/(1+eps1s))/2, eps1s, color='black', lw=2)
    #plt.plot(np.arcsin(((1-1/eps1s)/(1-eps1s**2))**.5), eps1s, color='red', lw=1)



## 2) Zero impedance mismatch
#plt.plot(beta_0, ones_like(beta_0), color='blue')
    
plt.xlim(0, max(freqs)/1e12)
plt.ylim(0, np.pi/2)
plt.xlabel('Frequency [THz]')
plt.ylabel(u"Incidence angle $\\beta_0$ [rad]"); #plt.ylim((-0.1,1.1)); plt.xlim(xlim)
plt.grid()
#11plt.show()
plt.title('Reflection amplitude (analytic model)')
plt.savefig('BetaF_%s_slab_%s_XX.png' % (polarisation, plot_type), bbox_inches='tight')
