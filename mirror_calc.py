#!/usr/bin/python
# coding: utf-8

from __future__ import division; from scipy.constants import *; from numpy import *

nS=1.46     ## refr. index of SiO2
nT=2.40     ## refr.index of TiO2
l0=550*nano  # wavelength in vacuum

print "a) Perpendicular incidence"
dS = pi/2/(nS*2*pi/l0*1)
dT = pi/2/(nT*2*pi/l0*1)
print "dS, dT", dS, dT
gammaS = nS
gammaT = nT
print "gammaS, gammaT", gammaS, gammaT
t = 2 / ((-gammaS/gammaT)**10 + (-gammaT/gammaS)**10)
print "10 bilayers -> Amplitude transmission: ", t, "(intensity t:", t**2, ")"


print
print
print "b) Nonperpendicular incidence"
cos_beta_SiO2 = sqrt(1-1/2/nS**2)
cos_beta_TiO2 = sqrt(1-1/2/nT**2)
dS = pi/2 / (nS * 2*pi/l0 * cos_beta_SiO2)
dT = pi/2 / (nT * 2*pi/l0 * cos_beta_TiO2)
print "dS, dT", dS, dT

print " b1) TE-polarisation"
gammaS = nS * cos_beta_SiO2
gammaT = nT * cos_beta_TiO2
print "gammaS, gammaT", gammaS, gammaT
t = 2 / ((-gammaS/gammaT)**10 + (-gammaT/gammaS)**10)
print "10 bilayers -> Amplitude transmission: ", t, "(intensity t:", t**2, ")"

print
print " b2) TM-polarisation (near Brewster angle)"
gammaS = nS / cos_beta_SiO2
gammaT = nT / cos_beta_TiO2
print "gammaS, gammaT", gammaS, gammaT
t = 2 / ((-gammaS/gammaT)**10 + (-gammaT/gammaS)**10)
print "10 bilayers -> Amplitude transmission: ", t, "(intensity t:", t**2, ")"
