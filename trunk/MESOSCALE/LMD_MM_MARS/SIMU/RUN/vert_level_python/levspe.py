#! /usr/bin/env python


import numpy as np
from ppclass import pp
import ppplot
from math import *
import matplotlib.pyplot as plt

param=np.loadtxt('paramlevspe')
hache=param[0]
psurf=param[1]
nlev=param[2]
altmax=param[3]
epsilon=param[4]
elong_cos=param[5]

print nlev, altmax, hache, psurf, epsilon, elong_cos

epsilon = epsilon / (nlev-1) / 100.
#print, 'epsilon',epsilon
#
exposant = np.pi/nlev/elong_cos
#
# alpha (plus grand ecart en km) est determine pour que max(altitudes)=altmax
#
alpha =  altmax / ( (nlev-1)/2. - np.sin(2.*exposant*(nlev-1))/4./exposant + epsilon*(nlev-1)**2/2. )
#alpha = 1.02799
print alpha
#
x=np.linspace(0,nlev-1,nlev)
#exit()
#print,'x',x

#
# CALC
#
#altitudes = alpha*x/2. - alpha*sin(2.*exposant*x)/4./exposant + epsilon*alpha*x**2/2.
altitudes = x*alpha/2. - alpha*np.sin(2.*exposant*x)/4./exposant + epsilon*alpha*x**2/2.
#print altitudes

#print altitudes
logpressions = np.log(psurf) - altitudes/hache
pressions=np.exp(logpressions)
ptop=psurf*np.exp(-altmax/hache)
#print ptop
etas=(pressions-ptop)/(psurf-ptop)
etas[nlev-1]=0.
#print etas
press=etas*(psurf-ptop)+ptop
#print press
pseudo=10.*np.log(psurf/press)
#print pseudo
#diff=pseudo - shift(pseudo,-1) & diff=-diff(0:nlev-2)
res=[]
res.append(0)
for i in range(int(nlev))[1:] :
	res.append(pseudo[i]-pseudo[i-1])
#print res

#save etas
np.savetxt('levels',etas)

plt.figure(figsize=(15, 15))
plt.subplot(221)
plt.plot(x, etas)
plt.xlabel('levels')
plt.ylabel('etas')
plt.grid()
#plt.title('a) NINO3 Sea Surface Temperature (seasonal)')
#plt.hold(False)

plt.subplot(222)
plt.plot(x, pseudo)
plt.xlabel('levels')
plt.grid()
plt.ylabel('pseudo-altitude (km)')

plt3 = plt.subplot(223)
plt.semilogy(x, press)
plt.xlabel('levels')
plt.grid()
plt.ylabel('pression (Pa)')

plt4 = plt.subplot(224)
plt.plot(x, res)
plt.xlabel('levels')
plt.grid()
plt.ylabel('resolution (km)')

plt.show()










