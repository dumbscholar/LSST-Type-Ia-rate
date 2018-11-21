import matplotlib.pyplot as plt
import scipy.integrate as integ
import numpy as np
import cosmology as cc

a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0) # change to 737 cosmology

h=0.7

z0=np.linspace(0.,4.,401)

def SNR(z):

	t=28.0/(1.+(1.+z)**2)

	f1 = lambda tD:tD**(-1.08)*(h*(0.0118+0.08*(np.sqrt(28.0/(t-tD)-1)-1))/(1+((np.sqrt(28.0/(t-tD)-1)-1)/3.3)**5.2))
	i1 = integ.quad(f1, 0.1, t)

	f2 = lambda tD:tD**(-1.08)
	i2 = integ.quad(f2, 0.1, 14.)

	SN=(0.04*0.032)*i1[0]/i2[0] #try to make it work with 1.28e-3
	
	return SN


N = np.zeros([401])
lz = np.zeros([401])

for i in range(401):
    N[i]+=SNR(z0[i])
    lz[i]=0.8*N[i]/(1+z0[i])

plt.subplot(1,2,1)
#plt.semilogy(z0,N,color='k')
plt.plot(z0,N*10000,color='k')
plt.xlabel('redshift, $z$',size=12)
plt.ylabel('$n_{Ia}$ ($10^{-4}$ yr$^{-1}$Mpc$^{-3}$)',size=12)

plt.subplot(1,2,2)
plt.plot(z0,lz*10000,color='k')
plt.xlabel('redshift, $z$',size=12)
plt.ylabel('peak luminosity function, ($d\Phi_{Ia}$/dM)$^*$ ($10^{-4}$ yr$^{-1}$Mpc$^{-3}$)',size=12)

plt.show()

