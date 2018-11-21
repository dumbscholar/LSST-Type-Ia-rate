import numpy as np
import matplotlib.pyplot as plt

z=np.linspace(0.0, 6.0, 61)
rate=(0.0118+0.08*z)*0.7/(1+(z/3.3)**5.2)
lz=np.linspace(0.0, 0.85, 18)
ratel=(0.0118+0.08*(10.**lz-1.))*0.7/(1+((10.**lz-1.)/3.3)**5.2)

plt.subplot(1,2,1)
plt.semilogy(z, rate,'k-')
plt.xlabel(r'redshift, $z$',size=12)
plt.ylabel(r'Star formation rate, $\rho_{SFR}$ (M$_{\odot}$yr$^{-1}$Mpc$^{-3}$)',size=12)

plt.subplot(1,2,2)
plt.semilogy(lz, ratel,'k-')
plt.xlabel(r'log $(1+z)$',size=12)
plt.ylabel(r'Star formation rate, $\rho_{SFR}$ (M$_{\odot}$yr$^{-1}$Mpc$^{-3}$)',size=12)

plt.subplots_adjust(wspace=0.3,hspace=None)
plt.show()
