import matplotlib.pyplot as plt
import scipy.integrate as integ
import numpy as np
import cosmology as cc

data=np.genfromtxt('lc.dat') #normal SN Ia with M_B=19.3; Nugent template for normal SN Ia

epoch=data[:,0]
Mu=data[:,1]
Mb=data[:,2]
Mv=data[:,3]
Mr=data[:,4]
Mi=data[:,5]

Lu=(3.828*10.**26)*10.**(0.4*(4.74-Mu)) #u-band luminosity
Lb=(3.828*10.**26)*10.**(0.4*(4.74-Mb)) #b-band luminosity
Lv=(3.828*10.**26)*10.**(0.4*(4.74-Mv)) #v-band luminosity
Lr=(3.828*10.**26)*10.**(0.4*(4.74-Mr)) #r-band luminosity
Li=(3.828*10.**26)*10.**(0.4*(4.74-Mi)) #i-band luminosity      

a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0) #change to 737 cosmology

h=0.7

z0=np.linspace(0.01,1.5,150) #avoiding zero redshift

def SNR(z):

    t=28.0/(1.+(1.+z)**2)

    f1 = lambda tD:tD**(-1.08)*(h*(0.0118+0.08*(np.sqrt(28.0/(t-tD)-1)-1))/(1+((np.sqrt(28.0/(t-tD)-1)-1)/3.3)**5.2))
    i1 = integ.quad(f1, 0.1, t)

    f2 = lambda tD:tD**(-1.08)
    i2 = integ.quad(f2, 0.1, 14.)

    SN=(0.04*0.032)*i1[0]/i2[0] #try to make it work with 1.28e-3
    
    return SN

nr = np.zeros([150])

n=0

for i in range(150):
    print z0[i]
    nr[i]=0.43633*(4/3)*np.pi*(a.Dcofz(z0[i]+0.01)**3-a.Dcofz(z0[i])**3)*SNR(z0[i])
    discard=0
    for j in range(89*5): #a fixed sample size (89*5 is the size of data) of Type Ia supernovae at some redshift
        k=np.random.choice(range(89)) #there are 89 elements in epoch array
        dlum=3.085678*10**22*a.Dlofz(z0[i]) #luminosity distance in metres
        Fu=Lu[k]/(4*np.pi*dlum**2) #100 = 10^2
        Fb=Lb[k]/(4*np.pi*dlum**2)
        Fv=Lv[k]/(4*np.pi*dlum**2)
        Fr=Lr[k]/(4*np.pi*dlum**2)
        Fi=Li[k]/(4*np.pi*dlum**2)
        sFu=Fu/((1+z0[i])*(420-300)) #filter ranges taken from Bessell90; take half-max ranges?
        sFb=Fb/((1+z0[i])*(560-360)) #specific flux or flux density in terms of per nm
        sFv=Fv/((1+z0[i])*(700-470)) #equally distributing the flux across the band for simplicity
        sFr=Fr/((1+z0[i])*(900-550))
        sFi=Fi/((1+z0[i])*(920-700))
        tfu=tfg=tfr=tfi=tfz=tfy=0.0        
        band=np.linspace(300*(1+z0[i]),920*(1+z0[i]),1000)
        flux=np.zeros([len(band)])
        for p in range(len(band)):
            if band[p]>=300*(1+z0[i]) and band[p]<=420*(1+z0[i]):
                flux[p]+=sFu*0.620*(1+z0[i])
            if band[p]>=360*(1+z0[i]) and band[p]<=560*(1+z0[i]):
                flux[p]+=sFb*0.620*(1+z0[i])
            if band[p]>=470*(1+z0[i]) and band[p]<=700*(1+z0[i]):
                flux[p]+=sFv*0.620*(1+z0[i])
            if band[p]>=550*(1+z0[i]) and band[p]<=900*(1+z0[i]):
                flux[p]+=sFr*0.620*(1+z0[i])
            if band[p]>=700*(1+z0[i]) and band[p]<=920*(1+z0[i]):
                flux[p]+=sFi*0.620*(1+z0[i])      
        for q in range(len(band)):
            if band[q]>=350 and band[q]<400: 
                tfu+=flux[q]
            if band[q]>=400 and band[q]<552:
                tfg+=flux[q]
            if band[q]>=552 and band[q]<691:
                tfr+=flux[q]
            if band[q]>=691 and band[q]<818:
                tfi+=flux[q]
            if band[q]>=818 and band[q]<922:
                tfz+=flux[q]
            if band[q]>=922 and band[q]<1060:
                tfy+=flux[q]
        mu=mg=mr=mi=mz=my=30.0 #arbitrary upper limit 
        if tfu!=0:
            mu=-26.74-2.5*np.log10(tfu/1361)   
        if tfg!=0:
            mg=-26.74-2.5*np.log10(tfg/1361)
        if tfr!=0:
            mr=-26.74-2.5*np.log10(tfr/1361)
        if tfi!=0:
            mi=-26.74-2.5*np.log10(tfi/1361)
        if tfz!=0:
            mz=-26.74-2.5*np.log10(tfz/1361)
        if tfy!=0:
            my=-26.74-2.5*np.log10(tfy/1361)
        #print mu, mg, mr, mi, mz, my
        if mu>24.9 and mg>26.2 and mr>26.4 and mi>25.7 and mz>25.0 and my>23.7:
            discard+=1
    print nr[i], discard*nr[i]/(89*5) #number of type Ia occuring in the LSST sky region at redshift & number detected
    nr[i]-=discard*nr[i]/(89*5)
    n+=nr[i]

plt.plot(z0,nr,color='k')
plt.xlabel('redshift, $z$',size=12)
plt.ylabel('number of SNe Ia detected per year by LSST',size=12)

print n

plt.show()
