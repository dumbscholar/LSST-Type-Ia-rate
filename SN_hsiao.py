import numpy as np
import matplotlib.pyplot as plt
import cosmology as cc
import scipy.integrate as integ

wls = 1.0e-7 #step size
scale = 1.145726382e41 #assuming peak B-band mag of -19.3

#print 'bolometric luminosity = ',wls*sum(flux*scale)*4*np.pi*1.0e7,' J/s' #from erg/s/cm/ster to J/s

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

for i in range(150): #index for redshift
    nr[i]=0.0002327*(4/3)*np.pi*(a.Dcofz(z0[i]+0.01)**3-a.Dcofz(z0[i])**3)*SNR(z0[i])
    print z0[i], nr[i]
    success=0

    flag=1
    fl=0
    for j in range(int(365*(1+z0[i]))): #index for day
        if flag==0:
           if (j-fl)%17!=0: #17th day; 17.5 days cadence
              continue
           else:
              flag=1
              fl+=1
        elif flag==1:
           if j%(17.5*2)!=0: #18th day; 17.5 days cadence
              continue
           else:
              flag=0
        
        for k in range(success,j):
            day=int(k/(1+z0[i]))-19 #look-back epoch OR apparent phase
            if day>=0 and day<10:
               data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day00"+str(day)+".dat")
            elif day>=10 and day<86:
               data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day0"+str(day)+".dat")
            elif day>=-19 and day<-9:
               data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day"+str(day)+".dat")
            elif day>=-9 and day<0:
               data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day-0"+str(-day)+".dat")
            else:
               continue

            wl=data[:,0]*(1+z0[i]) #length dilation
            flux=data[:,1]/(1+z0[i]) #time dilation
            wl_g=[]
            flux_g=[]
            m_g=30
            for l in range(len(wl)):
                if wl[l]>=4000 and wl[l]<=5520: #LSST g-band
                   wl_g.append(wl[l])
                   flux_g.append(flux[l])

            if sum(flux_g)!=0:
               M_g=71.1974-0.1-2.5*np.log10(wls*sum(flux_g)*scale*4*np.pi*1.0e7)
               m_g=M_g+5*np.log10(a.Dlofz(z0[i])*100000)-5
               
               if m_g<=26.2: #g-band limiting magnitude
                  success+=1
               #print M_g,m_g,success
    nr[i]=success/(365*(1+z0[i]))*nr[i]
    n+=nr[i]
    print success,nr[i]

   
plt.plot(z0,nr,color='k')
plt.xlabel('redshift, $z$',size=12)
plt.ylabel('number of SNe Ia detected per year by LSST',size=12)

print n

plt.show()
