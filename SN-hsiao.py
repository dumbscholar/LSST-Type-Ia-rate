import numpy as np
import matplotlib.pyplot as plt
import cosmology as cc
import scipy.integrate as integ

on = -19
off = 86

data0=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day000.dat")
wl_data=data0[:,0]
flux_data=np.zeros((wl_data.size,(off-on)))
for f in range(on,off):
    data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day%03d.dat" %f)
    flux_data[:,f+19]=data[:,1]

#print flux_data[:,0]

wls = 1.0e-7 #wavelength step size in cm
scale = 1.145726382e55 #assuming peak B-band mag of -19.3
erg_joule=1.0e-7 #erg-to-Joule conversion factor

#print 'bolometric luminosity = ',wls*sum(flux*scale)*4*np.pi*1.0e7,' J/s' #from erg/s/cm/ster to J/s

a = cc.cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,np.log10(8.0),1.0) #change to 737 cosmology

h=0.7 #hubble parameter

z0=np.linspace(0.01,2.0,200) #avoiding zero redshift; creating redshift array

def SNR(z): #calculating supernova rate as function of redshift

    t=28.0/(1.+(1.+z)**2)

    f1 = lambda tD:tD**(-1.08)*(h*(0.0118+0.08*(np.sqrt(28.0/(t-tD)-1)-1))/(1+((np.sqrt(28.0/(t-tD)-1)-1)/3.3)**5.2))
    i1 = integ.quad(f1, 0.1, t)

    f2 = lambda tD:tD**(-1.08)
    i2 = integ.quad(f2, 0.1, 14.)

    SN=(0.04*0.032)*i1[0]/i2[0] #try to make it work with 1.28e-3 instead of 0.04*0.032
    
    return SN

def detect(day,z,wls,scale,success,lumdis,data_flux,data_wl,k,cap):
            wl=data_wl*(1+z) #cosmological length dilation
            flux=data_flux/(1+z) #cosmological time dilation
            wl_g=[]
            flux_g=[]
            m_g=30.0 #zero value for g-band apparent magnitude
            g_adjust=0.1 #used in absolute magnitude formula and based on scaling of the Hsiao spectral template

            for l in range(len(wl)):
                if wl[l]>=4000 and wl[l]<=5520: #LSST g-band
                   wl_g.append(wl[l])
                   flux_g.append(flux[l])

            if sum(flux_g)!=0:
               M_g=71.1974-g_adjust-2.5*np.log10(wls*sum(flux_g)*scale*4*np.pi*erg_joule) #from wikipedia
               m_g=M_g+5*np.log10(lumdis*100000)-5 #distance modulus relation
               
               if m_g<=26.2: #g-band limiting apparent magnitude
                  success+=1
                  cap.append(k)
               #print day,M_g,m_g,success,k
            return success, cap

nr = np.zeros([200]) #initializing array for number of supernova (occuring) as function of redshift

n=0

zstep=0.01 #related to redshift array definition

for i in range(200): #index for redshift
    nr[i]=0.0002327*(4/3)*np.pi*(a.Dcofz(z0[i]+zstep)**3-a.Dcofz(z0[i])**3)*SNR(z0[i])/(h**3) #cosmology code distance is (h Mpc^-3)
    print z0[i], nr[i]
    success=0
    lumdis=a.Dlofz(z0[i])/h #in Mpc
    cap=[] #array of captured SNe Ia

    flag=1
    fl=0
    for j in range(int(365*(1+z0[i]))): #index for day; generating observation dates
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

        

        if j<=365:
          for k in range(0,j):
            day=int((j-k)/(1+z0[i])) #look-back epoch OR apparent phase
            if k in cap:
               continue
            if day>=86:
               continue
            data_flux=flux_data[:,day]
            success, cap=detect(day,z0[i],wls,scale,success,lumdis,data_flux,wl_data,k,cap) #calling function detect()
            
        else:
          for k in range(0,365):
            day=int((j-k)/(1+z0[i])) #look-back epoch OR apparent phase
            if k in cap:
               continue
            if day>=86:
               continue
            data_flux=flux_data[:,day]
            success, cap=detect(day,z0[i],wls,scale,success,lumdis,data_flux,wl_data,k,cap) #calling function detect()

    nr[i]=(success/365.0)*nr[i] #array for number of supernova (detected in g-band) as function of redshift
    n+=nr[i]
    print success,nr[i]

   
plt.plot(z0,nr,color='k')
plt.xlabel('redshift, $z$',size=12)
plt.ylabel('number of SNe Ia detected in 1 year \n in g-filter by LSST in a given FOV',size=12)

print n

plt.show()
