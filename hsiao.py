import numpy as np
import matplotlib.pyplot as plt

wls = 1.0e-7 #step size
scale = 1.145726382e41 #assuming peak B-band mag of -19.3

day=np.linspace(-19,85,106)

M_u=np.zeros([len(day)])
M_b=np.zeros([len(day)])
M_v=np.zeros([len(day)])
M_r=np.zeros([len(day)])
M_i=np.zeros([len(day)])

Mmax_u=0
Mmax_b=0
Mmax_v=0
Mmax_r=0
Mmax_i=0

for i in range(-19,86):
    if i>=0 and i<10:
       data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day00"+str(i)+".dat") # give proper path
    if i>=10 and i<86:
       data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day0"+str(i)+".dat") # give proper path
    if i>=-19 and i<-9:
       data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day"+str(i)+".dat") # give proper path
    if i>=-9 and i<0:
       data=np.genfromtxt("/home/dumbscholar/Astro_codes/supernu/analysis/hsiao07/day-0"+str(-i)+".dat",) # give proper path

    wl=data[:,0]
    flux=data[:,1]

    wl_u=data[200:320,0]
    flux_u=data[200:320,1]
    wl_b=data[260:460,0]
    flux_b=data[260:460,1]
    wl_v=data[370:600,0]
    flux_v=data[370:600,1]
    wl_r=data[450:800,0]
    flux_r=data[450:800,1]
    wl_i=data[600:820,0]
    flux_i=data[600:820,1]

    M_u[i+20]=71.1974-0.87-2.5*np.log10(wls*sum(flux_u*scale)*4*np.pi*1.0e7)
    if M_u[i+20]<Mmax_u:
       Mmax_u=M_u[i+20]
       day_u=i+20

    M_b[i+20]=71.1974-2.5*np.log10(wls*sum(flux_b*scale)*4*np.pi*1.0e7)
    if M_b[i+20]<Mmax_b:
       Mmax_b=M_b[i+20]
       day_b=i+20

    M_v[i+20]=71.1974-0.45-2.5*np.log10(wls*sum(flux_v*scale)*4*np.pi*1.0e7)
    if M_v[i+20]<Mmax_v:
       Mmax_v=M_v[i+20]
       day_v=i+20

    M_r[i+20]=71.1974-0.71-2.5*np.log10(wls*sum(flux_r*scale)*4*np.pi*1.0e7)
    if M_r[i+20]<Mmax_r:
       Mmax_r=M_r[i+20]
       day_r=i+20

    M_i[i+20]=71.1974-1.455-2.5*np.log10(wls*sum(flux_i*scale)*4*np.pi*1.0e7)
    if M_i[i+20]<Mmax_i:
       Mmax_i=M_i[i+20]
       day_i=i+20

#print 'bolometric luminosity = ',wls*sum(flux*scale)*4*np.pi*1.0e7,' J/s' #from erg/s/cm/ster to J/s

print(Mmax_u,day_u)
print(Mmax_b,day_b)
print(Mmax_v,day_v)
print(Mmax_r,day_r)
print(Mmax_i,day_i)

plt.plot(day+20,M_u,'m-',linewidth=0.75)
plt.plot(day+20,M_b,'b-',linewidth=0.75)
plt.plot(day+20,M_v,'g-',linewidth=0.75)
plt.plot(day+20,M_r,'r-',linewidth=0.75)
plt.plot(day+20,M_i,'c-',linewidth=0.75)
plt.gca().invert_yaxis()
plt.ylim(-15.0,-20.0)
plt.xlabel('Day')
plt.ylabel('Absolute Magnitude')
plt.show()
