from turtle import fd
import numpy as np
import math
I=np.array([1,0,0])
J=np.array([0,1,0])
K=np.array([0,0,1])
mu=398600
k=1
while k==1:
    choice=int(input('RV to OE (1) or OE to RV (2): '))
    if choice==1:
       k=0
       rvec1=float(input('Input I component of position vector in km: '))
       rvec2=float(input('Input J component of position vector in km: '))
       rvec3=float(input('Input K component of position vector in km: '))
       vvec1=float(input('Input I component of velocity vector in km/s: '))
       vvec2=float(input('Input J component of velocity vector in km/s: '))
       vvec3=float(input('Input K component of velocity vector in km/s: '))
       rvec=np.array([rvec1,rvec2,rvec3])
       vvec=np.array([vvec1,vvec2,vvec3])
       #rvec=np.array([-2436.45, -2436.45, 6891.0379])
       #vvec=np.array([5.088611, -5.088611, 0])
       r=((rvec[0]**2)+(rvec[1]**2)+(rvec[2]**2))**(1/2)
       v=((vvec[0]**2)+(vvec[1]**2)+(vvec[2]**2))**(1/2)
       hvec=np.cross(rvec,vvec)
       h=((hvec[0]**2)+(hvec[1]**2)+(hvec[2]**2))**(1/2)
       energy=.5*v**2-mu/r
       a=-mu/(2*energy)
       p=h**2/mu
       hhat=hvec/h
       i=math.acos(np.dot(hhat,K))
       rhat=rvec/r
       evec=np.cross(vvec,hvec)/mu-rhat
       e=((evec[0]**2)+(evec[1]**2)+(evec[2]**2))**(1/2)
       ehat=evec/e
       if np.dot(rvec,vvec)>=0:
          f=math.acos(np.dot(ehat,rhat)-.0000001)
       elif np.dot(rvec,vvec)<0:
          f=2*math.pi-math.acos(np.dot(ehat,rhat)-.0000001)
       ph=np.cross(K,hvec)
       phn=((ph[0]**2)+(ph[1]**2)+(ph[2]**2))**(1/2)
       nhat=np.cross(K,hvec)/phn
       if np.dot(ehat,K)>=0:
          w=math.acos(np.dot(nhat,ehat)-.0000001)
       elif np.dot(ehat,K)<0:
          w=2*math.pi-math.acos(np.dot(nhat,ehat)-.0000001)
       omega=math.atan2(nhat[1],nhat[0])
       print('a=',a)
       print('f=',f)
       print('e=',e)
       print('i=',i)
       print('w=',w)
       print('omega=',omega)
    elif choice==2:
       k=0
       a=float(input('Input semi major axis (a) in km: '))
       f=float(input('Input true anomoly (f) in radians: '))
       e=float(input('Input eccentricity (e): '))
       i=float(input('Input inclination (i) in radians: '))
       w=float(input('Input argument of periapsis (w) in radians: '))
       omega=float(input('Input RA of ascending nodes (capital omega) in radians: '))
       #a=7712.2
       #f=1.4901*10**(-8)
       #e=.001
       #i=1.1071
       #w=1.5708
       #omega=2.3562
       p=a*(1-e**2)
       if p==0:
          p=a
       h=(p*mu)**(1/2)
       r=p/(1+e*math.cos(f))
       rpqw=np.array([r*math.cos(f),r*math.sin(f),0])
       vpqw=(mu/h)*np.array([-math.sin(f),e+math.cos(f),0])
       v=((vpqw[0]**2)+(vpqw[1]**2)+(vpqw[2]**2))**(1/2)
       rot1=np.array([[math.cos(omega),-math.sin(omega),0],[math.sin(omega),math.cos(omega),0],[0,0,1]])
       rot2=np.array([[1,0,0],[0,math.cos(i),-math.sin(i)],[0,math.sin(i),math.cos(i)]])
       rot3=np.array([[math.cos(w),-math.sin(w),0],[math.sin(w),math.cos(w),0],[0,0,1]])
       rot4=np.matmul(rot1,rot2)
       rot=np.matmul(rot4,rot3)
       vvec=np.matmul(rot,vpqw)
       rvec=np.matmul(rot,rpqw)
       print('Velocity in ECI:',vvec)
       print('Position Vector in ECI:',rvec)
    else:
       print('Invalid input\n')
       


   
  