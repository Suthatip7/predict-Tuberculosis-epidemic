# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:40:22 2019

@author: User
"""

import scipy.integrate
import numpy
import matplotlib.pyplot as plt

def SIR_model(y,t,a,b,d,e,f,g,h,i,u0,u1,u2,u,r,j,k):
    S,EN,ER,IN,IR,R,I =y
    dS_dt=a+k*i*R-(u0+(1-r)*b*IN+(1-u)*b*IR)*S
    dEN_dt=(1-r)*b*S*IN-(u0+d)*EN
    dER_dt=(1-u)*b*S*IR-(u0+e)*ER
    dIN_dt=d*EN-(u0+u1+g+j*f)*IN
    dIR_dt=e*ER+j*f*IN-(u0+u2+h)*IR
    dR_dt=g*IN+h*IR-(u0+u*i)*R
    dI_dt=dIN_dt+dIR_dt


    return([dS_dt,dEN_dt,dER_dt,dIN_dt,dIR_dt,dR_dt,dI_dt])
    
#N0=5000
S0=5000.0
EN0=0.0
ER0=0.0
IN0=488.0
IR0=1.0
I0=IN0+IR0
#R0=N0-(S0+EN0+ER0+I0)
R0=1000.0
a=(S0+EN0+ER0+I0+R0)/(75.55*365.0)
r=0.0
u=0.5
b=1/(230000.0)
#r0=3.0
d=1/(50.0) #14-730 days
e=1/(20.0) #14-730 days
j=0.00025
f=1/(0.5*31.0) #<6month
g=1/(6.0*31.0) #ประมาณ 6month
h=1/(9*31.0) #9-24month 
k=0.05
i=1/(5.0*365.0)
u0=1/(76.55*365.0)
u1=1/(5.0*365.0)
u2=1/(3.0*365.0)
tb=[407,670,558,639,759,760,1484,1555,1397,1582,1397,1620,1236,1621,1425,1436,1618,1666,939,888,823,821,748,642]
mdr=[1,1.5,2,1,1,3,4,3,6,6,4,7,9,7,8,5,3,6,6,5,7,3,4,2]
tb_mdr=[408,671.5,560,640,760,763,1488,1558,1403,1588,1401,1627,1245,1628,1433,1441,1621,1672,945,893,830,824,752,644]
#b=(r0*u0*(u0+e)*(u0+u2+h))/((1-u)*a*e)
r2=((1-u)*a*b*e)/(u0*(u0+e)*(u0+u2+h))
S2u=-(u*a*b*e)/(r2*u0*(u0+e)*(u0+u2+h))
S2a=((1-u)*a*b*e)/(r2*u0*(u0+e)*(u0+u2+h))
S2b=((1-u)*a*b*e)/(r2*u0*(u0+e)*(u0+u2+h))
S2e=((1-u)*a*b*e)/(r2*((u0+e)**2)*(u0+u2+h))
S2u0=-((1-u)*a*b*e*(3*(u0**2)+2*u0*(h+e+u2)+e*(h+u2)))/(r2*u0*((u0+e)**2)*((u0+u2+h)**2))
S2u2=-(u2*(1-u)*a*b*e)/(r2*u0*(u0+e)*((u0+u2+h)**2))
S2h=-(h*(1-u)*a*b*e)/(r2*u0*(u0+e)*((u0+u2+h)**2))
r1=((1-r)*a*b*d)/(u0*(u0+d)*(u0+u1+g+j*f))
S1r=-(r*a*b*d)/(r1*u0*(u0+d)*(u0+u1+g+j*f))
S1a=((1-r)*a*b*d)/(r1*u0*(u0+d)*(u0+u1+g+j*f))
S1b=((1-r)*a*b*d)/(r1*u0*(u0+d)*(u0+u1+g+j*f))
S1d=((1-r)*a*b*d)/(r1*((u0+d)**2)*(u0+u1+g+j*f))
S1u0=-((1-r)*a*b*d*(3*(u0**2)+2*u0*(j*f+d+u1)+d*(j*f+g+u1)))/(r1*u0*((u0+d)**2)*((u0+u1+g+j*f)**2))
S1u1 =-(u1*(1-r)*a*b*d)/(r1*u0*(u0+d)*((u0+u1+g+j*f)**2))
S1g=-(g*(1-r)*a*b*d)/(r1*u0*(u0+d)*((u0+u1+g+j*f)**2))
S1j=-(j*f*(1-r)*a*b*d)/(r1*u0*(u0+d)*((u0+u1+g+j*f)**2))
S1f=-(j*f*(1-r)*a*b*d)/(r1*u0*(u0+d)*((u0+u1+g+k*j)**2))

t=numpy.linspace(0,5000,100000)
day=[0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690]

solution=scipy.integrate.odeint(SIR_model,[S0,EN0,ER0,IN0,IR0,R0,I0],t,args=(a,b,d,e,f,g,h,i,u0,u1,u2,u,r,j,k))
solution=numpy.array(solution)



#plt.figure(figsize=[6,4])
#plt.figure(1)
#plt.plot(t,solution[:,0],label="S(t)",color="magenta")
#plt.plot(t,solution[:,1],label="EN(t)",color="orange")
#plt.plot(t,solution[:,2],label="ER(t)",color="green")
#plt.plot(t,solution[:,3],label="IN(t)",color="red")
#plt.plot(t,solution[:,4],label="IR(t)",color="purple")
#plt.plot(t,solution[:,5],label="R(t)",color="brown")
#plt.plot(t,solution[:,6],label="I(t)",linestyle='--')
#plt.grid()
#plt.legend()
#plt.xlabel("Time")
#plt.ylabel("Proportions")
#plt.title("SEIR Model")
#plt.figure(2)
#plt.plot(t,solution[:,0],label="S(t)",color="magenta")
#plt.grid()
#plt.legend()
#plt.figure(3)
#plt.plot(t,solution[:,1],label="EN(t)",color="orange")
#plt.plot(t,solution[:,2],label="ER(t)",color="green")
#plt.grid()
#plt.legend()
plt.figure(4)
plt.scatter(day,tb,color="red")
plt.plot(day,tb,label="case",color="red")
plt.plot(t,solution[:,3],label="IN(t)",color="red",linestyle='--')
plt.grid()
plt.legend()
plt.title("the number of  normal infected individuals")
plt.xlabel("Time(day)")
plt.ylabel("Proportions")
#plt.figure(5)
#plt.plot(t,solution[:,5],label="R(t)",color="brown")
#plt.grid()
#plt.legend()
plt.figure(6)
plt.scatter(day,mdr,color="blue")
plt.plot(day,mdr,label="case",color="blue")
plt.plot(t,solution[:,4],label="IR(t)",color="blue",linestyle='--')
plt.grid()
plt.legend()
plt.title("the number of drug resistance infected individuals")
plt.xlabel("Time(day)")
plt.ylabel("Proportions")
plt.figure(7)
plt.scatter(day,tb,color="purple")
plt.plot(day,tb_mdr,label="case",color="purple")
plt.plot(t,solution[:,6],label="I(t)",linestyle='--',color="purple")
plt.grid()
plt.legend()
plt.title("the number of infected individuals")
plt.xlabel("Time(day)")
plt.ylabel("Proportions")
#plt.figure(8)
#plt.plot(t,solution[:,3],label="IN(t)",color="red",linestyle='--')
#plt.plot(t,solution[:,4],label="IR(t)",color="purple",linestyle='--')
#plt.grid()
#plt.legend()


plt.show()



print("reprodective 2 =",r2)
print("reprodective 1 =",r1)
print("Sensitivity analysis2 u=",S2u)
print("Sensitivity analysis2 a=",S2a)
print("Sensitivity analysis2 b=",S2b)
print("Sensitivity analysis2 e=",S2e)
print("Sensitivity analysis2 u0=",S2u0)
print("Sensitivity analysis2 u2=",S2u2)
print("Sensitivity analysis2 h=",S2h)
print("Sensitivity analysis1 r=",S1r)
print("Sensitivity analysis1 a=",S1a)
print("Sensitivity analysis1 b=",S1b)
print("Sensitivity analysis1 d=",S1d)
print("Sensitivity analysis1 u0=",S1u0)
print("Sensitivity analysis1 u1=",S1u1)
print("Sensitivity analysis1 g=",S1g)
print("Sensitivity analysis1 j=",S1j)
print("Sensitivity analysis1 f=",S1f)
