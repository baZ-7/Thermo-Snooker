# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 14:19:17 2022

@author: aoc21
File 5: Curve-fit plots

This code file contains the plots involving curvefitting
- plot (6) How changing ball radius affects equation of state
- plot (7) Histogram of ball speeds fitted with theoretical (Maxwell Boltzmann) 
- plot (8) Determining coefficients of Van Der Waals law analytically and graphically
"""
#%%
#imports
from tcp_Simulation import Simulation
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

#%%

#(6)How ball radius affect the equation of state
#pV graph
def inverse(x, i, c):
    inv = i*(1/x) + c
    return inv
pre1 = []
vol1 = []
for j in range(0, 10000):
    c = Simulation(10+j, 15, 1e-15, 1500)
    c.run(5, animate = False)
    for k in range(len(c.dt_c)):
        ss = [] #square speed
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    mss = np.sum(ss)
    pre1.append(mss*1*0.5/((10+j)*2*np.pi))
    vol1.append(((10+j)*2)*np.pi)
i = 1000 #guess for i
c = 0 #guess for c
guess3 = [i, c]
yt3, ycov3 = curve_fit(inverse,vol1,pre1,p0=guess3)
print('for r = 1e-15, i =', yt3[0], 'and c =', yt3[1])
t_pre1 = []
for i in range(len(vol1)):
    t_pre1.append(inverse(vol1[i],yt3[0],yt3[1]))
pl.plot(vol1, pre1, 'x', color='lightblue', label='raw data r=1e-15')
pre2 = []
vol2 = []
for j in range(0, 10000):
    c = Simulation(10+j, 15, 1, 1500)
    c.run(5, animate = False)
    for k in range(len(c.dt_c)):
        ss = [] #square speed
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    mss = np.sum(ss)
    pre2.append(mss*1*0.5/((10+j)*2*np.pi))
    vol2.append(((10+j)*2)*np.pi)
guess4 = [i, c]
yt4, ycov4 = curve_fit(inverse,vol2,pre2,p0=guess3)
print('for r = 1, i =', yt4[0], 'and c =', yt4[1])
t_pre2 = []
for i in range(len(vol2)):
    t_pre2.append(inverse(vol2[i],yt4[0],yt4[1]))
#(6)Plotting pressure against temperature
pl.plot(vol2, pre2, 'x', color='lightgray', label='raw data r=1')
pl.plot(vol1, t_pre1, color='blue', label='fit line r=1e-15')
pl.plot(vol2, t_pre2, color='black', label='fit line r=1')
pl.title('Pressure against Volume N=15')
pl.xlabel('Volume (m^2)')
pl.ylabel('Pressure (Pa)')
pl.legend(framealpha=1, frameon=True, title='ball radius')
pl.xlim(0, 3000)
pl.ylim(0, 6000)
pl.savefig('plot (6)')
pl.show()


#%%

#(7)plotting histogram of ball velocities
#defining the function
def mb (v, T, a):
    mb = v*np.exp(-(0.5*1*v**2)/((1.380649*10e-23)*T))*a
    return mb

#running simulation
c = Simulation(25, 100, 1, 5)
a7 = [] #arbitrary empty array
repeats1 = 10000
0 #running for this many frames
for i in range(0, repeats1): 
    c.run(1, animate = False)
    for i in range(len(c.dt_c)):
        a7.append(np.sqrt(c.ball[i].velocity[0]**2 + c.ball[i].velocity[1]**2))
        
bins = 100
h, xedges, patches = pl.hist(a7, bins=bins, color='dimgray', label='Velocity Distribution')
pl.xlabel('Ball Speeds')
pl.ylabel('Number of Balls')
bincenters = (xedges[:-1] + xedges[1:]) / 2

#plotting theoretical curve
a0 = 150
temp_guess = 6e22
guess = [temp_guess, a0] #initial guess
yt, ycov = curve_fit(mb, bincenters, h, p0=guess) #curve-fitting
pl.plot(bincenters, mb(bincenters, yt[0], yt[1]), color='red', label='Maxwell Boltzmann')
pl.legend(framealpha=1, frameon=True)
pl.savefig('plot (7)')
pl.show()

#%%

#(8)determining coefficients of Van Der Waals law
#defining the function
def vdw (p, V, N, a, b):
    T = (p*V + a*(N/V)**2 - p*N*b - (a*(N**3)*b)/(V**2))/(N*(1.380649*10e-23))
    return T


#analytically assuming ideal gas:
#    b = V/N
#and then a can be obtained graphically.
#using a 15 ball situation as an example:

c = Simulation(10, 15, 1, 5)
a = 0 #initial guess
b = 15*np.pi*(10)**2 #finding b deducing analytically
pre = []
temp = []
for j in range(0, 100):
    c = Simulation(10, 15, 1, j+0.1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    totalke = np.sum(ke) #total ke
    mss = np.sum(ss)
    temperature = totalke/(1*c.num*1.38065*10e-23) #temperature
    temp.append(temperature)
    pre.append(mss*1*0.5/(10*2*np.pi))
guess2 = [(10*2)*np.pi, 15, a, b]
yt2, ycov2 = curve_fit(vdw, pre, temp, p0=guess2) #curve-fitting
print('a =',yt2[2],'+-',ycov2[2][2],'and b =',yt2[3],'+-',ycov2[3][3])
t_temp = []
for i in range(len(pre)):
    t_temp.append(vdw(pre[i], (10*2)*np.pi, 15, yt2[2], yt2[3]))
a8,b8 = np.polyfit(t_temp, pre, 1)
fitpre = np.multiply(t_temp,a8) + b8
pl.plot(t_temp, pre, 'x', label = 'fitted data')
pl.plot(t_temp, fitpre, label = 'best fit line')
pl.title('Van Der Waals')
pl.xlabel('Temperature (K)')
pl.ylabel('Pressure (Pa)')
pl.legend(framealpha=1, frameon=True)
pl.savefig('plot (8)')
pl.show()








    