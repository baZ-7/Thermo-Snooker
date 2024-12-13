# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:02:13 2022

@author: aoc21
>>The Thermo Snooker Computing Project<<
File 4: Plots

This code file contains the following plots appended on the report:
- plot (1) Histogram of ball distance from container centre
- plot (2) Histogram of inter-ball separation
- plot (5) Pressure against temperature

"""
#%%
#imports
from tcp_Simulation import Simulation
import matplotlib.pyplot as pl

#%% 

#(1)Histogram of ball distance from container centre
c = Simulation(26, 100, 1, 5)
a1 = [] #arbitrary empty array
repeats1 = 100000 #running for this many frames
for i in range(0, repeats1): 
    c.run(1, animate = False)
    for i in range(len(c.dt_c)):
        a1.append(np.sqrt(c.ball[i]._position[0]**2 + c.ball[i]._position[1]**2))
    

#(1)plotting the histogram for ball positions from the center of the container
pl.hist(a1, bins=100, color='orange',edgecolor='black')
pl.title('Histogram of ball distance from container center')
pl.xlabel('Distance from Center')
pl.ylabel('Number of Balls')
pl.xlim(0,25)
pl.savefig('plot (1)')
pl.show()


#%%

#(2)Histogram of inter-ball separation
c = Simulation(26, 100, 1, 5)
a2 = [] #arbitrary empty array
repeats2 = 20000 #running for this many frames
for i in range(0, repeats2):
    c.run(1, animate = False)
    for i in range(len(c.dt_c)):
        for j in range(len(c.dt_c)):
            if i != j:
                a2.append(np.sqrt((c.ball[i]._position[0]-c.ball[j]._position[0])**2 + (c.ball[i]._position[1]-c.ball[j]._position[1])**2) - 2)
#to find to inter-ball separation, the -2 at the end to shift to graph to start at 0

#(2)plotting the histogram for inter-ball separations
pl.hist(a2, bins=100, color='orange', edgecolor='black')
pl.title('Histogram of inter-ball separations')
pl.xlabel('Separation')
pl.ylabel('Number of Balls')
pl.xlim(0,50)
pl.savefig('plot (2)')
pl.show()


#%%

#(5.1)Pressure against temperature changing volume
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(10, 15, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((10)*2*np.pi))
pl.plot(temp, pre, color='red', label='r=10')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 15, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='purple', label='r=20')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(30, 15, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((30)*2*np.pi))
pl.plot(temp, pre,color='blue', label='r=30')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(40, 15, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((40)*2*np.pi))
pl.plot(temp, pre, color='green', label='r=40')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(50, 15, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((50)*2*np.pi))
pl.plot(temp, pre, color='orange', label='r=50')
#(5)Plotting pressure against temperature
pl.title('Pressure against Temperature N=15')
pl.xlim(0,0.4e27)
pl.xlabel('Temperature (K)')
pl.ylabel('Pressure (Pa)')
pl.legend(framealpha=1, frameon=True, title='container radius')
pl.savefig('plot (5a)')
pl.show()


#%%

#(5.2)Pressure against temperature changing number of balls
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 10, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='purple', label='N=10')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 20, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='green', label='N=20')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 30, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='cyan', label='N=30')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 40, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='red', label='N=40')
pre = []
temp = []
for j in range(0, 1000):
    c = Simulation(20, 50, 1, j+1)
    c.run(10, animate = False)
    for k in range(len(c.dt_b)):
        ke = [] #ke
        ss = [] #square speed
        ke.append(c.ball[k].ke())
        ss.append(c.ball[k].velocity[0]**2+c.ball[k].velocity[1]**2)
    temp.append(np.sum(ke)/(1*c.num*1.38065*10e-23))
    pre.append(np.sum(ss)*1*0.5/((20)*2*np.pi))
pl.plot(temp, pre, color='orange', label='N=50')
#(5)Plotting pressure against temperature
pl.title('Pressure against Temperature r=20')
pl.ylim(0,5000)
pl.xlim(0,0.8e27)
pl.xlabel('Temperature (K)')
pl.ylabel('Pressure (Pa)')
pl.legend(framealpha=1, frameon=True, title='number of balls')
pl.savefig('plot (5b)')
pl.show()

#%%

