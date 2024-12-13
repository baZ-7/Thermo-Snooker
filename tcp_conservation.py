# -*- coding: utf-8 -*-
"""
Created on Sat Dec  3 13:30:02 2022

@author: aoc21
>>The Thermo Snooker Computing Project<<
File 3: Conservation Tests

This code file contains the tests for conservation:
- testing for the conservation of kinetic energy and momentum through time
- also contains the plot (3) & (4) which corresponds to KE and momentum against time

"""
#%%
#imports
from tcp_Simulation import Simulation

#%%

#(3,4) Conservation Tests for (3) System Kinetic Energy & (4) Momentum
c = Simulation(20, 50, 1, 10)
syske = [] #array of system KE
mo = [] #array of mo
t = [] #array of time
dta = 0 #dt average
arb = 0 #arbitrary constant
repeats3 = 100000 #running for this many frames
for i in range(0, repeats3):
    c.run(1, animate=False)
    kel = [] #empty KE list
    mel = [] #empty momentum list
    for k in range(0, len(c.dt_c)):
        kel.append(c.ball[k].ke())
        mel.append((c.cim))
        
    syske.append(np.sum(kel)) #system ke
    mo.append(np.sum(mel)) #average change in momentum

    tmin = 0 #arbitrary time constant
    if c.b_min <= c.c_min:
        tmin = c.b_min
    else:
        tmin = c.c_min 
    if arb == 0:
        t.append(tmin) #appending first value
        arb = arb + 1
    else:
        t.append(tmin + t[arb-1]) #appending time values
        arb = arb + 1 #for the accumulation of dt
        
#%%        

#(3)plotting system KE against time
a3,b3 = np.polyfit(t, syske, 1)
fitke = np.multiply(t,a3) + b3
pl.plot(t,syske,'x', label='raw data',color='lightgray')
pl.plot(t,fitke, label='best fit line',color='red')
pl.yticks([])
pl.title('System Kinetic Energy against Time')
pl.xlabel('Time (t)')
pl.ylabel('System Kinetic Energy (J)')
pl.legend(framealpha=1, frameon=True)
pl.savefig('plot (3)')
pl.show()


#%%

#(4)plotting total momentum against time
a4,b4 = np.polyfit(t, mo, 1)
pl.plot(t,mo,'x', label='raw data',color='lightgray')
fitmo = np.multiply(t,a4) + b4
pl.plot(t,fitmo, label='best fit line',color='red')
pl.title('Momentum against Time')
pl.xlabel('Time (t)')
pl.ylabel('Momentum (kgms^-1)')
pl.legend(framealpha=1, frameon=True)
pl.savefig('plot (4)')
pl.show()
