# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 13:41:18 2022

@author: aoc21
>>The Thermo Snooker Computing Project<<
File 2: Simulation Class

This code file contains the Simulation class:
For creating a Simulation of the thermodynamic snooker situation 
- using the Ball class from file 1
Running this file will output an animation characterised by the parameters at the bottom of the file.
(if the simulation section in 'un-docstringed')
"""

#%%
#imports
from tcp_Ball import Ball 
import numpy as np
import pylab as pl 
import random

#%%
#Class Simulation
class Simulation:
    def __init__(self, Ra, n, ra, v):
        '''
        initiating the Simulation
         ===== variables =====
        Ra --- radius of the container
        n --- number of balls
        ra --- radius of the ball
        v --- velocity of ball i.e. from -v to v
        ======================
        '''
        #generating the container
        self.container = Ball(R=Ra, con=True)
        self.num = n
        #setting a array of ball objects
        self.ball = []
        #radius of ball
        self.ra = ra
        #radius of container
        self.Ra = Ra
        #velocity
        self.v = v
        #change in momentum
        self.cim = 0
        #initial momentum
        self.inimo = 0
        #generating a random assortment of ball objects
        #using polar coordinates r and t for theta 
        r = []
        t = []
        vx = []
        vy = []
        #generating positions for the balls in a systematic ring formation
        ballcount = self.num #counting number of balls
        ringcount = 2 #counting number of rings
        while ballcount > 0:
            #using triangular numbers
            a = (ringcount*(ringcount+1)*0.5) #arbitrary constant
            for i in range(0,int(a)):
                r.append((ringcount-1)*(self.ra*2+1))
                t.append((i+1)*((2*np.pi)/a))
            ballcount = ballcount - a
            ringcount = ringcount + 1
            b = (ringcount*(ringcount+1)*0.5) #arbitrary constant
            if ballcount < b:
                for i in range(0,int(b)):
                    r.append((ringcount-1)*(self.ra*2+1))
                    t.append((i+1)*((2*np.pi)/b))
            else:
                pass
        
        for i in range(self.num):
            vx.append(random.uniform(-v, v))
            vy.append(random.uniform(-v, v))
          
        #computing the random data into the ball objects
        for i in range(self.num):
            a = Ball(m=1, R=self.ra, r=([r[i]*np.cos(t[i]),r[i]*np.sin(t[i])]),v=([vx[i],vy[i]]),con=False)
            self.ball.append(a)    
        #setting initial variables
        self.dt_b = []
        self.dt_c = []
        self.b_i = 0
        self.b_j = 0
        self.c_i = 0
        
        #generating a list of dt values for ball-container collisions
        for i in range(self.num):
            self.dt_c.append(self.ball[i].time_to_collision(self.container))
        #generating an array of dt values for ball-to-ball collisions
        for i in range(self.num):
            a = []
            for j in range(self.num):
                if i != j:
                    a.append(self.ball[i].time_to_collision(self.ball[j]))
                else:
                    a.append(1000)
            self.dt_b.append(a)
        
        #finding minimum dt of ball-container collisions
        self.c_min = 1000
        for i in range(len(self.dt_c)):
            if 0 <= self.dt_c[i] and self.dt_c[i] < self.c_min:
                self.c_min = self.dt_c[i]
                self.c_i = i
        #finding minimum dt of ball-ball collisions
        self.b_min = 1000
        for i in range(len(self.dt_b)):
            for j in range(len(self.dt_b[i])):
                if (0 <= self.dt_b[i][j] and self.dt_b[i][j] < self.b_min):
                    self.b_min = self.dt_b[i][j]
                    self.b_i = i
                    self.b_j = j
    
    def next_collision(self):
        '''
        Finding the next collision to happen
        split into two situations of ball-ball and ball-container collisions
        '''
        if (0 <= self.b_min and self.b_min < self.c_min): #if next collision is ball-ball
            #self.cim = 0
            inim = self.ball[self.b_i].me() + self.ball[self.b_j].me()
            self.ball[self.b_i]._position = self.ball[self.b_i].move(self.b_min) #moving ball 1
            self.ball[self.b_j]._position = self.ball[self.b_j].move(self.b_min) #moving ball 2
            self.ball[self.b_i].collide(self.ball[self.b_j]) #colliding
            finm = self.ball[self.b_i].me() + self.ball[self.b_j].me()
            self.cim = finm - inim #change in momentum
            
            #find updated position for the rest of the balls
            for i in range(len(self.dt_b)):
                if i == self.b_i or i == self.b_j:
                    pass
                else:
                    self.ball[i]._position = self.ball[i].move(self.b_min)
            
            #finding the updated dt for the rest of the balls
            for i in range(len(self.dt_b)):
                if i != self.b_i and i != self.b_j:
                    #find new dt for ball-wall collisions for rest of the balls
                    self.dt_c[i] = self.dt_c[i] - self.b_min
                else:
                    #find new dt for ball-wall collisions for the collided balls
                    self.dt_c[i] = self.ball[i].time_to_collision(self.container)
                    
                for j in range(len(self.dt_b[i])):
                    if i == j:
                        self.dt_b[i][j] = 1000
                    else:
                        if j != self.b_j and i != self.b_i and i != self.b_j and j != self.b_i:
                            self.dt_b[i][j] = self.dt_b[i][j] - self.b_min
                            #finding the new dt of the two balls collided
                        else:
                            self.dt_b[i][j] = self.ball[i].time_to_collision(self.ball[j])

            
            #finding new b_min
            self.b_min = 1000
            for i in range(len(self.dt_b)):
                for j in range(len(self.dt_b[i])):
                    if (0 <= self.dt_b[i][j] and self.dt_b[i][j] < self.b_min):
                        self.b_min = self.dt_b[i][j]
                        self.b_i = i
                        self.b_j = j
            
            #finding new c_min
            self.c_min = 1000
            for i in range(len(self.dt_c)):
                if 0 <= self.dt_c[i] < self.c_min:
                    self.c_min = self.dt_c[i]
                    self.c_i = i
        
        else: #if next collision is ball-container
            inim = self.ball[self.c_i].me()
            self.ball[self.c_i]._position = self.ball[self.c_i].move(self.c_min)
            self.ball[self.c_i].collide(self.container)
            finm = self.ball[self.c_i].me()
            self.cim = finm - inim #change in momentum
            
            #find updated position for the rest of the balls
            for i in range(len(self.dt_b)):
                if i != self.c_i:
                    self.ball[i]._position = self.ball[i].move(self.c_min)
            
            #finding the updated dt for the rest of the balls
            #finding the new dt of the two balls collided
            for i in range(len(self.dt_c)):
                if i != self.c_i:
                    #new dt_c for rest of the balls
                    self.dt_c[i] = self.dt_c[i] - self.c_min
                else:
                    #new dt_c for the collided ball
                    self.dt_c[i] = self.ball[i].time_to_collision(self.container)
                for j in range(len(self.dt_c)):
                    if i == j:
                        self.dt_b[i][j] = 1000
                    else:
                        if i != self.c_i:
                            #new dt_b for rest of the balls
                            self.dt_b[i][j] = self.dt_b[i][j] - self.c_min
                            
                        else:
                            #new dt_b for the collided ball
                            self.dt_b[i][j] = self.ball[i].time_to_collision(self.ball[j])
            
            #finding new c_min
            self.c_min = 1000
            for i in range(len(self.dt_c)):
                if 0 <= self.dt_c[i] < self.c_min:
                    self.c_min = self.dt_c[i]
                    self.c_i = i
                    
            #finding new b_min
            self.b_min = 1000
            for i in range(len(self.dt_b)):
                for j in range(len(self.dt_b[i])):
                    if (0 <= self.dt_b[i][j] and self.dt_b[i][j] < self.b_min):
                        self.b_min = self.dt_b[i][j]
                        self.b_i = i
                        self.b_j = j
                        
    
    def run(self, num_frames, animate=False):
        
        bpatch = []
                    
        if animate:
            f = pl.figure()
            ax = pl.axes(xlim=(-self.container.radius, self.container.radius), ylim=(-self.container.radius, self.container.radius))
            cpatch = self.container.get_patch()
            ax.add_artist(cpatch)
            for i in range(self.num):
                bpatch.append(self.ball[i].get_patch())
                ax.add_patch(bpatch[i])
            bpatch[i].center = self.ball[i]._position
            self.next_collision()
    
        #(5)for creating a plot of pressure against temperature
        self.pre = []
        self.temp = []
        #(7)for creating a histogram of ball velocities
        self.histvals3 = []

        for frame in range(num_frames):
            if animate:
                for i in range(len(self.dt_c)):
                    bpatch[i].center = self.ball[i]._position
            self.next_collision()

            if animate:
                pl.pause(0.01)
        if animate:
            pl.show()

#%%
#Running the simulation
'''
c = Simulation(10, 15, 1, 5)
c.run(100, animate=True)
'''