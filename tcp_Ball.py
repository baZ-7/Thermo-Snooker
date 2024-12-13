# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 10:00:52 2022

@author: Andrew Chim aoc21
>>The Thermo Snooker Computing Project<<
File 1: Ball Class

This code file contains the Ball class:
For creating a 'ball' object and a 'container' object for the simulation
    
"""
#%%
#imports
import numpy as np
import pylab as pl 

#%%
#Class Ball
class Ball:
    
    def __init__(self, m=0, R=0.0, r=np.array([0.0,0.0]), v=np.array([0.0,0.0]), con=False):
        '''
        initiating the class Ball
        ===== variables =====
        m --- mass, mass of the object, set as 0 initiated but set as 1 throughout the simulation
        R --- radius, radius of the object, set as 0 when initiated and can be varied
        r --- position 2D-vector, set as 0 when initiated, hidden variable for protection
        v --- velocity 2D-vector, set as 0 when initiated
        con --- determining whether the object is a container, i.e. false = not a container & true = container
        =====================
        '''
        self.radius = R
        self.mass = m
        self._position = r
        self.velocity = v
        #setting the object as a container
        self.container = con 
        #container velocity 
        self.con_vel = np.array([0.0,0.0])
        if self.container == True:
            #setting the container a mass of ~infinity
            #self.mass = 2**63
            #setting initial collision point
            self.col_point = ([0,0])
            
        
        else:
            return
            
        
    def pos(self):
        return self._position
    
    def vel(self):
        return self.velocity
    
    def move(self, dt):
        r_new = ([np.multiply(self.velocity[0],dt) + self._position[0],np.multiply(self.velocity[1],dt) + self._position[1]])
        return r_new
    
    def time_to_collision(self, other):
        '''
        calculate the time until the next collision between this ball and another one,or the container.
        using the equation provided
        '''
        if other.container == False: #ball-ball collisions
            R_di = self.radius + other.radius #sum of radius
            r_di = ([self._position[0] - other._position[0],self._position[1] - other._position[1]]) #difference in position
            v_di = ([self.velocity[0] - other.velocity[0],self.velocity[1] - other.velocity[1]]) #difference in velocity
            rv = np.dot(r_di,v_di) 
            vv = np.dot(v_di,v_di)
            rr = np.dot(r_di,r_di)
            if vv == 0:
                dt = 1000
            else:
                if (rv)**2 - vv*((rr) - (R_di)**2) >= 0:    
                    dt = (- rv - np.sqrt((rv)**2 - vv*((rr) - (R_di)**2)))/(vv)
                else:
                    dt = 1000 #arbitrary large value
                if dt < 0:
                    dt = 1000 #arbirary large value
           
        else: #ball-container collisions
            #case when the object is indeed a container:
            Ra = - self.radius + other.radius
            #finding unit vector form of velocity
            unit_v = (self.velocity/np.sqrt(self.velocity[0]**2 + self.velocity[1]**2))
            a = unit_v[0]**2 + unit_v[1]**2 #a,b,c of the quadratic formula
            b = 2*(unit_v[0]*self._position[0]+unit_v[1]*self._position[1])
            c = self._position[0]**2 + self._position[1]**2 - Ra**2
            arb1 = (-b+np.sqrt(b**2-4*a*c))/2*a #arbitrary constant
            #finding collision point
            self.col_point = arb1*unit_v + self._position
            #time to collision
            r = np.sqrt((self.col_point[0] - self._position[0])**2 + (self.col_point[1] - self._position[1])**2)
            v = np.sqrt(self.velocity[0]**2 + self.velocity[1]**2)
            dt = r / v
        return dt
    
    def collide(self, other):
        '''
        make the changes to the velocities of the ball and the other one due to a collision.
        separate the components of the velocity of the two balls:
        parallel and perpendicular to the line of centers between the two balls
        compute parallel component: 1D Newtonian mechanics collision equation
        compute perpendicular component: unchanged
        '''
        inim = 0 #initial momentum
        finm = 0 #final momentum
        
        if other.container == False: #ball-ball collisions
            r_di = ([self._position[0] - other._position[0],self._position[1] - other._position[1]]) #difference in position
            rm = np.sqrt([(self._position[0] - other._position[0])**2 + (self._position[1] - other._position[1])**2]) #mod of r_di
            r_diu = r_di/rm #difference in position unit vector form
            a_par = np.multiply(([np.dot(self.velocity,r_di), np.dot(self.velocity,r_di)])/rm,r_diu)
            #finding component of velocity of the ball parallel to the line of centers
            b_par = np.multiply(([np.dot(other.velocity,r_di), np.dot(other.velocity,r_di)])/rm,r_diu)
            #finding component of velocity of the other ball parallel to the lines of centers
            a_per = ([self.velocity[0] - a_par[0], self.velocity[1] - a_par[1]]) #finding the perpendicular component
            b_per = ([other.velocity[0] - b_par[0], other.velocity[1] - b_par[1]])
        
            #using the 1-D Newtonian equation
            a = np.array(a_par)#np.sqrt(np.dot(a_par,a_par)) #finding magnitude of parallel component
            b = np.array(b_par)#np.sqrt(np.dot(b_par,b_par))
            ma = self.mass
            mb = other.mass 
            v_1_x = ((ma - mb)*a[0]/(ma + mb)) + ((2*mb)*b[0]/(ma + mb))
            v_2_x = ((2*ma)*a[0]/(ma + mb)) + ((mb - ma)*b[0]/(ma + mb))
            v_1_y = ((ma - mb)*a[1]/(ma + mb)) + ((2*mb)*b[1]/(ma + mb))
            v_2_y = ((2*ma)*a[1]/(ma + mb)) + ((mb - ma)*b[1]/(ma + mb))
            
            self.velocity = ([v_1_x + a_per[0],v_1_y + a_per[1]]) 
            other.velocity = ([v_2_x + b_per[0],v_2_y + b_per[1]]) 
            
            self.cim = 0 #no change in momentum for ball to ball
        
        else: #ball-container collisions
            ini_velocity = self.velocity
            Ra = - self.radius + other.radius
            r_diu = self.col_point/np.sqrt(self.col_point[0]**2 + self.col_point[1]**2)
            a_par = np.dot(self.velocity,r_diu)*r_diu
            a_per = self.velocity - a_par
            self.velocity = - a_par + a_per
            
    
    def ke(self):
        '''
        calculating the KE for a given velocity
        0.5*m*v**2
        '''
        v = np.sqrt(self.velocity[0]**2 + self.velocity[1]**2)
        return 0.5*self.mass*(v**2)
    
    def me(self):
        '''
        calculating the momentum for a given velocity
        '''
        me = np.sqrt(self.velocity[0]**2 + self.velocity[1]**2)*self.mass
        return me
    
    def get_patch(self):
        '''
        obtaining patch for the ball and container
        '''
        if self.container == False: #ball
            return pl.Circle(self._position, self.radius, fc='r')
        else: #container
            return pl.Circle([0., 0.], self.radius, ec='b', fill=False, ls='solid')

#%%

        
    
    
