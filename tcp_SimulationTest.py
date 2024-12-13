# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 23:07:42 2022

@author: aoc21
"""

from tcp_Simulation import Simulation
c = Simulation(10,2,1,5)
c.run(100, animate=True)