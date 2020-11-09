# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 11:14:39 2020

@author: amill
"""

from sympy import *
import numpy as np 
import random 
import matplotlib.pyplot as plt 
import math as math
import time 
import matplotlib.pyplot as plt 

lx = 20
#size of simulation box for MC 
simbox_size = np.array([lx, lx])
p = 0.1
sigma = 1
epsilon = 2
T = 1
Nt = 50 #number of time steps 

#number of particles in or box. They need to have density 0.95. Density is mass over volume 
number_particles = int(p*lx*lx)

#array which randomly assigns the coordiantes of each particle 
particle_coords = [(random.uniform(0, lx), random.uniform(0, lx))]

#radius function
def radius(x, y):
    r = math.sqrt((x**2)+(y**2))
    return r 

#while loop that generates particles and gives them coordinates iif they dont overlap within a certain radius, THIS INTIALISES THE POSITIONS 
#NEED TO DO BOUNDARY CONDITIONS 
counter = 1
while counter < number_particles:
    #automatically sets if particles are overlapping to false 
    overlapping = False
    #generates random x and y coordinates for a possible new particle
    x = random.uniform(0, lx)
    y = random.uniform(0, lx)
    #for loop to check if the new particles radius is below the value we need it to be (sigma)
    for i,j in particle_coords:
        if abs(radius(x-i, y-j)) < sigma:
            overlapping = True 
    #if the particles are not overlapping then it adds its coordiates to the list and increases the counter 
    if overlapping == False:
        particle_coords.append([x, y])
        counter += 1
    
array_particle_coords = np.array(particle_coords)


#setting up random veloctiy of particle. PART OF INTIIALISING THE SYSTEM 
particle_velocities = []
for i in range(len(array_particle_coords)): 
    v_x = random.uniform(-1, 1)
    v_y = random.uniform(-1, 1)
    particle_velocities.append([v_x, v_y])

particle_velocities = np.array(particle_velocities)

#trying to find average velocities across all particles in x direction 
v_x_sum = 0
for i in particle_velocities:
    v_x_sum += i[0]
    
v_x_avg = v_x_sum/len(particle_velocities)

#trying to find average velocities across all particles in y direction 
v_y_sum = 0
for i in particle_velocities:
    v_y_sum += i[1]
    
v_y_avg = v_y_sum/len(particle_velocities)

#which average do they want us to use for the F_s function 
v_avg_square = ((v_y_sum)**2 + (v_x_sum)**2)/(2*len(particle_coords))

#equation for the drift factor given in the assessment 
f_s = math.sqrt((3*T)/(v_avg_square))

#adjusting the particle velocity by the factor given in the lecture to remove intiial dift 
adjusted_initial_velocities = []
for i in particle_velocities:
    adjusted_velocity_x = (i[0]-v_x_avg) * f_s
    adjusted_velocity_y = (i[1]-v_y_avg) * f_s
    adjusted_initial_velocities.append([adjusted_velocity_x, adjusted_velocity_y])

#constants needed for lennart jones 
r = 1
sigma = 1
alpha = 1 

#lennart jones potential equation 
U = 4*epsilon*(((sigma/r)**12)+alpha*((sigma/r)**6))

#differential of lennart jones, just to see if I can do it 
r = Symbol('r')
U = 4*epsilon*(((sigma/r)**12)+alpha*((sigma/r)**6))
Uprime = U.diff(r) #this is olllllldddddd news now with the variable being reassigned below.     

#for loop to calculate the force interactions between different particles
#value of r for which force basically becomes zero 
#sum of the forces in x direction 
force_sum_x = 0
#sum of forces in y direction 
force_sum_y = 0
r_cutoff = 2.5
#list of particle forces in x and y direction 
force_x = np.zeros(number_particles)
force_y = np.zeros(number_particles)
for i in range(number_particles):
    for j in range(i+1, number_particles): 
        atom1 = particle_coords[i]
        atom2 = particle_coords[j]
        #boundary conditions 
        #distance in the x coords of atoms 
        dx = atom1[0] - atom2[0]
        #distance in the y coords of atoms 
        dy = atom1[1] - atom2[1]
        #if the value of dx is smaller when we pass it through the boundary then change it 
        if lx - dx%lx < dx:
            dx = abs(lx-(dx%lx))
        #same witht eh y condition 
        if lx - dy%lx < dy:
            dy = abs((lx-(dy%lx)))
        #if the radius is smaller than the cut off (lennard jones potential is basically 0 after 2.5)
        if ((dx**2) + (dy**2))**0.5 < r_cutoff:
            r = (dx**2 + dy**2)**0.5
            #differential of lennard jones potential 
            F = -4*epsilon*(((12*(sigma**12))/r**13) - (((alpha*6*(sigma**6))/r**7)))
            #unit values for each direction so we can seperate force into x and y units in the list. 
            x_unit = (atom1[0]-atom2[0])/r
            y_unit = (atom1[1]-atom2[1])/r
            force_x[i] += F*x_unit
            force_y[i] += F*y_unit
            force_x[j] -= F*x_unit
            force_y[j] -= F*y_unit
    
m=1 
#kinetic energy of the system 
sum_velocity_x = 0
for i in adjusted_initial_velocities:
    sum_velocity_x += i[0]**2 

sum_velocity_y = 0
for i in adjusted_initial_velocities:
    sum_velocity_y += i[1]**2

total_kinetic = 0.5*m*((sum_velocity_x + sum_velocity_y))     

#temperature for the system
boltzmann = 1.0
#euqation for temperature of a system 
T_new = (2*total_kinetic)/(3*boltzmann)

#potential energy of the system 
potential_energy = 0
for i in range(number_particles):
    for j in range(number_particles):
        if i == j:
            continue
        else: 
            atom1 = particle_coords[i]
            atom2 = particle_coords[j]
            #boundary conditions 
            #distance in the x coords of atoms 
            dx = atom1[0] - atom2[0]
            #distance in the y coords of atoms 
            dy = atom1[1] - atom2[1]
            #if the value of dx is smaller when we pass it through the boundary then change it 
            if lx - dx%lx < dx:
                dx = abs(lx-(dx%lx))
            #same witht eh y condition 
            if lx - dy%lx < dy:
                dy = abs((lx-(dy%lx)))
            #if the radius is smaller than the cut off (lennard jones potential is basically 0 after 2.5)
            if ((dx**2) + (dy**2))**0.5 < r_cutoff:
                r = (dx**2 + dy**2)**0.5
                #equation for lennard jones potential energy  
                U_potential = 4*epsilon*((sigma/r)**12 - alpha*((sigma/r)**6))
                potential_energy += U_potential
 

#total energy of the system 
total_energy = total_kinetic + potential_energy

#constants needed for iterations, N is number of time steps and dt is the change in time
num_time_steps = 100
dt = 0.01

#these are the arrays we will add stuff to. We can use hstack method to add a new array to the end of an array 
particles = array_particle_coords
velocities = np.array(adjusted_initial_velocities)
forces = np.array([force_x, force_sum_y])

x_pos_particle_1 = []
for i in particle_coords:
    x_pos_particle_1.append(i[0])

x_pos_particle_1 = np.array(x_pos_particle_1)

y_pos_particle_1 = []
for i in particle_coords:
    y_pos_particle_1.append(i[1])

y_pos_particle_1 = np.array(y_pos_particle_1)

x_velocity = []
for i in velocities:
    x_velocity.append(i[0])

x_velocity = np.array(x_velocity)

y_velocity = []
for i in velocities:
    y_velocity.append(i[1])

y_velocity = np.array(y_velocity)

x_pos_particle_2 = []
for i in range(len(x_pos_particle_1)):
    new_x_position = x_pos_particle_1[i] + x_velocity[i]*dt + (0.5*force_x[i]*(dt**2))/m
    x_pos_particle_2.append(new_x_position)

y_pos_particle_2 = []
for i in range(len(y_pos_particle_1)):
    new_y_position = y_pos_particle_1[i] + y_velocity[i]*dt + (0.5*force_y[i]*(dt**2))/m
    y_pos_particle_2.append(new_y_position)

#===================================================================================================================
#ALL BETWEEN THE EQUALS SIGNS NEEDS TO BE PUT IN A WHILE LOOP TO ITERATE OVER IT
#have a look at the puesdocode but this is the first stage we need to put in a loop 
counter = 2 #this is because we have already calculated the first couple of positions before, so technically the first run of the loop is on the third time step 
#this is a frequency counter asked for in the paper 
f_log = 5 
while counter < num_time_steps:
    positions_array_for_force_calc = np.array((x_pos_particle_2, y_pos_particle_2)).T
    force_x = np.zeros(number_particles)
    force_y = np.zeros(number_particles)
    for i in range(number_particles):
        for j in range(number_particles):
            if i == j:
                continue
            else: 
                atom1 = positions_array_for_force_calc[i]
                atom2 = positions_array_for_force_calc[j]
                #boundary conditions 
                #distance in the x coords of atoms 
                dx = atom1[0] - atom2[0]
                #distance in the y coords of atoms 
                dy = atom1[1] - atom2[1]
                #if the value of dx is smaller when we pass it through the boundary then change it 
                if lx - dx%lx < dx:
                    dx = abs(lx-(dx%lx))
                #same witht eh y condition 
                if lx - dy%lx < dy:
                    dy = abs((lx-(dy%lx)))
                #if the radius is smaller than the cut off (lennard jones potential is basically 0 after 2.5)
                if ((dx**2) + (dy**2))**0.5 < r_cutoff:
                    r = (dx**2 + dy**2)**0.5
                    #differential of lennard jones potential 
                    F = -4*epsilon*(((12*(sigma**12))/r**13) - (((alpha*6*(sigma**6))/r**7)))
                    #unit values for each direction so we can seperate force into x and y units in the list. 
                    x_unit = (atom1[0]-atom2[0])/r
                    y_unit = (atom1[1]-atom2[1])/r
                    force_x[i] += F*x_unit
                    force_y[i] += F*y_unit
    
    #then we will use this function to calculate the new positions based off the 2 previous lists. 
    x_pos_particle_3 = []
    y_pos_particle_3 = []
    for i in range(number_particles):
        new_x_position = 2*x_pos_particle_2[i] - x_pos_particle_1[i] + (force_x[i]/m)*(dt**2)
        x_pos_particle_3.append(new_x_position)
        new_y_position = 2*y_pos_particle_2[i] - y_pos_particle_1[i] + (force_y[i]/m)*(dt**2)
        y_pos_particle_3.append(new_y_position)
        
    #then we create the new velocity 
    for i in range(number_particles):
        new_velocity_x = (x_pos_particle_3[i] - x_pos_particle_1[i])/(2*dt)
        x_velocity[i] = new_velocity_x
        new_velocity_y = (y_pos_particle_3[i] - y_pos_particle_1[i])/(2*dt)
        y_velocity[i] = new_velocity_y
    
    #kinetic energy of the system 
    sum_velocity_x = 0
    for i in x_velocity:
        sum_velocity_x += i**2 
    
    sum_velocity_y = 0
    for i in y_velocity:
        sum_velocity_y += i**2
    
    total_kinetic = 0.5*m*((sum_velocity_x + sum_velocity_y))     
    
    #temperature for the system
    boltzmann = 1.0
    #euqation for temperature of a system 
    T_new = (2*total_kinetic)/(3*boltzmann)
    
    #potential energy of the system 
    potential_energy = 0
    for i in range(number_particles):
        for j in range(number_particles):
            if i == j:
                continue
            else: 
                atom1 = positions_array_for_force_calc[i]
                atom2 = positions_array_for_force_calc[j]
                #boundary conditions 
                #distance in the x coords of atoms 
                dx = atom1[0] - atom2[0]
                #distance in the y coords of atoms 
                dy = atom1[1] - atom2[1]
                #if the value of dx is smaller when we pass it through the boundary then change it 
                if lx - dx%lx < dx:
                    dx = abs(lx-(dx%lx))
                #same witht eh y condition 
                if lx - dy%lx < dy:
                    dy = abs((lx-(dy%lx)))
                #if the radius is smaller than the cut off (lennard jones potential is basically 0 after 2.5)
                if ((dx**2) + (dy**2))**0.5 < r_cutoff:
                    r = (dx**2 + dy**2)**0.5
                    #equation for lennard jones potential energy  
                    U_potential = 4*epsilon*((sigma/r)**12 - alpha*((sigma/r)**6))
                    potential_energy += U_potential
    
    #total energy of the system 
    total_energy = total_kinetic + potential_energy
    counter += 1
    if counter%f_log == 0:
        print('KE =', total_kinetic)
        print('PE =', potential_energy)
        print('Total Energy =', total_kinetic + potential_energy)
        print('The new Temperature of the system is', T_new)
    
    #then we need to reassign the names of the lists so they get iterated between over and over again 
    x_pos_particle_1 = x_pos_particle_2 
    y_pos_particle_1 = y_pos_particle_2 
    x_pos_particle_2 = x_pos_particle_3 
    y_pos_particle_2 = y_pos_particle_3
    
#===================================================================================================================

#this sets up our scatter plot which will show the initial state of the particles, this may need to change soon though as the particles evolve. 
x = array_particle_coords[:,0]
y = array_particle_coords[:,1]

plt.plot(x,y, 'o')
              