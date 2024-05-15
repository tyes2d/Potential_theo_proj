#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:05:25 2024

@author: gauthier
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import os


def vortYposition(y_real, z_real, N_imag, L, circ):
    
    if type(y_real)==float: y_real=np.array([y_real])
    if type(z_real)==float: z_real=np.array([z_real])

    y=np.zeros((len(y_real), 2*N_imag)) #Position of the vortex
    z=np.zeros((len(z_real), 2*N_imag)) #Position of the vortex
    Circ_images=np.zeros_like(y)        #Circulation of the vortex
    
    for j in range(len(y_real)):
        y_temp=np.array(())
        for i in range(N_imag):
            y_left=-2*i*L+(-1)**i*y_real[j]
            y_right=2*i*L+(-1)**i*y_real[j]
            y_temp=np.append(y_left, y_temp) #left-side vortex
            y_temp=np.append(y_temp,y_right) #right-side vortex
            Circ_images[j,int(len(Circ_images[j,:])/2)+i]=(-1)**i*circ[j]
            z[j,:]=-1*z_real[j]

        Circ_images[j,:int(len(Circ_images[j,:])/2)]=np.flip(Circ_images[j,int(len(Circ_images[j,:])/2):])
        print(len(y_temp))
        y[j]=y_temp
        
    print('\n\n Position of the vortices\n', y)    
    print('\n\n Circulation of the vortices\n', Circ_images)    
    
    return y, z, Circ_images
        

def veloc_field(y_imag, z_imag, y_real, z_real, circulation, circulation_real, y_mesh, z_mesh):

    '''U_theta=Gamma/(2.pi.(r-r_vortex)) --> Uy=Σ_i(-Gamma_i/(2.pi)*(z-zc_i)/((z-zc_i)^2+(y-yc_i)^2)
                                             Uz=Σ_i(Gamma_i/(2.pi)*(y-yc_i)/((z-zc_i)^2+(y-yc_i)^2)'''

    N_realVortex=len(circulation_real) #number of real vortex

    Uy=np.zeros_like(y_mesh)
    Uz=np.zeros_like(y_mesh)

    G=circulation/(2*math.pi) #reduced circulation for more clear reading

    for i in range(len(y_mesh[:,0])):
        for j in range(len(y_mesh[0,:])):

            #Adding the contribution of the real vortices
            Uy[i,j]=np.sum(circulation_real*(z_mesh[i,j]-z_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2))
            Uz[i,j]=np.sum(circulation_real*(y_mesh[i,j]-y_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2))
            
            #Then the contribution of the image of the real vortices
            for k in range(N_realVortex): #supposed to be the number of vortices
                Uy[i,j]+=np.sum(-1*G[k,:]*(z_mesh[i,j]-z_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))
                Uz[i,j]+=np.sum(G[k,:]*(y_mesh[i,j]-y_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))


    # Uy=np.sum(circulation_real*(z_mesh-z_real)/((y_mesh-y_real)**2+(z_mesh-z_real)**2))
    # Uz=np.sum(circulation_real*(y_mesh-y_real)/((y_mesh-y_real)**2+(z_mesh-z_real)**2))

    # for k in range(N_realVortex): #supposed to be the number of vortices
    #     Uy+=np.sum(G[k,:]*(z_mesh-z_imag[k,:])/((y_mesh-y_imag[k,:])**2+(z_mesh-z_imag[k,:])**2))
    #     Uz+=np.sum(G[k,:]*(y_mesh-y_imag[k,:])/((y_mesh-y_imag[k,:])**2+(z_mesh-z_imag[k,:])**2))




    return Uy, Uz

'''Problem is the triple loop, a bit heavy '''



    #for k in range(len(position[:,0])):

