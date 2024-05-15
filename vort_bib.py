"""
Created on Thu May  2 12:05:25 2024

@author: gauthier
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import os


######### Coordinates Y & Z of the wall image vortex ###########
def vortYposition(y_real, z_real, N_imag, L, Circ_real):
    
    if type(y_real)==float: y_real=np.array([y_real])
    if type(z_real)==float: z_real=np.array([z_real])

    y_images=np.zeros((len(y_real), 2*N_imag)) #Position of the vortex
    z_images=np.zeros((len(z_real), 2*N_imag)) #Position of the vortex
    Circ_images=np.zeros_like(y_images)        #Circulation of the vortex
    
    for j in range(len(y_real)):
        y_temp=np.array(())
        for i in range(1,N_imag+1):
            y_left=-2*i*L+(-1)**i*y_real[j]
            y_right=2*i*L+(-1)**i*y_real[j]
            y_temp=np.append(y_left, y_temp) #left-side vortex
            y_temp=np.append(y_temp,y_right) #right-side vortex
            Circ_images[j,int(len(Circ_images[j,:])/2)+i-1]=(-1)**i*Circ_real[j]
            z_images[j,:]=z_real[j]

        Circ_images[j,:int(len(Circ_images[j,:])/2)]=np.flip(Circ_images[j,int(len(Circ_images[j,:])/2):])
        print(len(y_temp))
        y_images[j]=y_temp
        
    print('\n\n Position of the vortices\n', y_images)    
    print('\n\n Circulation of the vortices\n', Circ_images)    
    
    return y_images, z_images, Circ_images


######### Coordinates Z of the ground images of the real and wall-imaged vortices ###########
def vortZposition(z_real,z_images, Circ_real, Circ_images):
    
    # Initialize the z and circulation arrays for ground images for both the real vortices and their lateral images
    if (type(z_images) == list) :
        z_ground_images=np.zeros((len(z_real)))
        Circ_ground_images=np.zeros((len(z_real)))
        y_ground_images=np.zeros((len(z_real)))
        
        z_ground_images[0:len(z_real)]= -1*z_real   # Vertical images of the real vortices
        Circ_ground_images[0:len(z_real)]= -1*Circ_real

    else: 
        z_ground_images=np.zeros((len(z_real)) + len(z_images[0,:])*len(z_images[:,0]))
        Circ_ground_images=np.zeros((len(z_real)) + len(z_images[0,:])*len(z_images[:,0]))
        
        z_ground_images[0:len(z_real)]= -1*z_real   # Vertical images of the real vortices
        Circ_ground_images[0:len(z_real)]= -1*Circ_real

        for i in range(len(z_real)):                # Vertical images of the horizontal images 
            z_ground_images[len(z_real)+i*len(z_images[0,:]) : len(z_real)+(i+1)*len(z_images[0,:]) ] = -1*z_images[i,:]
            Circ_ground_images[len(z_real)+i*len(z_images[0,:]) : len(z_real)+(i+1)*len(z_images[0,:]) ] = -1*Circ_images[i,:]      
            

    print('\n\n Position of the ground images\n', z_ground_images)
    print('\n\n Circulation of the ground images\n', Circ_ground_images)
    
    return z_ground_images, Circ_ground_images
        

def veloc_field(y_imag, z_imag, y_real, z_real, circulation, circulation_real, y_mesh, z_mesh, eps):

    '''U_theta=Gamma/(2.pi.(r-r_vortex)) --> Uy=Σ_i(-Gamma_i/(2.pi)*(z-zc_i)/((z-zc_i)^2+(y-yc_i)^2)
                                             Uz=Σ_i(Gamma_i/(2.pi)*(y-yc_i)/((z-zc_i)^2+(y-yc_i)^2)'''

    N_realVortex=len(circulation_real) #number of real vortex
    Uy=np.zeros_like(y_mesh)
    Uz=np.zeros_like(y_mesh)

    G=circulation/(2*math.pi) #reduced circulation for more clear reading
    G_real=circulation_real/(2*math.pi)   
    for i in range(len(y_mesh[:,0])):
        for j in range(len(y_mesh[0,:])):

            #Adding the contribution of the real vortices
            Uy[i,j]=np.sum(G_real*(z_mesh[i,j]-z_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2))\
                    +np.sum(-1*G_real*(z_mesh[i,j]+z_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]+z_real)**2))
            Uz[i,j]=np.sum(G_real*(y_mesh[i,j]-y_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2))\
                    +np.sum(-1*G_real*(y_mesh[i,j]-y_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]+z_real)**2))
            
            #Then the contribution of the lateral images
            for k in range(N_realVortex): #supposed to be the number of vortices
                Uy[i,j]+=np.sum(-1*G[k,:]*(z_mesh[i,j]-z_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))\
                        +np.sum(G[k,:]*(z_mesh[i,j]+z_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]+z_imag[k,:])**2))
                Uz[i,j]+=np.sum(G[k,:]*(y_mesh[i,j]-y_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))\
                        +np.sum(-1*G[k,:]*(y_mesh[i,j]-y_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]+z_imag[k,:])**2))

    # for i in range(len(y_mesh[:,0])):
    #     for j in range(len(y_mesh[0,:])):

    #         #Adding the contribution of the real vortices

    #         #Bad condition, does not take into account the eps radius back and forth


    #         Uy[i,j]=np.sum((z_mesh[i,j]**2+y_mesh[i,j]**2>=z_real**2+y_real**2)*G_real*(z_mesh[i,j]-z_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2)\
    #                       +(z_mesh[i,j]**2+y_mesh[i,j]**2<z_real**2+y_real**2)*G_real*(z_mesh[i,j]-z_real)/())\
    #                 +np.sum(-1*G_real*(z_mesh[i,j]+z_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]+z_real)**2))
    #         Uz[i,j]=np.sum(G_real*(y_mesh[i,j]-y_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]-z_real)**2))\
    #                 +np.sum(-1*G_real*(y_mesh[i,j]-y_real)/((y_mesh[i,j]-y_real)**2+(z_mesh[i,j]+z_real)**2))
            
    #         #Then the contribution of the lateral images
    #         for k in range(N_realVortex): #supposed to be the number of vortices
    #             Uy[i,j]+=np.sum(-1*G[k,:]*(z_mesh[i,j]-z_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))\
    #                     +np.sum(G[k,:]*(z_mesh[i,j]+z_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]+z_imag[k,:])**2))
    #             Uz[i,j]+=np.sum(G[k,:]*(y_mesh[i,j]-y_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]-z_imag[k,:])**2))\
    #                     +np.sum(-1*G[k,:]*(y_mesh[i,j]-y_imag[k,:])/((y_mesh[i,j]-y_imag[k,:])**2+(z_mesh[i,j]+z_imag[k,:])**2))


 #Come up with a core size for each vortices, additionnal feature ?

    return Uy, Uz




