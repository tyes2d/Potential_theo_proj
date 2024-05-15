# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import os
import vort_bib


#======== DOMAIN =========

L=0.6    #Horizontal boundary
Ny=10
H=0.8    #Vertical boundary
Nz=10
h=0.16  #hub hieght
Cy=5 #dilatation factopr for the domain

y=np.linspace(-Cy*L,Cy*L,Cy* Ny+1)
z=np.linspace(H,-H,2* Nz+1)

y_mesh,z_mesh=np.meshgrid(y,z)


#=========Vortex============

i_v=np.argmin(np.abs(y_mesh[0,:]-0))
j_v=np.argmin(np.abs(z_mesh[:,0]-h))


y_real=np.array([0]) # position of the real vortices
z_real=np.array([0.5]) # position of the real vortices
circ_real=np.array([1]) # circulation of the real vortices 

###################################################################################
# Default value to trigger a condition: if type(circ/z_images)==list --> No walls #
###################################################################################
circ_images=[]
y_images=[]
z_images=[]



#getting the position and circulation of the image vortices
y_images, z_images, circ_images=vort_bib.vortYposition(y_real, z_real,3,L, circ_real)
z_ground_images, circ_ground_images = vort_bib.vortZposition(z_real,z_images, circ_real, circ_images)
# Getting the velocity indued by all the vortices
Uy, Uz=vort_bib.veloc_field(y_images, z_images, y_real, z_real ,circ_images, circ_real, y_mesh, z_mesh)

