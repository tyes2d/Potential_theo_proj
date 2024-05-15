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

#getting the position and circulation of the image vortices
Y_Vortex, Z_Vortex, Circ_Vortex=vort_bib.vortYposition(y_real, z_real,3,L, circ_real)
Z_images, Circ_images = vort_bib.vortZposition(z_real,Z_Vortex, circ_real, Circ_Vortex)
# Getting the velocity indued by all the vortices
# Uy, Uz=vort_bib.veloc_field(Y_Vortex, Z_Vortex, y_real, z_real ,Circ_Vortex, circ_real, y_mesh, z_mesh)

