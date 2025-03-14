
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import os
import vort_bib


#======== DOMAIN =========

L=0.6    #Horizontal boundary
Ny=100
H=0.8    #Vertical boundary
Nz=25
h=0.16  #hub hieght
D=0.157
Cy=5 #dilatation factopr for the domain
U_inf=10

y=np.linspace(-Cy*L,Cy*L,Cy* Ny+1)
z=np.linspace(-H,H,2* Nz+1)
y_streamplot=np.linspace(-Cy*L,Cy*L,Cy* Ny+1)
z_streamplot=np.linspace(-H,H,2* Nz+1)

# ############### TO BE COMMENTED IF REGULAR MESH ##################
# y=np.tan(y)
# y = -Cy*L + (y - np.min(y)) * (2*Cy*L) / (np.max(y) - np.min(y))
# z=np.tan(z)
# z = -1*H + (z - np.min(z)) * (2*H) / (np.max(z) - np.min(z))
# ##################################################################


y_mesh,z_mesh=np.meshgrid(y,z)


#=========Vortex============

i_v=np.argmin(np.abs(y_mesh[0,:]-0))
j_v=np.argmin(np.abs(z_mesh[:,0]-h))


y_real=np.array([0, 0]) # position of the real vortices
z_real=np.array([h+D/4, h-D/4]) # position of the real vortices
circ_real=np.array([1, -1]) # circulation of the real vortices 

###################################################################################
# Default value to trigger a condition: if type(circ/z_images)==list --> No walls #
###################################################################################
circ_images=[]
y_images=[]
z_images=[]
eps=0.01 # core size of the vortex, for display purposes

print('====================================================================')
print('==================== Wind energy course project ====================')
print('=============== Derek Micheletto, Gauthier Leclercq ================')
print('====================================================================')

#getting the position and circulation of the image vortices
y_images, z_images, circ_images=vort_bib.vortYposition(y_real, z_real,3,L, circ_real)
z_ground_images, circ_ground_images = vort_bib.vortZposition(z_real,z_images, circ_real, circ_images)
# Getting the velocity indued by all the vortices
Uy, Uz=vort_bib.veloc_field(y_images, z_images, y_real, z_real ,circ_images, circ_real, y_mesh, z_mesh, eps)



############# Figures ###############

plt.figure(1, figsize=(15,4))
plt.contourf(y,z,np.sqrt(Uy**2+Uz**2) ,40) # Warning, log scale
plt.streamplot(y_streamplot,z_streamplot,Uy,Uz, density=(10,5), color='white', linewidth=0.5) # Warning, log scale
plt.hlines(0,-L,L,colors='k',linestyle='solid',linewidth=2)
plt.vlines(-L,0,H,colors='k',linestyle='solid',linewidth=2)
plt.vlines(L,0,H,colors='k',linestyle='solid',linewidth=2)
plt.xlim(left=-0.5*L, right=0.5*L)
plt.ylim(bottom=0, top=0.5*H)
plt.colorbar()



# plt.figure(1, figsize=(15,4))
# plt.contourf(y,z,np.sqrt(Uy**2+Uz**2),40, norm='log') # Warning, log scale
# plt.hlines(0,-L,L,colors='k',linestyle='solid',linewidth=2)
# plt.vlines(-L,0,H,colors='k',linestyle='solid',linewidth=2)
# plt.vlines(L,0,H,colors='k',linestyle='solid',linewidth=2)
# plt.colorbar()



# plt.figure(2,figsize=(15,4))

# Uy_quiver=Uy
# Uz_quiver=Uz
# Umagn=np.sqrt(Uy**2+Uz**2)

# threshold=20

# Uy_quiver[Umagn>threshold]=0
# Uz_quiver[Umagn>threshold]=0

# plt.quiver(y,z, Uy,Uz)



######### TIME EVOLUTION ###########

tf,Nt=10,1000
dt=tf/Nt
tolerance_detectionY, tolerance_detectionZ=L/Ny, H/(Nz)

for i in range(Nt):
	index_update=np.zeros([2,len(y_real)])

	for k in range(len(y_real)): 	#to find the index of the real vorteicces
		index_updateY=np.where(np.abs(y_mesh-y_real[k])<=tolerance_detectionY)[1][0]
        #if(np.where(np.abs(z_mesh-z_real[k])<=tolerance_detectionZ).size==0):break
		index_updateZ=np.where(np.abs(z_mesh-z_real[k])<=tolerance_detectionZ)[0][0]
		y_real[k]=y_mesh[0,index_updateY]+Uy[0,index_updateY]*i*dt
		z_real[k]=z_mesh[index_updateZ,0]+Uz[index_updateZ,0]*i*dt

	y_images, z_images, circ_images=vort_bib.vortYposition(y_real, z_real,3,L, circ_real)
	z_ground_images, circ_ground_images = vort_bib.vortZposition(z_real,z_images, circ_real, circ_images)
	# Getting the velocity induced by all the vortices
	Uy, Uz=vort_bib.veloc_field(y_images, z_images, y_real, z_real ,circ_images, circ_real, y_mesh, z_mesh, eps)



	plt.figure(i+100, figsize=(15,4))
	plt.title('X position: '+str(U_inf*i*dt))
	plt.contourf(y,z,np.sqrt(Uy**2+Uz**2) ,np.linspace(np.min(np.sqrt(Uy**2+Uz**2)), np.max(np.sqrt(Uy**2+Uz**2)), num=50), extend="both") # Warning, log scale
	#plt.streamplot(y_streamplot,z_streamplot,Uy,Uz, density=(10,5), color='white', linewidth=0.5) # Warning, log scale
	plt.hlines(0,-L,L,colors='k',linestyle='solid',linewidth=2)
	plt.vlines(-L,0,H,colors='k',linestyle='solid',linewidth=2)
	plt.vlines(L,0,H,colors='k',linestyle='solid',linewidth=2)
	# plt.xlim(left=-0.5*L, right=0.5*L)
	# plt.ylim(bottom=0, top=0.5*H)
	plt.colorbar()
	print (i)



# for i in range(Nt):
# 	for j in range(len(y_real)):
# 		z_real[j]=Uz[y_real==y_mesh[j] and z_real==z_mesh[j]]*i*dt
# 	y_images, z_images, circ_images=vort_bib.vortYposition(y_real, z_real,3,L, circ_real)
# 	z_ground_images, circ_ground_images = vort_bib.vortZposition(z_real,z_images, circ_real, circ_images)
# 	# Getting the velocity indued by all the vortices
# 	Uy, Uz=vort_bib.veloc_field(y_images, z_images, y_real, z_real ,circ_images, circ_real, y_mesh, z_mesh, eps)



