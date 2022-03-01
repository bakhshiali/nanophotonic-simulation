from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
def energies(kx,ky,a,t):
    return t*np.sqrt(1+4*np.power(np.cos(np.sqrt(3)*kx*a*0.5),2)+4*np.cos(np.sqrt(3)*kx*a*0.5)*np.cos(3*ky*a*0.5))


fig = plt.figure(1)
#ax = plt.axes(projection='3d')
ax = fig.add_subplot(1, 1, 1, projection='3d')

a=1.42#/np.sqrt(3) #graphene carbon carbon distance in nm (0.142)
t=3.033#e.v
points=100
kx=np.linspace(-2, 2, points)
ky=np.linspace(-2, 2, points)
kx,ky=np.meshgrid(kx, ky)
energy=energies(kx,ky,a,t)

b1=[(2*np.pi)/(3*a),(2*np.pi*np.sqrt(3))/(3*a)]
b1value=np.sqrt(np.power(b1[0],2)+np.power(b1[1],2))
b2=[(2*np.pi)/(3*a),-1*(2*np.pi*np.sqrt(3))/(3*a)]

def hexagon(center_x,center_y,radius,height_z):
    z = np.linspace(9.5, height_z, 1)
    theta = np.linspace(0, 2*np.pi, 7)#n-1
    theta_grid, z_grid=np.meshgrid(theta, z)
    y_grid = radius*np.cos(theta_grid) + center_x
    x_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid
r1=np.sqrt(3)*(b1value*0.5)

Xc,Yc,Zc = hexagon(0,0,r1,10)
ax.plot_wireframe(Xc,Yc,Zc, alpha=1)
Xc1,Yc1,Zc1 = hexagon(0,0,b1value*0.5,10)#inner hexagon
ax.plot_wireframe(Yc1,Xc1,Zc1,alpha=1)

surf=ax.plot_surface(kx, ky, energy, rstride=1, cstride=1,lightsource=None,
                cmap=plt.get_cmap('rainbow'),  clim=(0,10),facecolor= None,
				antialiaseds= None,offsets= None,alpha=0.7)#,alpha=0.7)#edgecolor='none',
surf2=ax.plot_surface(kx, ky, -energy,rstride=1, cstride=1,lightsource=None,
                cmap=plt.get_cmap('rainbow_r'),  clim=(-10,0),facecolor= None,
				antialiaseds= None,offsets= None,alpha=0.7)

cset = ax.contour(kx, ky, energy, zdir='z', offset=-10, cmap=plt.get_cmap('rainbow'))
cset = ax.contour(kx, ky, energy, zdir='x', offset=-2.5, cmap=plt.get_cmap('rainbow'))
cset = ax.contour(kx, ky, energy, zdir='y', offset=2.5, cmap=plt.get_cmap('rainbow'))
cset = ax.contour(kx, ky, -energy, zdir='x', offset=-2.5, cmap=plt.get_cmap('rainbow'))
cset = ax.contour(kx, ky, -energy, zdir='y', offset=2.5, cmap=plt.get_cmap('rainbow'))

fig.colorbar(cset, shrink=0.5, aspect=5)
#fig.colorbar(surf2, shrink=0.5, aspect=5)
ax.set_title('Ali Bakhshi\nGraphene Energy Spectrum(1st Nearest Neighbor)');

ax.set_xlabel('kx (1/Å)')
ax.set_xlim(-2, 2)
ax.set_ylabel('ky (1/Å)')
ax.set_ylim(-2, 2)
ax.set_zlabel('energy (ev)')
ax.set_zlim(-10, 10)

#2D structure
k1=[(-1*b1value)/(2*np.cos(30*(np.pi/180))),0]
gamma=[0,0]
m=[0,b1value/2]
k2=[-1*k1[0]*np.cos(60*(np.pi/180)),m[1]]

ax.scatter(k1[0],k1[1],0,color='black')
ax.scatter(gamma[0],gamma[1],0,color='black')
ax.scatter(m[0],m[1],0,color='black')
ax.scatter(k2[0],k2[1],0,color='black')

K1Gamma=np.linspace(k1[0], 0, points)
GammaM=np.linspace(0,m[1],points)
MK2=np.linspace(0,k2[0],points)
energyK1Gamma=energies(K1Gamma,0,a,t)
energyGammaM=energies(0,GammaM,a,t)
energyMK2=energies(MK2,m[1],a,t)
energyspace=np.concatenate((energyK1Gamma,energyGammaM,energyMK2),axis=0)
kspace=np.linspace(k1[0], k2[0], 3*points)
kspace2=np.concatenate((K1Gamma,GammaM,MK2),axis=0)
kspace2.sort()
#ax.plot(kspace, energyspace, zs=-12, zdir='z', label='curve in (x,y)')

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1, 1, 1)

ax2.set_title('band structure of monolayer graphene (1st Nearest Neighbor)')
ax2.plot(kspace, energyspace)
ax2.plot(kspace, -energyspace)
ax2.set_ylabel('Energy (ev)', {'color': 'C0', 'fontsize': 12})
ax2.set_xlabel('k space', {'color': 'C0', 'fontsize': 12})

plt.xticks((k1[0], -k2[0],m[0],k2[0]), ('$K$', r'$\Gamma$', '$M$','$K$'))#, color='k', size=12)
plt.show()
