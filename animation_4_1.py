#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 4,3
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

#FLAG_save = False
FLAG_save = True

r = 1 # radius of circle
def compute_phase(k,t):
    global freq
    global phi0
    return phi0[k]+freq[k]*t

def circle(k,t):
    phi_k = compute_phase(k,t)
    #print("k=%d t=%f freq=%f    phi=%f"%(k,t,freq[k],phi_k))
    #phi_k = phi*(1.0+k0.1)
    return r*np.cos(phi_k), r*np.sin(phi_k)

Time = 15.0
dt = 0.1
num_time = int(Time/dt)
num_phases = 10
num_phasepairs = num_phases*(num_phases-1)//2
phi0 = np.zeros(num_phases) 
freq = np.random.normal(1.0,0.03,num_phases)
#freq = np.zeros(num_phases) 
#freq[0:5]=1.0
#freq[5:10]=1.5
phi0 = np.random.normal(2.0,0.1,num_phases)
#phi0[0] = 0.0
#phi0[1] = 0.5
#phi0[2] = 1.0
#phi0[3] = 2.0

# create a figure with an axes
fig, ax = plt.subplots()
# set the axes limits
ax.axis([-1.5,1.5,-1.5,1.5])
# set equal aspect such that the circle is not shown as ellipse
ax.set_aspect("equal")
# create a point in the axes
p = []
for k in range(num_phasepairs):
    point, = ax.plot(1,0, marker="o")
    p.append(point)
#point1, = ax.plot(1,0, marker="*")
#p.append(point1)
#point2, = ax.plot(1,0, marker=".")
#p.append(point2)
#point2, = ax.plot(1,0, marker=".")
#p.append(point2)

#point1, = ax.plot(0,1, marker="*")
num = 100
label_color = ['k','r','g','b','y']
x0 = np.zeros(num)
y0 = np.zeros(num)
x = np.zeros(num_phases)
y = np.zeros(num_phases)
for k in range(num):
    phi = 2*np.pi*float(k)/float(num)
    x0[k] = np.cos(phi)
    y0[k] = np.sin(phi)
ax.plot(x0,y0,'k',alpha=0.2)
for pp in p: 
    pp.set_data([],[])


dphase = np.zeros((num_phasepairs,num_time))
for i in range(num_time):
    t = i * dt
    counter = 0
    for k in range(num_phases-1):
        phi_k = compute_phase(k,t)
        for l in range(num_phases-k-1):
            ll = k + l + 1
            phi_l = compute_phase(ll,t) 
            dphase[counter,i]=phi_k-phi_l
            counter += 1
# Updating function, to be repeatedly called by the animation
def update(n):
    global dphase
    # obtain point coordinates 
    str = ''
    for kk in range(num_phasepairs):
        #x[kk],y[kk]= circle(kk,t)
        x=np.cos(dphase[kk,n])
        y=np.sin(dphase[kk,n])
        # set point's coordinates
        str += '     %f %f '%(x,y)
        p[kk].set_data([x],[y])
        #p[kk].set_color(label_color[kk])
    #print(str+'\n')
    return p

# create animation with 10ms interval, which is repeated,
# provide the full circle (0,2pi) as parameters
anim = FuncAnimation(fig, update, interval=100, blit=True, repeat=False,
                    frames=range(num_time))
if FLAG_save==True:
    f = r"./animation.gif" 
    writergif = PillowWriter(fps=30) 
    anim.save(f, writer=writergif)
plt.show()
