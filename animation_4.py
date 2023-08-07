#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 4,3
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

#FLAG_save = False
FLAG_save = True

r = 1 # radius of circle
def circle(k,t):
    global freq
    global phi0
    phi_k = phi0[k]+freq[k]*t
    #print("k=%d t=%f freq=%f    phi=%f"%(k,t,freq[k],phi_k))
    #phi_k = phi*(1.0+k0.1)
    return r*np.cos(phi_k), r*np.sin(phi_k)

Time = 15.0
dt = 0.1
num_time = int(Time/dt)
num_phases = 10
phi0 = np.zeros(num_phases) 
#freq = np.zeros(num_phases) 
freq = np.random.normal(1.0,0.1,num_phases)
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
for k in range(num_phases):
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
    
# Updating function, to be repeatedly called by the animation
def update(t):
    global dt
    # obtain point coordinates 
    str = ''
    for kk in range(num_phases):
        x[kk],y[kk]= circle(kk,t)
        # set point's coordinates
        str += '     %f %f '%(x[kk],y[kk])
        p[kk].set_data([x[kk]],[y[kk]])
        #p[kk].set_color(label_color[kk])
    #print(str+'\n')
    return p

# create animation with 10ms interval, which is repeated,
# provide the full circle (0,2pi) as parameters
num_loop = 1

anim = FuncAnimation(fig, update, interval=100, blit=True, repeat=False,
                    frames=np.linspace(0,num_loop*Time,num_loop*num_time, endpoint=False))
if FLAG_save==True:
    f = r"./animation.gif" 
    writergif = PillowWriter(fps=30) 
    anim.save(f, writer=writergif)
plt.show()
