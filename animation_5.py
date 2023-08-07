#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 4,3
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter

#FLAG_save = False
FLAG_save = True

class phases():
    def __init__(self):
        num_times = 0
        #self.in_label = './unknown_timeseries_new.dat'
        self.in_label = './unknown_timeseries_phases.dat'
        
    def read_in(self):
        data0 = np.loadtxt(self.in_label)
        self.num_times = np.shape(data0)[0]
        self.num_phases = 2#np.shape(data0)[1]-1
        self.Time = data0[-1:0]
        self.data = data0[:,1:self.num_phases+1]
        self.dt = data0[1,0]-data0[0,0]
        
P = phases()
P.read_in()

print("num_phases: %d"%(P.num_phases))
r = 1 # radius of circle
def circle(k,t):
    global freq
    global phi0
    phi_k = phi0[k]+freq[k]*t
    #print("k=%d t=%f freq=%f    phi=%f"%(k,t,freq[k],phi_k))
    #phi_k = phi*(1.0+k0.1)
    return r*np.cos(phi_k), r*np.sin(phi_k)



Time = P.Time
dt = P.dt
num_time = P.num_times
num_phases = P.num_phases
phi0 = np.zeros(num_phases) 
#freq = np.random.normal(1.0,0.00,num_phases)
#freq = np.zeros(num_phases) 
#freq[0:5]=1.0
#freq[5:10]=1.5
#phi0 = np.random.normal(2.0,0.1,num_phases)
#phi0[0] = 0.0
#phi0[1] = 0.5
#phi0[2] = 1.0
#phi0[3] = 2.0


fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False)

# create a figure with an axes

#fig,ax = plt.subplots()

#ax = fig.add_subplot()
# set the axes limits
ax.axis([-1.5,1.5,-1.5,1.5])
# set equal aspect such that the circle is not shown as ellipse
ax.set_aspect("equal")
# create a point in the axes
p = []
for k in range(num_phases):
    point, = ax.plot(1,0, marker="o")
    p.append(point)
for pp in p: 
    pp.set_data([],[])
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
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

    
# Updating function, to be repeatedly called by the animation
def update(n):
    global P
    t = n*P.dt
    # obtain point coordinates 
    str = 't=%f '%t
    print(np.shape(P.data))
    for kk in range(num_phases):
        phi = P.data[n,kk]
        xx = np.cos(phi)
        yy = np.sin(phi)
        x[kk] = xx
        y[kk] = yy
        # set point's coordinates
        str += '   %f    '%(phi)
        p[kk].set_data([x[kk]],[y[kk]])
        #p[kk].set_color(label_color[kk])
    print(str)
    tlabel = 'time:%f'%t
    time_text.set_text(tlabel)
    if t<0.3 or t>0.6:
        time_text.set_color('k')
    else:
        time_text.set_color('r')
    return p[0],p[1],time_text

# create animation with 10ms interval, which is repeated,
# provide the full circle (0,2pi) as parameters
num_loop = 1

anim = FuncAnimation(fig, update, interval=60, blit=True, repeat=False,
                    frames=range(num_loop*num_time))
if FLAG_save==True:
    f = r"./animation.gif" 
    writergif = PillowWriter(fps=30) 
    anim.save(f, writer=writergif)
plt.show()
