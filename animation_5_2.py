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
        self.num_phases = (np.shape(data0)[1]-1)//2
        self.Time = data0[-1:0]
        self.data = data0[:,1:2*self.num_phases+1]
        self.dt = data0[1,0]-data0[0,0]
     
    def simulate(self):
        self.Time = 10.0
        self.dt = 0.1
        self.num_times = int(self.Time/self.dt)
        self.num_phases = 1
        self.data = np.zeros((self.num_times,self.num_phases*2))
        self.data_R = np.zeros((self.num_times,self.num_phases*2))
        self.time = np.zeros(self.num_times)
        self.x = np.zeros((self.num_times,self.num_phases*2))
        freq1 = 1.0
        freq2 = 1.0
        self.phi0 = np.zeros(2)
        self.phi0[0]=0.5
        self.phi0[1]=1.5
        for i in range(self.num_times):
            t = i*self.dt
            self.time[i]=t
            for k in range(self.num_phases):
                R1 = np.random.normal(1.0,0.1)
                phi1 = self.phi0[k*2]+freq1*t + np.random.normal(0.0,0.1)
                R2 = np.random.normal(1.0,0.1)
                phi2 = self.phi0[k*2+1]+freq2*t + np.random.normal(0.0,0.1)
                self.data[i,k*2]=phi1
                self.data[i,k*2+1]=phi2
                self.data_R[i,k*2]=R1
                self.data_R[i,k*2+1]=R2
                self.x[i,0]=R1*np.cos(phi1)
                self.x[i,1]=R1*np.sin(phi1)
                
    def plot_ts(self):
        fig = plt.figure()
        ax = fig.add_subplot(211)
        ax.plot(self.time,self.x[:,0])
        ax = fig.add_subplot(212)
        ax.plot(self.time,self.x[:,1])
        plt.show()
        
        ndata = np.zeros(self.num_times)
        ndata1 = np.cos(1.0*self.time)+np.random.normal(0.0,1.0,self.num_times)
        ndata2 = np.cos(1.0*self.time)+np.random.normal(0.0,1.0,self.num_times)
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
        ax.plot(self.time,ndata1)
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.plot(self.time,ndata2)
        plt.show()

FLAG_dphase = 1

P = phases()
#P.read_in()
P.simulate()
P.plot_ts()

print("num_phases: %d"%(P.num_phases))

Time = P.Time
dt = P.dt
num_time = P.num_times
if FLAG_dphase==1:
    num_phases = P.num_phases
else:
    num_phases = P.num_phases*2
phi0 = np.zeros(num_phases) 

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False)

ax.axis([-1.5,1.5,-1.5,1.5])
ax.set_aspect("equal")
p = []
if FLAG_dphase==1:
    for k in range(num_phases):
        point, = ax.plot(1,0, marker="o")
        p.append(point)
else:
    for k in range(num_phases*2):
        point, = ax.plot(1,0, marker="o")
        p.append(point)
for pp in p: 
    pp.set_data([],[])
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

num = 100
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
def update_dphase(n):
    global P
    t = n*P.dt
    # obtain point coordinates 
    str = 't=%f '%t
    #print(np.shape(P.data))
    for kk in range(num_phases):
        phi = P.data[n,2*kk+1]-P.data[n,2*kk]
        xx = np.cos(phi)
        yy = np.sin(phi)
        x[kk] = xx
        y[kk] = yy
        # set point's coordinates
        str += '   %f    '%(phi)
        p[kk].set_data([x[kk]],[y[kk]])
        #p[kk].set_color(label_color[kk])
        #print(str)
    tlabel = 'time:%f'%t
    time_text.set_text(tlabel)
    if t<0.3 or t>0.6:
        time_text.set_color('k')
    else:
        time_text.set_color('r')
    if num_phases == 1:
        return p[0],time_text
    if num_phases == 2:
        return p[0],p[1],time_text
    else:
        return p[0],p[1],\
                p[2],p[3],\
                p[4],p[5],\
                p[6],p[7],\
                p[8],p[9],\
                p[10],p[11],\
                p[12],p[13],\
                p[14],p[15],\
                p[16],p[17],\
                p[18],p[19],time_text

def update_phase(n):
    global P
    t = n*P.dt
    # obtain point coordinates 
    str = 't=%f '%t
    #print(np.shape(P.data))
    for kk in range(num_phases):
        phi = P.data[n,kk]
        #xx = P.data_R[n,kk]*np.cos(phi)
        #yy = P.data_R[n,kk]*np.sin(phi)
        xx = np.cos(phi)
        yy = np.sin(phi)
        x[kk] = xx
        y[kk] = yy
        # set point's coordinates
        str += '   %f    '%(phi)
        p[kk].set_data([x[kk]],[y[kk]])
        #p[kk].set_color(label_color[kk])
        #print(str)
    tlabel = 'time:%f'%t
    time_text.set_text(tlabel)
    if t<0.3 or t>0.6:
        time_text.set_color('k')
    else:
        time_text.set_color('r')
    if num_phases == 2:
        return p[0],p[1],time_text
    else:
        return p[0],p[1],\
                p[2],p[3],\
                p[4],p[5],\
                p[6],p[7],\
                p[8],p[9],\
                p[10],p[11],\
                p[12],p[13],\
                p[14],p[15],\
                p[16],p[17],\
                p[18],p[19],time_text
    
# create animation with 10ms interval, which is repeated,
# provide the full circle (0,2pi) as parameters
num_loop = 1

if FLAG_dphase == 1:
    anim = FuncAnimation(fig, update_dphase, interval=60, blit=True, repeat=False,
                    frames=range(num_loop*num_time))
else:
    anim = FuncAnimation(fig, update_phase, interval=60, blit=True, repeat=False,
                    frames=range(num_loop*num_time))

if FLAG_save==True:
    f = r"./animation.gif" 
    writergif = PillowWriter(fps=30) 
    anim.save(f, writer=writergif)
plt.show()
