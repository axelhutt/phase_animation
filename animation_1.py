#!/usr/bin/env python3


import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


fig = plt.figure(figsize=(8,8))
ax = plt.subplots()
#plt.plot([-1,1],[-1,1],'k')
circle = plt.Circle((0,0),1.0)
ax.add_patch(circle)
phi = 6.28*np.random.random(10)
print(phi)
#index=count()
def animate(i):
    #plt.Circle(( 0.0 , 0.0 ), 1.0,fill = False )
    if i<10:
        #y.append(random.randint(2,20))
        x = np.cos(phi[i])
        y = np.sin(phi[i])
        #plt.style.use("ggplot")  
        print(i,phi[i])
        plt.plot(x,y,marker='o')
    
anim = FuncAnimation(fig, animate, interval=5)
#plt.show()
if FLAG_save==True:
    f = r"./animation.gif" 
    writergif = PillowWriter(fps=30) 
    anim.save(f, writer=writergif)
plt.show()

    