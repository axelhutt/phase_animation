#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

import analysis_routines as ar

class simulation(ar.tfmethods):
    
    def __init__(self):
        self.T = 0.9
        self.dt = 0.001
        self.num_time = int(self.T/self.dt)
        
        self.w1 = 2.0*np.pi*40.0
        self.w2 = 2.0*np.pi*30.0
        
        K0 = 300.0
        self.K = np.zeros(self.num_time)
        coupling_start = self.num_time//3
        coupling_end   = 2*self.num_time//3
        self.K[coupling_start:coupling_end]=K0
        
        fR = 0.0/self.T
        self.R = np.zeros(self.num_time)+1.0
        for i in range(self.num_time):
            time = self.dt*float(i)
            self.R[i]=0.5*(np.cos(2.0*np.pi*fR*time)+1.0)*0.7+0.3
        
        self.num_trials = 20
        self.phi = np.zeros((self.num_trials,2,self.num_time))
        self.phi[:,0,0] = 0.0
        self.phi[:,1,0] = np.pi/2.0
        
        self.x = np.zeros((self.num_trials,2,self.num_time))
        
        self.t = 0
        
        self.fmin = 20.0
        self.fmax = 55.0
        self.num_scales = 40
        
    def map_phi_x(self):
        self.x[self.trial,0,self.t] = self.R[self.t]*np.sin(self.phi[self.trial,0,self.t])+0.1*np.random.normal(0.0,1.0)
        self.x[self.trial,1,self.t] = self.R[self.t]*np.sin(self.phi[self.trial,1,self.t])+0.1*np.random.normal(0.0,1.0)
        
    def model(self,phi,dphi):
        
        phi1 = phi[0]
        phi2 = phi[1]
        dphi[0] = self.w1+0.5*self.K[self.t]*np.sin(phi1-phi2)
        dphi[1] = self.w2-0.5*self.K[self.t]*np.sin(phi1-phi2)
        
        return dphi
    
    def integration(self):
        
        dphi = np.zeros(2)
        
        for trial in range(self.num_trials):
            self.trial = trial
            self.phi[trial,0,0] = np.random.random()*2.0*np.pi
            self.phi[trial,1,0] = np.random.random()*2.0*np.pi#np.pi/2.0
            #print("trial: %f %f"%(self.phi[trial,0,0],self.phi[trial,1,0]))
            dphi[:]       = 0.0
            for i in range(self.num_time-1):
                self.t = i
                dphi = self.model(self.phi[trial,:,i],dphi)
                self.phi[trial,:,i+1]=self.phi[trial,:,i]+self.dt*dphi[:]
                self.map_phi_x()
        
        #plt.figure(0)
        #for k in range(self.num_trials):
        #plt.plot(range(self.num_trials),self.phi[:,0,0])
        #plt.show()
        
        
    def compute_tfplot(self):
        import math
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ax = plt.subplot(2,2,1)            
        WT1 = self.Wavelettransform_Morlet(self.x[0,0,:])
        
        absWT = np.abs(WT1)
        plt.imshow(absWT, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        #plt.ylabel("frequency")
        #plt.xlabel("time")
        #plt.title("T-F plot")
        forceAspect(ax,aspect=1)
        
        ax2 = plt.subplot(2,2,2)            
        WT2 = self.Wavelettransform_Morlet(self.x[0,1,:])
        absWT2 = np.abs(WT2)
        plt.imshow(absWT2, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        #plt.ylabel("frequency")
        #plt.xlabel("time")
        #plt.title("T-F plot")
        forceAspect(ax2,aspect=1)
        
        ax = plt.subplot(2,2,3)            
        phaseWT1 = np.arctan(np.imag(WT1)/np.real(WT1))
        plt.imshow(phaseWT1, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        #plt.ylabel("frequency")
        #plt.xlabel("time")
        #plt.title("T-F plot")
        forceAspect(ax,aspect=1)
        
        ax2 = plt.subplot(2,2,4)
        phaseWT2 = np.arctan(np.imag(WT2)/np.real(WT2))
        plt.imshow(phaseWT2, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        #plt.ylabel("frequency")
        #plt.xlabel("time")
        #plt.title("T-F plot")
        forceAspect(ax2,aspect=1)
    
        fig0 = plt.figure(4)
        ax = plt.subplot(1,1,1)            
        dphase = phaseWT1-phaseWT2
        
        plt.imshow(dphase, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        forceAspect(ax,aspect=1)

        plt.show()
    
    def compute_waveletcoherence(self):
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        self.Compute_waveletcoherence(self.x[:,0,:],self.x[:,1,:])
        print("plot coherence...")
        fig = plt.figure(0)
        ax2 = plt.subplot(1,1,1)            
        plt.imshow(self.waveletcoherence, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        plt.ylabel("frequency")
        plt.xlabel("time")
        plt.title("wavelet coherence")
        forceAspect(ax2,aspect=1)
        plt.show()
    
    def compute_phasecoherence(self):
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        self.FLAG_THR = 1
        results = self.Compute_phasecoherence(self.x[0,0,:],self.x[0,1,:])
        
        fig = plt.figure(1)
        ax2 = plt.subplot(1,1,1)            
        plt.imshow(results[2], interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        plt.ylabel("frequency")
        plt.xlabel("time")
        plt.title("phase coherence plot")
        forceAspect(ax2,aspect=1)

        fig = plt.figure(2)
        ax2 = plt.subplot(1,1,1)            
        plt.imshow(results[1], interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        plt.ylabel("frequency")
        plt.xlabel("time")
        plt.title("power plot")
        forceAspect(ax2,aspect=1)
        
        plt.show()
        
    def compute_spectralcoherence(self):
        import scipy.signal as ss
        
        fs      = 1.0/self.dt
        df      = 1.0#(self.fmax-self.fmin)/float(self.num_scales)
        T       = int(1.0/(df*self.dt))
        if T>self.num_time:
            T=self.num_time
        overlap = int(T*0.8)
        print("fs=%f nperseg=%d noverlap=%d"%(fs,T,overlap))
        self.freq,self.SC = ss.coherence(self.x[0,0,:],self.x[0,1,:],fs=fs,nperseg=T,noverlap=overlap)
        print("shape of SC:",np.shape(self.SC))
        fig = plt.figure(3)
        ax2 = plt.subplot(1,1,1)            
        plt.imshow(self.SC, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
                    extent=[0,self.T,self.fmin,self.fmax])
        plt.colorbar()
        plt.ylabel("frequency")
        plt.xlabel("time")
        plt.title("spectral coherence")
        forceAspect(ax2,aspect=1)
        
        plt.show()
        
    
    def plot_ts(self):
        
        times = np.linspace(0,self.T,self.num_time)
        fig = plt.figure(1)
        ax2 = plt.subplot(3,1,1)
        plt.plot(times,self.K,'k')
        ax2 = plt.subplot(3,1,2)
        plt.plot(times,self.x[0,0,:],'k',times,self.x[0,1,:],'r')
        ax2 = plt.subplot(3,1,3)
        plt.plot(times,self.phi[0,0,:]-self.phi[0,1,:],'k')
        
        plt.show()
        
    def write_out(self):
        label = './unknown_timeseries_new.dat'
        f = open(label,'w+')
        times = np.linspace(0,self.T,self.num_time)
        for i in range(self.num_time):
            str = '%f %f    %f\n'%(times[i],self.x[0,0,i],self.x[0,1,i])
            f.write(str)
        f.close()
        print("written %s"%label)

        label = './unknown_timeseries_phases.dat'
        f = open(label,'w+')
        times = np.linspace(0,self.T,self.num_time)
        
        for i in range(self.num_time):
            str = '%f  '%times[i]
            for k in range(self.num_trials):
                str += '  %f %f'%(self.phi[k,0,i],self.phi[k,1,i])
            str += '\n'
            f.write(str)
        f.close()
        print("written %s"%label)
        
S = simulation()
S.integration()
#S.compute_spectralcoherence()
S.plot_ts()
#S.compute_tfplot()
#S.compute_phasecoherence()
#S.compute_waveletcoherence()
S.write_out()
