#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

class tfmethods():
    def Wavelettransform_Morlet(self,y):
        from scipy.fft import fft, ifft
        
        def FT_realmotherwavelet(k):
            nu = k*self.df
            return np.sqrt(np.pi)*np.exp(-0.25*(self.sigma-2.0*np.pi*nu)**2)
        
        def FT_complexmotherwavelet(k):
            nu = k*self.df
            return np.sqrt(np.pi)*(np.exp(-0.25*(self.sigma-2.0*np.pi*nu)**2)-np.exp(-0.25*self.sigma**2)*np.exp(-np.pi**2*nu**2))
        
        def FT_motherwavelet(k):
            #return FT_realmotherwavelet(k)
            return FT_complexmotherwavelet(k)
        
        num_time = np.size(y)
        
        self.sigma = 10.0
        center_frequency = self.sigma/(2.0*np.pi)
        self.df    = 1.0/self.T
        
        #self.num_scales = 20 # number of scales
        self.freqs = np.linspace(self.fmin,self.fmax,self.num_scales)
        a = center_frequency/(self.freqs)
        DFT_convolved = np.zeros(num_time,dtype=complex)
        WT = np.zeros((self.num_scales,num_time),dtype=complex)
        
        #print("all prepared......")
        DFT = fft(y)
        #print("size of y:",np.shape(y))
        #print("size of DFT:",np.shape(DFT))
        nums2 = self.num_scales//2
        for l in range(self.num_scales):
            #if l%nums2 == 0:
            #    print("l=%d (%d)"%(l,self.num_scales))
            for k in range(num_time//2):
                DFT_convolved[k] = DFT[k]*FT_motherwavelet(a[l]*k)
                kk = k+num_time//2
                DFT_convolved[kk] = DFT[kk]*FT_motherwavelet(a[l]*(-num_time//2+k))
            WT[l,:] = ifft(DFT_convolved)
            #print("scale %d WT:"%l,WT[l,0:10]," DFT:",DFT[0:10])
        return WT	
    
    
    def Compute_waveletcoherence(self,y1,y2):
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define time step
        dt=self.dt
        
        # now downsampling 
        ntimes=self.num_time
        data=np.zeros((self.num_trials,ntimes,2))
        data[:,:,0]=y1.copy()
        data[:,:,1]=y2.copy()
        del y1
        del y2
        
        timemax=ntimes*dt
        endtime=timemax
        
        self.num_scales = 80
        dff=(self.fmax-self.fmin)/float(self.num_scales)
        
        res01 = self.Wavelettransform_Morlet(data[0,:,0])
        shape = np.shape(res01)
        del res01
        res1 = np.zeros((self.num_trials,shape[0],shape[1]),dtype=complex)
        res2 = np.zeros((self.num_trials,shape[0],shape[1]),dtype=complex)
        cross_wt = np.zeros((shape[0],shape[1]),dtype=complex)
        cross0_wt = np.zeros((shape[0],shape[1]),dtype=complex)
        sum1 = np.zeros((shape[0],shape[1]),dtype=complex)
        sum2 = np.zeros((shape[0],shape[1]),dtype=complex)
        sum01 = np.zeros((shape[0],shape[1]),dtype=complex)
        sum02 = np.zeros((shape[0],shape[1]),dtype=complex)
        cwt  = np.zeros((shape[0],shape[1]))
        print("shape:",shape)
        print("compute statistics....")
        for trial in range(self.num_trials):
            res1[trial,:,:]=self.Wavelettransform_Morlet(data[trial,:,0])
            res2[trial,:,:]=self.Wavelettransform_Morlet(data[trial,:,1])
            print("trial: %d  "%trial)
            #print("trial: %d  "%trial," ",res1[trial,1,1]," ",res2[trial,1,1])
            for i in range(shape[0]):
                for j in range(shape[1]):
                    cross0_wt[i,j] = res1[trial,i,j]*np.conj(res2[trial,i,j])
                    cross_wt[i,j] += cross0_wt[i,j]/float(self.num_trials)
                    sum01[i,j] = res1[trial,i,j]*np.conj(res1[trial,i,j])
                    sum1[i,j] += sum01[i,j]/float(self.num_trials)
                    sum02[i,j] = res2[trial,i,j]*np.conj(res2[trial,i,j])
                    sum2[i,j] += sum02[i,j]/float(self.num_trials)
                    cwt[i,j]  = np.abs(cross_wt[i,j])**2/(sum1[i,j]*sum2[i,j])
            #print("trial:",trial,"  cross:",cross_wt[1,1]," sum1:",sum1[1,1]," sum2:",sum2[1,1])
            #cwt = np.power(np.abs(cross_wt)/(np.abs(sum1)*np.abs(sum2)),2)
            #print("trial:",trial,"  cross:",np.abs(cross0_wt[1:3,400])," sum1:",np.abs(sum01[1:3,400])," sum2:",np.abs(sum02[1:3,400]),"  cwt:",cwt[1:3,400])
            #print("trial:",trial,"  cwt:",cwt[:,420]," ",np.shape(cwt))
        self.waveletcoherence = cwt
        
    
    
    def Compute_phasecoherence(self,y1,y2):
        import numpy.ma as ma
        #import matplotlib.pyplot as plt
        #import scipy.fftpack
        #import matplotlib.mlab as mlab
        #from cwt_modules_Bergner import *
        
        def forceAspect(ax,aspect=1):
            im = ax.get_images()
            extent =  im[0].get_extent()
            ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)
            
        ## define tiem step
        dt=self.dt
        
        # now downsampling 
        counter = 0
        counter_sample = 0
        #print ('downsample data...')
        ntimes=self.num_time
        data=np.zeros((self.num_trials,ntimes,2))
        data[:,:,0]=y1.copy()
        data[:,:,1]=y2.copy()
        del y1
        del y2
        dataavg = np.mean(data,axis=2)
        
        #print ('\n perform time-frequency analysis....')
        
        timemax=ntimes*dt
        endtime=timemax

        #self.fmin=0.5
        #self.fmax=30.0#fs/2.0 # Nyquist
        self.num_scales = 80
        dff=(self.fmax-self.fmin)/float(self.num_scales)
        
        powerthreshold=0.992
        #print("compute wavelet transform.....")
        
        
        res01 = self.Wavelettransform_Morlet(data[0,:,0])
        #res2 = self.Wavelettransform_Morlet(data[0,:,1])
        shape = np.shape(res01)
        del res01
        res1 = np.zeros((self.num_trials,shape[0],shape[1]),dtype=complex)
        res2 = np.zeros((self.num_trials,shape[0],shape[1]),dtype=complex)
        
        for trial in range(self.num_trials):
            res1[trial,:,:]=self.Wavelettransform_Morlet(data[trial,:,0])
            res2[trial,:,:]=self.Wavelettransform_Morlet(data[trial,:,1])
        dangle = np.angle(res1)-np.angle(res2)
        power = np.mean((np.abs(res1)+np.abs(res2))/2,axis=0)
        
        #print("shape of power:",np.shape(power))
        #plt.figure(10)
        #plt.imshow(dangle, interpolation='nearest', cmap=plt.cm.jet,origin='lower', \
        #    extent=[0,self.T,self.fmin,self.fmax])
        #plt.show()
        
        
        #poweravg=np.abs(self.Wavelettransform_Morlet(dataavg))
        maxpower_freq=np.squeeze(np.max(power,axis=0))
        minpower_freq=np.squeeze(np.min(power,axis=0))
        diffpower_freq=maxpower_freq-minpower_freq
        
        ##### plot result
        #print ('compute PLV.....')
        
        PLV_d 		= 0.0
        PLV_d_var 	= 0.0
        f_d_min 		= 0.5
        f_d_max 		= 4.0
        
        PLV_t 		= 0.0
        PLV_t_var 	= 0.0
        f_t_min 		= 4.0
        f_t_max 		= 8.0
        
        PLV_a 		= 0.0
        PLV_a_var 	= 0.0
        f_a_min 		= 8.0
        f_a_max 		= 12.0
        
        PLV_b 		= 0.0
        PLV_b_var 	= 0.0
        f_b_min 		= 12.0
        f_b_max 		= 25.0
        
        PLV_g 		= 0.0
        PLV_g_var 	= 0.0
        f_g_min 		= 25.0
        f_g_max 		= 29.0
        
        num_freqs = np.shape(dangle)[1]
        
        num_f = int((f_d_max-f_d_min)/dff)
        k0 = int(f_d_min/dff)
        #print("d: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_t_max-f_t_min)/dff)
        k0 = int(f_t_min/dff)
        #print("t: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_a_max-f_a_min)/dff)
        k0 = int(f_a_min/dff)
        #print("a: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_b_max-f_b_min)/dff)
        k0 = int(f_b_min/dff)
        #print("b: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        num_f = int((f_g_max-f_g_min)/dff)
        k0 = int(f_g_min/dff)
        #print("g: min=%d num=%d  num_freqs=%d"%(k0,num_f,num_freqs))
        
        windowduration=100
        windowstep = 1 ### !=1: to be implemented
        imax_essential = (ntimes-windowduration)//windowstep
        #print ntimes-windowduration
        PLV = np.zeros((num_freqs,imax_essential))
        helpcos=np.zeros(num_freqs)
        helpsin=np.zeros(num_freqs)
        #print("shape of PLV:",np.shape(PLV)," and helpcos:",np.shape(helpcos))
        
        imax2 = imax_essential//5
        power_m = power.copy()
        #plt.figure(10)
        #plt.imshow(power)
        #plt.show()
        for i in range(imax_essential):
            #if i%imax2 == 0:
            #    print("%d(%d)"%(i,imax_essential))
            ii = i*windowstep
            helpcos[:]=0.0
            helpsin[:]=0.0
            for w in range(windowduration):
                for trial in range(self.num_trials):
                    help = np.squeeze(dangle[trial,:,ii+w])
                    helpcos=helpcos+np.cos(help)/(windowduration*float(self.num_trials))
                    helpsin=helpsin+np.sin(help)/(windowduration*float(self.num_trials))
            help=np.sqrt(helpcos*helpcos+helpsin*helpsin)
            #print("shape of PLV:",np.shape(PLV)," and help:",np.shape(help))
            ### threshold PLV by average power
            if self.FLAG_THR == 1:
                help1=ma.masked_where((np.squeeze(power[:,i])-minpower_freq[i])/diffpower_freq[i]<powerthreshold,help) 
                PLV[:,i]=help1.filled(0)
                help2=ma.masked_where((np.squeeze(power[:,i])-minpower_freq[i])/diffpower_freq[i]<powerthreshold,power[:,i]) 
                power_m[:,i]=help2.filled(0)
                #print(help1)
                #print(help2,'\n')
                #print(help1.filled(0),help2.filled(0))
                #print("thresholding !!!!!!!!!!!!!!!!")
            else:
                PLV[:,i]=help
                power_m[:,i] = power[:,i]
            # delta
            num_f = int((f_d_max-f_d_min)/dff)
            k0 = int(f_d_min/dff)
            for k in range(num_f):
                PLV_d		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_d_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # theta
            num_f = int((f_t_max-f_t_min)/dff)
            k0 = int(f_t_min/dff)
            for k in range(num_f):
                PLV_t		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_t_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # alpha
            num_f = int((f_a_max-f_a_min)/dff)
            k0 = int(f_a_min/dff)
            for k in range(num_f):
                PLV_a		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_a_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # beta
            num_f = int((f_b_max-f_b_min)/dff)
            k0 = int(f_b_min/dff)
            for k in range(num_f):
                PLV_b		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_b_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
            # gamma
            num_f = int((f_g_max-f_g_min)/dff)
            k0 = int(f_g_min/dff)
            for k in range(num_f):
                PLV_g		+=PLV[k0+k,i]/float(num_f*(ntimes-windowduration))
                PLV_g_var	+=PLV[k0+k,i]**2/float(num_f*(ntimes-windowduration))
                
        PLV_d_var = PLV_d_var - PLV_d**2
        PLV_t_var = PLV_t_var - PLV_t**2
        PLV_a_var = PLV_a_var - PLV_a**2
        PLV_b_var = PLV_b_var - PLV_b**2
        PLV_g_var = PLV_g_var - PLV_g**2
        #print("gPLV: delta=%f(%f) theta=%f(%f) alpha=%f(%f) beta=%f(%f) gamma=%f(%f)  number of values=%d"%(\
        #        PLV_d,np.sqrt(PLV_d_var),\
        #        PLV_t,np.sqrt(PLV_t_var),\
        #        PLV_a,np.sqrt(PLV_a_var),\
        #        PLV_b,np.sqrt(PLV_b_var),\
        #        PLV_g,np.sqrt(PLV_g_var),
        #        num_f*(ntimes-windowduration)\
        #    ))
        
        results = []
        results.append([PLV_d,np.sqrt(PLV_d_var),\
            PLV_t,np.sqrt(PLV_t_var),\
            PLV_a,np.sqrt(PLV_a_var),\
            PLV_b,np.sqrt(PLV_b_var),\
            PLV_g,np.sqrt(PLV_g_var),\
            num_f*(ntimes-windowduration)])
        #results.append(poweravg[:,0:ntimes-windowduration])
        results.append(power_m[:,0:ntimes-windowduration])
        results.append(PLV[:,0:ntimes-windowduration])
        
#       fign=plt.figure()
#       ax = fign.add_subplot(211)
#       plt.imshow(poweravg[:,0:ntimes-windowduration], interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-dt*(ntimes-windowduration), self.fmin, self.fmax])
#       plt.colorbar()
#       forceAspect(ax,aspect=2)
#       
#       ax = fign.add_subplot(212)
#       plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, self.fmin, self.fmax])
#       #plt.imshow(PLV, interpolation='nearest', cmap=plt.cm.jet,origin='lower', extent=[0,endtime-windowduration*dt, fmin, fmax])
#       forceAspect(ax,aspect=2)
#       plt.colorbar()
#       
#       plt.show()
#       
        return results
    
        
        
    def Compute_powerspectrum(self,y):
        deltaf = 0.1
        T=10.0/deltaf
        fs = 1.0/self.dt
        num_perseg   = 1*int(T/self.dt)
        
        num_windows =50
        num_onewindow = self.num_time//num_windows 
        fmin_d = 0.1
        fmax_d = 4.0
        fmin_t = 4.0
        fmax_t = 8.0
        fmin_a = 8.0
        fmax_a = 12.0
        fmin_b = 12.0
        fmax_b = 25.0
        power_d     = 0
        power_t     = 0
        power_a     = 0
        power_b     = 0
        power_d_var = 0.0
        power_t_var = 0.0
        power_a_var = 0.0
        power_b_var = 0.0
        
        power_mean = 0.0
        power_var = 0.0
        df_n = fs/float(num_perseg)
        #print("--------- df_n=%f"%df_n)
        for num in range(num_windows):
            y_window = y[num*num_onewindow:(num+1)*num_onewindow-1]
            # Welch PSD
            [power,freqs_n]=mlab.psd(y, NFFT=num_perseg,Fs=fs, detrend=mlab.detrend_none, window=mlab.window_hanning, noverlap=num_perseg*0.9, pad_to=None, sides='onesided', scale_by_freq=True)
            
            fmin_num=int(fmin_d/df_n)
            fmax_num=int(fmax_d/df_n)
            power_d += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_d_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_t/df_n)
            fmax_num=int(fmax_t/df_n)
            power_t += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_t_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_a/df_n)
            fmax_num=int(fmax_a/df_n)
            power_a += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_a_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            fmin_num=int(fmin_b/df_n)
            fmax_num=int(fmax_b/df_n)
            power_b += np.mean(power[fmin_num:fmax_num])/float(num_windows)
            power_b_var += np.mean(power[fmin_num:fmax_num])**2/float(num_windows)
            
            power_mean += power/float(num_windows)
            power_var  += np.power(power,2)/float(num_windows)
            
        power_d_var = power_d_var-power_d**2
        power_t_var = power_t_var-power_t**2
        power_a_var = power_a_var-power_a**2
        power_b_var = power_b_var-power_b**2
        power_var = power_var-np.power(power_mean,2)
        
        print("power_d=%d(%g)  power_t=%d(%g)  power_a=%d(%g)  power_b=%d(%g)  "%(\
                power_d,power_d_var,\
                power_t,power_t_var,\
                power_a,power_a_var,\
                power_b,power_b_var))  
        
        num_freqs_n=np.shape(freqs_n)[0]
        df_n=freqs_n[2]-freqs_n[1]
        
        fmin_n = 0.05
        fmax_n = 60.0
        fmin_u_num_n=int(fmin_n/df_n)
        fmax_u_num_n=int(fmax_n/df_n)
        freqs_p_n  = freqs_n[fmin_u_num_n:fmax_u_num_n]
        p_true     = power_mean[fmin_u_num_n:fmax_u_num_n]
        
        
        results = [] 
        r1 = [power_d,power_d_var,power_t,power_t_var,power_a,power_a_var,power_b,power_b_var]
        results.append(r1)
        results.append(freqs_p_n)
        results.append(p_true)
        return results
    
    
    