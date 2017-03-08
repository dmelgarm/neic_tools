# -*- coding: utf-8 -*-

class neic_catalog:
    '''
    Read an manipulate the USGS NEIC finite fault database
    '''
    
    def __init__(self,path_to_files,catalog_file,percent_cutoff=0.1,get_stfs=True,
            quiet=True,stf_dt=0.1,stf_tmax=250,moment_percent=0.95,percent_cutoff_vrup=0.3):
        '''
        Initalize the class
        '''
        
        from numpy import genfromtxt
        
        self.percent_cutoff=percent_cutoff
        self.path_to_files=path_to_files
        self.catalog_file=catalog_file
        self.percent_cutoff_vrup=percent_cutoff_vrup
        self.Nevents=len(genfromtxt(catalog_file))
        self.get_events_metadata()
        self.get_all_mags_moments()
        self.get_all_slip_and_rise_times(quiet=quiet)
        self.get_all_rakes()
        self.define_event_type()
        if get_stfs==True:
            self.get_all_stfs(stf_dt,stf_tmax)
            self.get_event_durations(moment_percent=moment_percent)
        
    
    
    def get_events_metadata(self):
        ''' 
        read ID lat lon,for all events
        '''
        
        from numpy import genfromtxt
        
        self.IDs=genfromtxt(self.catalog_file,usecols=6,dtype='S')
        self.origin_times=genfromtxt(self.catalog_file,usecols=0,dtype='S')
        self.epicenters=genfromtxt(self.catalog_file,usecols=[1,2,3])
        self.event_class=genfromtxt(self.catalog_file,usecols=10,dtype='S')
        self.names=genfromtxt(self.catalog_file,usecols=12,dtype='S')

    def get_all_mags_moments(self,quiet=True):
        '''
        Get all mangitudes and moments from the slip inversions
        '''
        
        from numpy import zeros,log10
        
        
        self.moments=zeros(self.Nevents)
        self.raw_moments=[]
        for k in range(self.Nevents):
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            self.moments[k]=self.get_moment(fault)
            self.raw_moments.append(fault[:,10]/1.e7)
        
        self.magnitudes=(2./3)*(log10(self.moments)-9.1)


    def get_all_slip_and_rise_times(self,quiet=True):
        '''
        Get all mean rise times
        '''
        
        from numpy import zeros,where
        
        self.mean_rise_times=zeros(self.Nevents)
        self.max_rise_times=zeros(self.Nevents)
        self.mean_slip=zeros(self.Nevents)
        self.max_slip=zeros(self.Nevents)
        self.time_start=[]
        self.time_up=[]
        self.time_down=[]
        self.subfault_distance_to_hypocenter=[]
        self.subfault_rupture_velocity=[]
        self.mean_rupture_velocity=zeros(self.Nevents)
        self.std_rupture_velocity=zeros(self.Nevents)
        self.mean_pulse_lengths=zeros(self.Nevents)
        
        for k in range(self.Nevents):
            
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            self.mean_rise_times[k]=self.get_mean_rise_time(fault)
            self.max_rise_times[k]=self.get_max_rise_time(fault)
            self.mean_slip[k]=self.get_mean_slip(fault)
            self.max_slip[k]=self.get_max_slip(fault)
            
            #time stuff
            time_start,time_up,time_down=self.get_asymetric_times(fault)
            self.time_start.append(time_start)
            self.time_up.append(time_up)
            self.time_down.append(time_down)
            
            #distance stuff
            distance=self.get_distance_to_hypo(fault)
            self.subfault_distance_to_hypocenter.append(distance)
            
            #get the implied rupture velocity toe ach subfault
            rupture_velocity=distance/time_start
            i=where(time_start==0)[0]
            rupture_velocity[i]=0
            self.subfault_rupture_velocity.append(rupture_velocity)
            
            #And the mean velocity for that model
            vr,std_vr=self.get_mean_rupture_velocity(fault,rupture_velocity,distance)
            self.mean_rupture_velocity[k]=vr
            self.std_rupture_velocity[k]=std_vr
            
            #Get the estimated pulse width
            pulse_length=self.get_mean_pulse_length(fault,vr)
            self.mean_pulse_lengths[k]=pulse_length
            
            
            
    def get_all_stfs(self,stf_dt,stf_tmax):
        
        ''' 
        Get source time functions for all the events
        '''

        self.stf_times=[]
        self.stf_moment_rates=[]
        for k in range(self.Nevents):
            moment=self.raw_moments[k]
            tstart=self.time_start[k]
            tup=self.time_up[k]
            tdown=self.time_down[k]
            t,stf=self.get_source_time_function(moment,tstart,tup,tdown,stf_dt,stf_tmax) 
            self.stf_times.append(t)
            self.stf_moment_rates.append(stf) 
            
    
    def get_event_durations(self,moment_percent=0.95):
        '''
        get event duration defiend as the time when moment_percent of total moment is achieved
        '''
        from numpy import zeros,where
        from scipy.integrate import cumtrapz,trapz
        
        self.event_durations=zeros(self.Nevents)
        self.centroid_times=zeros(self.Nevents)
        for k in range(self.Nevents):
            
            #Calcualte moment as a function of time and total moment
            moment=cumtrapz(self.stf_moment_rates[k],self.stf_times[k],initial=0)
            Mo=trapz(self.stf_moment_rates[k],self.stf_times[k])
            
            #When is threshodl exceeded?
            i=where(moment>moment_percent*Mo)[0]
            duration=self.stf_times[k][i[0]]
            centroid_time=(self.stf_moment_rates[k]*self.stf_times[k]).sum()/self.stf_moment_rates[k].sum()
            self.event_durations[k]=duration
            self.centroid_times[k]=centroid_time
            
            
    def get_distance_to_hypo(self,fault):
        '''
        Get straight line distance from subfault center to hypocenter
        '''   
        
        from numpy import argmin,size,ones
        from pyproj import Geod
        
        #Projection object for distances
        p=Geod(ellps='WGS84')
        
        #Firs find the coordiantes of the hypocenter
        time_start=fault[:,7]
        i=argmin(time_start)
        if size(i)>1: #Hypocenter is ambiguos
            'ERROR: Too many possible hypocenters'
        else:
            hypo=fault[i,0:3]
            # get horizontal distances
            az,baz,dist_h=p.inv(ones(len(fault))*hypo[1],ones(len(fault))*hypo[0],fault[:,1],fault[:,0])
            dist_h=dist_h/1000.
            #vertical dsitances
            dist_v=abs(ones(len(fault))*hypo[2]-fault[:,2])
            #total distance
            distance=(dist_h**2+dist_v**2)**0.5
            
        return distance
            
            
    def get_all_rakes(self,quiet=True):
        '''
        Get all mean rake directions
        '''
        
        from numpy import zeros
        
        self.mean_rakes=zeros(self.Nevents)
        for k in range(self.Nevents):
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            self.mean_rakes[k]=self.get_mean_rake(fault)
            

    def define_event_type(self):
        '''
        Based on the rake call the event, thrust, normal or strike slip
        '''
        
        from numpy import zeros
        
        self.event_type=zeros(len(self.magnitudes)).astype('string')
        for k in range(len(self.magnitudes)):    
            rake=self.mean_rakes[k]
            if rake>45 and rake<135:
                self.event_type[k]='T'
            elif rake >-135 and rake <-45:
                self.event_type[k]='N'
            else:
                self.event_type[k]='S'

                
                                
    def run_regression(self,fix_exponent=True,Niter=100e3,burn=20e3,prior='normal',
            select_event_type=None,inversion_type='linear',dependent_variable='rise_time'):
        '''
        Run regression for new scaling coefficients, functional form is:
        
            trise = A*Mo^k
        '''
        
        from numpy import log10,ones,log,where
        from numpy.linalg import lstsq
        import pymc as pm   
        
        if select_event_type==None:
            rise_time=self.mean_rise_times
            durations=self.event_durations
            centroid_times=self.centroid_times
            pulse_lengths=self.mean_pulse_lengths
            moment=self.moments        
        elif select_event_type=='i': #megathrust events
            i=where(self.event_class==select_event_type)[0]
            rise_time=self.mean_rise_times[i]
            durations=self.event_durations[i]
            centroid_times=self.centroid_times[i]
            pulse_lengths=self.mean_pulse_lengths[i]
            moment=self.moments[i]
        elif select_event_type=='u': #uper plate events
            i=where(self.event_class==select_event_type)[0]
            rise_time=self.mean_rise_times[i]
            durations=self.event_durations[i]
            centroid_times=self.centroid_times[i]
            pulse_lengths=self.mean_pulse_lengths[i]
            moment=self.moments[i]
        elif select_event_type=='l': #mantle events
            i=where(self.event_class==select_event_type)[0]
            rise_time=self.mean_rise_times[i]
            durations=self.event_durations[i]
            centroid_times=self.centroid_times[i]
            pulse_lengths=self.mean_pulse_lengths[i]
            moment=self.moments[i]
        elif select_event_type=='n/a': #non subduction
            i=where(self.event_class==select_event_type)[0]
            rise_time=self.mean_rise_times[i]
            durations=self.event_durations[i]
            centroid_times=self.centroid_times[i]
            pulse_lengths=self.mean_pulse_lengths[i]
            moment=self.moments[i]
        else:
            print 'ERROR: unknown event class'
            return
        
        #How many left over?
        Nevents=len(moment)
        
        if inversion_type=='linear':
            #G matrix and d vector   
            if fix_exponent==True: #force k=1/3 in regression
                if dependent_variable=='rise_time':
                    d=log10(rise_time)-(1./3)*log10(moment)
                elif dependent_variable=='duration':
                    d=log10(durations)-(1./3)*log10(moment)
                elif dependent_variable=='centroid_time':
                    d=log10(centroid_times)-(1./3)*log10(moment)
                elif dependent_variable=='pulse_length':
                    d=log10(pulse_lengths)-(1./3)*log10(moment)
                G=ones((Nevents,1))
                
                #Run regression          
                coefficients,a,b,c=lstsq(G,d)
                
                ##Determine A and k
                A=10**coefficients
                k=1./3
            else:
                if dependent_variable=='rise_time':
                    d=log10(rise_time)
                elif dependent_variable=='duration':
                    d=log10(durations)
                elif dependent_variable=='centroid_time':
                    d=log10(centroid_times)
                elif dependent_variable=='pulse_length':
                    d=log10(pulse_lengths)
                                                    
                G=ones((Nevents,2))
                G[:,1]=log10(moment)
                
                #Run regression          
                coefficients,a,b,c=lstsq(G,d)
                
                ##Determine A and k
                A=10**coefficients[0]
                k=coefficients[1]
                
            return A,k  
        
        elif inversion_type=='bayesian':
            
            if prior=='normal':
                A = pm.distributions.Normal('A', mu = -5, tau = 1)
                k = pm.distributions.Normal('k', mu = 1./3, tau = 0.15)

            elif prior=='uninformative':            
                A = pm.Uniform('A', -10, 0)
    
                @pm.stochastic(observed=False)
                def k(value=0):
                    return -1.5 * log(1 + value ** 2)
            
            @pm.stochastic(observed=False)
            def sigma(value=1):
                return -log(abs(value))
            
            # Define the form of the model and likelihood
            @pm.deterministic
            def y_model(x=log10(moment), A=A, k=k):
                return A + k * x
            
            if dependent_variable=='rise_time':
                y = pm.Normal('y', mu=y_model, tau=1. / sigma ** 2, observed=True, value=log10(rise_time))
            elif dependent_variable=='duration':
                y = pm.Normal('y', mu=y_model, tau=1. / sigma ** 2, observed=True, value=log10(durations))
            elif dependent_variable=='centroid_time':
                y = pm.Normal('y', mu=y_model, tau=1. / sigma ** 2, observed=True, value=log10(centroid_times))
            elif dependent_variable=='pulse_length':
                y = pm.Normal('y', mu=y_model, tau=1. / sigma ** 2, observed=True, value=log10(pulse_lengths))
            else:
                print 'ERROR: Unknown dependent variable type'
                return
                
            # package the full model in a dictionary
            model = dict(A=A, k=k, sigma=sigma,
                        y_model=y_model, y=y)

            # run the basic MCMC
            mcmc = pm.MCMC(model)
            mcmc.sample(iter=Niter, burn=burn)
            
            #Unwrap coeficients
            A=10**(mcmc.stats('A')['A']['mean'])
            k=mcmc.stats('k')['k']['mean']
            
            return A,k,mcmc
            
            
            

    def get_source_time_function(self,moment,tstart,tup,tdown,dt,tmax):
        '''
        For an array of rupture start times, time ups and time down build entire STF
        '''
        
        from numpy import cos,pi,arange,where,zeros,roll
        from scipy.integrate import trapz
        
        t=arange(0,tmax,dt)
        stf=zeros(len(t))
        for k in range(len(moment)):
            
            #Check that nothing is zero
            if tup[k]==0 or tdown[k]==0:
                s=zeros(len(t))
            else:
                s1=(1./(tup[k]+tdown[k]))*(1-cos((pi*t)/tup[k]))
                i=where(t>tup[k])[0]
                s1[i]=0
                s2=(1./(tup[k]+tdown[k]))*(1+cos((pi*(t-tup[k]))/tdown[k]))
                i=where(t<=tup[k])[0]
                s2[i]=0 
                i=where(t>tup[k]+tdown[k])[0]
                s2[i]=0
                #add the two 
                s=s1+s2
                #shift to right onset time
                dN=int(tstart[k]/dt)
                s=roll(s,dN)
                #rescale to the correct moment
                area=trapz(s,t)
                scale=moment[k]/area
                s=s*scale
            stf+=s
        
        #Sanity check
        Mo=trapz(stf,t)
        if (Mo/moment.sum()>1.02) or (Mo/moment.sum()<0.98):
            print 'ERROR moment discrepancy is larger than 2%'
        else:
            return t,stf
                        

    def get_mean_rupture_velocity(self,fault,rupture_velocity,distance):
        '''
        Filter by percent cutoff slip and get WEIGHTED mean rupture velocity
        '''
        from numpy import where,average,sum
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        slip=fault[:,3]
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff_vrup*peak_slip)[0]
        distance=distance[i]
        rupture_velocity=rupture_velocity[i]
        slip=slip[i]
        
        #Mask out subfaults too close to hypocenter
        #i=where(distance<self.mask_distance)[0]
        #rupture_velocity=rupture_velocity[i]
        #slip=slip[i]
        
        #Eliminate vrup=0
        i=where(rupture_velocity>0)[0]
        distance=distance[i]
        rupture_velocity=rupture_velocity[i]
        slip=slip[i]
        
        #mean_rupture_velocity=rupture_velocity.mean()
        #std_rupture_velocity=rupture_velocity.std()
        
        weights=distance
        mean_rupture_velocity=average(rupture_velocity,weights=weights)
        #And get weighted std
        N=len(weights)
        A=sum(weights*(rupture_velocity-rupture_velocity.mean())**2)
        B=((N-1.)/N)*weights.sum()
        std_rupture_velocity=(A/B)**0.5
        
        return mean_rupture_velocity,std_rupture_velocity                                                
                                                                        

    def get_mean_pulse_length(self,fault,vr_mean):
        '''
        Filter by percent cutoff slip and get WEIGHTED mean rupture velocity
        '''
        from numpy import where,average,sum
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        rise_times=fault[i,8]+fault[i,9]
        pulse_lengths=rise_times*vr_mean
        
        return pulse_lengths.mean()
        

        
        return mean_rupture_velocity,std_rupture_velocity,distance      
          
                    
                                                                                                                                                                                                                                   
    def get_all_rise_time_pdfs(self,bins=10,quiet=True):
        '''
        Get all pdfs for rise times
        '''
    
        from numpy import zeros
        
        self.pdf_frequency=zeros((self.Nevents,bins))
        self.pdf_bin_edges=zeros((self.Nevents,bins+1))
        
        for k in range(self.Nevents):
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            h,b=self.get_rise_time_pdf(fault,bins=bins)
            self.pdf_frequency[k,:]=h
            self.pdf_bin_edges[k,:]=b     
        

    def read_neic_param(self,fault_file):
        '''
        Parse the param text format
        '''
        
        
        from numpy import r_,array,expand_dims
        
        f=open(fault_file)
        #get number of segments
        line=f.readline()
        n_seg=int(line.split()[-1])
        segments_read=0
        
        fault_out=array([])
        kread=0
        while True:
            #read fault segment header
            line=f.readline()
            if line=='':
                break
            Dx=float(line.split()[6].replace('km',''))
            Dy=float(line.split()[10].replace('km',''))
            Nsubfaults=int(line.split()[4])*int(line.split()[8])
            #skip junk lines
            for k in range(8):
                line=f.readline()
            #read segment data
            for k in range(Nsubfaults):
                line=f.readline()
                if kread==0:
                    fault_out=expand_dims(array(line.split()).astype('float'),0)
                else:
                    fault_out=r_[fault_out,expand_dims(array(line.split()).astype('float'),0)]   
                kread+=1
            #Done
            if line=='':
                break
        return fault_out
        
        
    def get_moment(self,fault):
        '''
        Get moment directly from param file
        '''
        
        moment=fault[:,10].sum()/1.e7
        
        return moment
        
        
        
    def get_mean_rake(self,fault):
        '''
        Get mean rake for a fault model
        '''
        
        from numpy import where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        
        #get mean rise time of those faults
        j=where(fault[:,4]>180)
        fault[j,4]=fault[j,4]-360
        mean_rake=fault[i,4].mean()
        
        return mean_rake
                
        
    def get_mean_slip(self,fault):
        '''
        Get mean rise time for a fault model
        '''
        
        from numpy import where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        
        #get mean rise time of those faults
        mean_slip=(fault[i,3]).mean()/100.
        
        return mean_slip
        
    def get_max_slip(self,fault):
        '''
        Get mean rise time for a fault model
        '''
        
        from numpy import where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        
        #get mean rise time of those faults
        max_slip=(fault[i,3]).max()/100.
        
        return max_slip
        
    def get_asymetric_times(self,fault):
        '''
        Get up and down times of the assymetric cosine function
        '''
        time_start=fault[:,7]
        time_up=fault[:,8]
        time_down=fault[:,9]
        return time_start,time_up,time_down
        
                    
    def get_mean_rise_time(self,fault):
        '''
        Get mean rise time for a fault model
        '''
        
        from numpy import where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        
        #get mean rise time of those faults
        mean_rise_time=(fault[i,8]+fault[i,9]).mean()
        
        return mean_rise_time


    def get_max_rise_time(self,fault):
        '''
        Get max rise time for a fault model
        '''
        
        from numpy import where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
        
        #get mean rise time of those faults
        max_rise_time=(fault[i,8]+fault[i,9]).max()
        
        return max_rise_time        
                        
    
    def get_rise_time_pdf(self,fault,bins=10):
        '''
        get the pdf of rise times for an event
        '''
        
        from numpy import histogram,where
        
        #Find peak slip
        peak_slip=fault[:,3].max()
        
        #Find subfaults alrger than percent_cutoff of peak slip
        i=where(fault[:,3]>self.percent_cutoff*peak_slip)[0]
         
        hist,bin_edges=histogram((fault[i,8]+fault[i,9]),bins=bins)
        
        return hist,bin_edges
    
    
    

    
        