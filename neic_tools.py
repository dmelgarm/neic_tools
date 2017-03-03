
class neic_catalog:
    '''
    Read an manipulate the USGS NEIC finite fault database
    '''
    
    def __init__(self,path_to_files,catalog_file,percent_cutoff=0.1,quiet=True):
        '''
        Initalize the class
        '''
        
        from numpy import genfromtxt
        
        self.percent_cutoff=percent_cutoff
        self.path_to_files=path_to_files
        self.catalog_file=catalog_file
        self.Nevents=len(genfromtxt(catalog_file))
        self.get_events_metadata()
        self.get_all_mags_moments()
        self.get_all_rise_times(quiet=quiet)
        self.get_all_rakes()
        self.define_event_type()
        
    
    
    def get_events_metadata(self):
        ''' 
        read ID lat lon,for all events
        '''
        
        from numpy import genfromtxt
        
        self.IDs=genfromtxt(self.catalog_file,usecols=6,dtype='S')
        self.origin_times=genfromtxt(self.catalog_file,usecols=0,dtype='S')
        self.epicenters=genfromtxt(self.catalog_file,usecols=[1,2,3])
        

    def get_all_mags_moments(self,quiet=True):
        '''
        Get all mangitudes and moments from the slip inversions
        '''
        
        from numpy import zeros,log10
        
        
        self.moments=zeros(self.Nevents)
        for k in range(self.Nevents):
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            self.moments[k]=self.get_moment(fault)
        
        self.magnitudes=(2./3)*(log10(self.moments)-9.1)


    def get_all_rise_times(self,quiet=True):
        '''
        Get all mean rise times
        '''
        
        from numpy import zeros
        
        self.mean_rise_times=zeros(self.Nevents)
        self.max_rise_times=zeros(self.Nevents)
        for k in range(self.Nevents):
            #read fault info
            fault_file=self.path_to_files+self.IDs[k]+'.param'
            if quiet==False:
                print 'Reading '+fault_file
            fault=self.read_neic_param(fault_file)
            self.mean_rise_times[k]=self.get_mean_rise_time(fault)
            self.max_rise_times[k]=self.get_max_rise_time(fault)
            
            
            
            
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
    
        