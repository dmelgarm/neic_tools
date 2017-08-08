from numpy import genfromtxt
from mudpy import gmttools
from string import rjust
from mudpy import gmttools

home='/Users/dmelgar/fakequakes/'
project_name='tohoku_usgs'
run_name='tohoku'
for k in range(20):
    num=rjust(str(k),6,'0')
    f=genfromtxt(home+project_name+'/output/ruptures/'+run_name+'.'+num+'.rupt')
    
    fout=open(home+project_name+'/output/ruptures/'+run_name+'.'+num+'.usgs','w')
    fout.write('#Total number of fault_segments=     3\n')
    fout.write('#Fault_segment =   1 nx(Along-strike)=  25 Dx= 25.00km ny(downdip)=   5 Dy= 16.60km\n')
    fout.write('#Boundary of Fault_segment     1. EQ in cell 12,3. Lon: 142.3700   Lat: 38.3200\n')
    fout.write('#Lon.  Lat.  Depth\n')
    fout.write('143.80890       40.66050       13.96590\n')
    fout.write('141.60339       35.33500       13.96590\n')
    fout.write('140.75380       35.55160       34.83410\n')
    fout.write('142.95930       40.87710       34.83410\n')
    fout.write('143.80890       40.66050       13.96590\n')
    fout.write('#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo\n')
    for i in range(0,125):
        slip=f[i,9]
        mo=slip*f[i,13]*25e3*16e3*1e7
        slip=slip*100
        line='%.5f\t%.5f\t%.2f\t%.2f\t90.00\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\n' % (f[i,2],f[i,1],f[i,3],slip,f[i,4],f[i,5],f[i,12],f[i,7]/2.,f[i,7]/2.,mo)
        fout.write(line)  
        
    fout.write('Fault_segment =   2 nx(Along-strike)=  25 Dx= 25.00km ny(downdip)=   5 Dy= 16.60km\n')
    fout.write('Boundary of Fault_segment     2. EQ in cell 12,3. Lon: 142.3700   Lat: 38.3200\n')
    fout.write('#Lon.  Lat.  Depth\n')
    fout.write('144.73550       40.43200        2.27270\n')
    fout.write('142.45410       35.10650        2.27270\n')
    fout.write('141.55321       35.32850       13.49400\n')
    fout.write('143.83450       40.65400       13.49400\n')
    fout.write('144.73550       40.43200        2.27270\n')
    fout.write('#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo\n')
    
    for i in range(125,250):
        slip=f[i,9]
        mo=slip*f[i,13]*25e3*16e3*1e7
        slip=slip*100
        line='%.5f\t%.5f\t%.2f\t%.2f\t90.00\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\n' % (f[i,2],f[i,1],f[i,3],slip,f[i,4],f[i,5],f[i,12],f[i,7]/2.,f[i,7]/2.,mo)
        fout.write(line)  

    fout.write('#Fault_segment =   3 nx(Along-strike)=  25 Dx= 25.00km ny(downdip)=   3 Dy= 16.60km\n')
    fout.write('#Boundary of Fault_segment     3. EQ in cell 12,3. Lon: 142.3700   Lat: 38.3200\n')
    fout.write('#Lon.  Lat.  Depth\n')
    fout.write('142.93370       40.88180       37.21840\n')
    fout.write('140.64619       35.55630       37.21840\n')
    fout.write('140.14529       35.67950       54.21530\n')
    fout.write('142.43269       41.00500       54.21530\n')
    fout.write('142.93370       40.88180       37.21840\n')
    fout.write('#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo\n')       

    for i in range(250,325):
        slip=f[i,9]
        mo=slip*f[i,13]*25e3*16e3*1e7
        slip=slip*100
        line='%.5f\t%.5f\t%.2f\t%.2f\t90.00\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\n' % (f[i,2],f[i,1],f[i,3],slip,f[i,4],f[i,5],f[i,12],f[i,7]/2.,f[i,7]/2.,mo)
        fout.write(line)         
        
        
fout.close()