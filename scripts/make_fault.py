import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from numpy import array,c_,arange,ones,savetxt

fout='/Users/dmelgar/USGSFF/tohoku_usgs.fault'
path_to_files='/Users/dmelgar/USGSFF/param/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,percent_cutoff_vrup=0.3,duration_cutoff=0.9)

segs=neic.segment_data[100]
Dx=25e3
Dy=16e3

#No, lon, lat, depth(km), strike, dip, type, rise_time(s), length(m), width(m)
lon=[]
lat=[]
depth=[]
strike=[]
dip=[]
rise_time=[]
length=[]
width=[]
for i in range(len(segs)):
    fault=segs[i]
    for j in range(len(fault)):
        lon.append(fault[j,1])
        lat.append(fault[j,0])
        depth.append(fault[j,2])
        rise_time.append(0.5)
        strike.append(fault[j,5])
        dip.append(fault[j,6])
        length.append(Dx)
        width.append(Dy)
   
num=arange(1,len(lon)+1,1) 
tp=ones(len(num))
lon=array(lon)
lat=array(lat)
depth=array(depth)
strike=array(strike)
dip=array(dip)
rise_time=array(rise_time)
length=array(length)
width=array(width)
out=c_[num,lon,lat,depth,strike,dip,tp,rise_time,length,width]

savetxt(fout,out,fmt='%d\t%.6f\t%.6f\t%.4f\t%.2f\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\t')