from bs4 import BeautifulSoup
import urllib2
import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from obspy.core import UTCDateTime
from numpy import array,where

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2)


for k in range(len(neic.IDs)):
    evid=neic.IDs[k]
    url="https://earthquake.usgs.gov/earthquakes/eventpage/us"+evid
    page = urllib2.urlopen(url)
    soup = BeautifulSoup(page.read())
    title=array(soup.find('title').decode().split())
    position=where(title=='-')[0]+1
    region=''
    for i in range(position,len(title)):
        if i==0:
            region=title[position]
        elif i==len(title)-1:
            region=region+'-'+title[i].replace('</title>','')
        else:
            region=region+'-'+title[i]

    ev_string=str(UTCDateTime(neic.origin_times[k]).year)+'-M'+str(neic.magnitudes[k])[0:3]+region
    print ev_string

    