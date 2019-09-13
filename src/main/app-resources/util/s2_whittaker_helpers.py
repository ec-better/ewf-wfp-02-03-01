#from __future__ import absolute_import, division, print_function
import os
import sys 
sys.path.append('/'.join([os.environ['_CIOP_APPLICATION_PATH'], 'util']))
sys.path.append('../util')
import numpy as np
import gdal
import osr
from urlparse import urlparse
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from whittaker import ws2d, ws2doptv, ws2doptvp
from itertools import chain
import cioppy
import array
import geopandas as gpd

def hello_world():
    
    print 'Hello World!'

def get_pipeline_results(pipeline_parameters, search_params):
    
    ciop = cioppy.Cioppy()
    
    if not 'cat' in search_params:
        # add cat out key
        search_params['cat'] = 'out'
    
    creds = '{}:{}'.format(pipeline_parameters['username'], 
                           pipeline_parameters['api_key'])
    
    search = gpd.GeoDataFrame(ciop.search(end_point=pipeline_parameters['end_point'],
                                      params=search_params,
                                      output_fields='link:results',
                                      model='GeoTime',
                                      creds=creds))
    
    fields = 'title,identifier,self,enclosure,cat,cc,wkt,updated,startdate'
    search_result_params = []

    df = pd.DataFrame()

    for index, row in search.iterrows():

        end_point = row['link:results']

        temp_df = pd.DataFrame.from_dict(ciop.search(end_point=end_point,
                                                     params=search_result_params,
                                                     output_fields=fields, 
                                                     model='EOP', 
                                                     creds=creds))

        df = df.append(temp_df, ignore_index=True)
        
    df = df.merge(df.apply(lambda row: analyse_row(row), axis=1), 
              left_index=True,
              right_index=True)
        
    return df

def get_sub_tiles(data_pipeline_results, pipeline_parameters, tiling_factor):
    
    sub_tiles = pd.DataFrame()

    for index, entry in data_pipeline_results.iterrows():

        print entry.jday

        src_ds = gdal.Open(get_vsi_url(entry.enclosure, 
                                       pipeline_parameters['username'], 
                                       pipeline_parameters['api_key']))


        step_x = src_ds.RasterXSize / tiling_factor
        step_y = src_ds.RasterYSize / tiling_factor

        for x in range(0, src_ds.RasterXSize / step_x):

                cols = step_x
                start_x = x * step_x 


                for y in range(0, src_ds.RasterYSize / step_y):

                    temp_dict = dict()

                    rows = step_y
                    start_y = y * step_y

                    #print 'tile_{}_{}'.format(x, y) 

                    #vsi_mem = '/vsimem/t.tif'
                    #ds_src = src_ds

                    #gdal.Translate(vsi_mem, 
                    #               ds_src,
                    #               srcWin=[start_x, start_y, cols, rows],
                    #               bandList=[bands['B04'], bands['B08'], bands['SCL']])

                    #ds_mem = gdal.Open(vsi_mem)

                    #if ds_mem is None:
                    #    raise

                    #geo_transform = ds_mem.GetGeoTransform()
                    #projection = ds_mem.GetProjection()

                    #temp_dict['geo_transform'] = [geo_transform]
                    #temp_dict['projection'] = projection
                    temp_dict['sub_tile'] = 'tile_{}_{}'.format(x, y)
                    temp_dict['start_x'] = start_x
                    temp_dict['start_y'] = start_y
                    temp_dict['cols'] = cols
                    temp_dict['rows'] = rows
                    temp_dict['self'] = entry.self
                    temp_dict['title'] = entry.title
                    temp_dict['day'] = entry.day
                    temp_dict['jday'] = entry.jday
                    temp_dict['enclosure'] = entry.enclosure

                    #for index, band in enumerate(['B04', 'B08', 'SCL']):
                    #    # read the data
                    #    temp_dict[band] = np.array(ds_mem.GetRasterBand(index + 1).ReadAsArray())

                    #temp_dict['MASK'] = ((temp_dict['SCL'] == 2) | (temp_dict['SCL'] == 4) | (temp_dict['SCL'] == 5) | (temp_dict['SCL'] == 6) | (temp_dict['SCL'] == 7) | (temp_dict['SCL'] == 10) | (temp_dict['SCL'] == 11)) & (temp_dict['B08'] + temp_dict['B04'] != 0)

                    #temp_dict['NDVI'] = np.where(temp_dict['MASK'], (temp_dict['B08'] - temp_dict['B04']) / (temp_dict['B08'] + temp_dict['B04']).astype(np.float), np.nan)

                    #for band in ['B04', 'B08', 'SCL']:
                    #    temp_dict.pop(band, None)    

                    pd.Series(temp_dict)

                    sub_tiles = sub_tiles.append(pd.Series(temp_dict), ignore_index=True)  

                    #ds_mem = None
    print('Done!')
    return sub_tiles

def get_vsi_url(enclosure, user, api_key):
    
    
    parsed_url = urlparse(enclosure)

    url = '/vsicurl/%s://%s:%s@%s/api%s' % (list(parsed_url)[0],
                                            user, 
                                            api_key, 
                                            list(parsed_url)[1],
                                            list(parsed_url)[2])
    
    return url 

def analyse_row(row):
    
    series = dict()
    
    series['day'] = row['title'][11:19]
    series['jday'] = '{}{}'.format(datetime.datetime.strptime(series['day'], '%Y%m%d').timetuple().tm_year,
                                   "%03d"%datetime.datetime.strptime(series['day'], '%Y%m%d').timetuple().tm_yday)
    
    return pd.Series(series)


def analyse_subtile(row, parameters):
    
    series = dict()
    
    src_ds = gdal.Open(get_vsi_url(row.enclosure, 
                                   parameters['username'], 
                                   parameters['api_key']))
    
    bands = dict()

    for band in range(src_ds.RasterCount):

        band += 1

        bands[src_ds.GetRasterBand(band).GetDescription()] = band 
        
    vsi_mem = '/vsimem/t.tif'
   
    gdal.Translate(vsi_mem, 
                   src_ds,
                    srcWin=[row.start_x, row.start_y, row.cols, row.rows])
    
    ds_mem = gdal.Open(vsi_mem)
    
    if ds_mem is None:
        raise

    # get the geocoding for the sub-tile
    series['geo_transform'] = [ds_mem.GetGeoTransform()]
    series['projection'] = ds_mem.GetProjection()
    
    for band in ['B04', 'B08', 'SCL']:
        # read the data
        series[band] = np.array(ds_mem.GetRasterBand(bands[band]).ReadAsArray())
    
    series['MASK'] = ((series['SCL'] == 2) | (series['SCL'] == 4) | (series['SCL'] == 5) | (series['SCL'] == 6) | (series['SCL'] == 7) | (series['SCL'] == 10) | (series['SCL'] == 11)) & (series['B08'] + series['B04'] != 0)
    
    series['NDVI'] = np.where(series['MASK'], (series['B08'] - series['B04']) / (series['B08'] + series['B04']).astype(np.float), np.nan)
    
    # remove the no longer needed bands
    for band in ['B04', 'B08', 'SCL']:
        series.pop(band, None)    
    
    return pd.Series(series)

class DateHelper(object):
    """Helper class for handling dates in temporal interpolation."""

    def __init__(self, rawdates, rtres, stres, start=None, nupdate=0):
        """Creates the date lists from input.

        Args:
             rawdates: list of dates from raw file(s)
             rtres: raw temporal resolution
             stres: smooth temporal resolution
             start: start date for custom interpolation
             nupdate: number of points in time to be updated in file (backwards)
            """

        if start:
            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime('%Y%j')
            tdiff = (fromjulian(stop) - fromjulian(rawdates[0])).days
            self.daily = [(fromjulian(rawdates[0]) + datetime.timedelta(x)).strftime('%Y%j') for x in range(tdiff+1)]
            self.target = [self.daily[x] for x in range(self.daily.index(start), len(self.daily), stres)]
            self.target = self.target[-nupdate:]
        else:
            yrmin = int(min([x[:4] for x in rawdates]))
            yrmax = int(max([x[:4] for x in rawdates]))
            daily_tmp = [y for x in range(yrmin, yrmax+2, 1) for y in tvec(x, 1)]
            stop = (fromjulian(rawdates[-1]) + datetime.timedelta(rtres)).strftime('%Y%j')
            self.daily = daily_tmp[daily_tmp.index(rawdates[0]):daily_tmp.index(stop)+1]

            if stres == 5:
                target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in pentvec(x)]
            elif stres == 10:
                target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in dekvec(x)]
            else:
                target_temp = [y for x in range(yrmin, yrmax+1, 1) for y in tvec(x, stres)]
            target_temp.sort()

            for sd in self.daily:
                if sd in target_temp:
                    start_target = sd
                    del sd
                    break
            for sd in reversed(self.daily):
                if sd in target_temp:
                    stop_target = sd
                    del sd
                    break
            self.target = target_temp[target_temp.index(start_target):target_temp.index(stop_target)+1]
            self.target = self.target[-nupdate:]

    def getDV(self, nd):
        """Gets an array of no-data values in daily timesteps.

        Args:
            nd: no-data value

        Returns:
            numpy array with no-data values in daily steps
        """

        return np.full(len(self.daily), nd, dtype='double')

    def getDIX(self):
        """Gets indices of target dates in daily no-data array.

        Returns:
            list with indices of target dates in no-data array
        """

        return [self.daily.index(x) for x in self.target]

def fromjulian(x):
    """Parses julian date string to datetime object.

    Args:
        x: julian date as string YYYYJJJ

    Returns:
        datetime object parsed from julian date
    """

    return datetime.datetime.strptime(x, '%Y%j').date()

def tvec(yr, step):
    """Create MODIS-like date vector with given timestep.

    Args:
        yr: year
        step: timestep

    Returns:
        list with dates
    """

    start = fromjulian('{}001'.format(yr)) + datetime.timedelta()
    tdiff = fromjulian('{}001'.format(yr+1)) - start
    tv = [(start + datetime.timedelta(x)).strftime('%Y%j') for x in range(0, tdiff.days, step)]
    return tv

def dekvec(yr):
    """Create dekadal date vector for given year with fixed days.

    Args:
        yr: year

    Returns:
        list of dates
    """

    return([
        datetime.datetime.strptime(str(yr)+y+x, '%Y%m%d').date().strftime('%Y%j')
        for x in ['05', '15', '25'] for y in [str(z).zfill(2)
                                              for z in range(1, 13)]
    ])


def plot(y, dts,z=None,z_asy=None):
    plt.close()
    xax = [fromjulian(x) for x in dts]
    plt.figure(figsize=(15,8))
    plt.ylim(0,1)
    plt.plot(xax,y,label='y')

    try:
        plt.plot(xax,z,label='z')
    except ValueError:
        pass
    
    try:
        plt.plot(xax,z_asy,label='z_asy')
    except ValueError:
        pass
    
    plt.legend()
    plt.show()
    
def whittaker(ts, dates):
    
    mask = np.ones(len(dates), dtype=bool)

    mask[np.argwhere(np.isnan(ts))] = False
    
    if not np.any(mask.astype(int)):
        
        loptvp = -999
        zvp = 0
        return zvp, loptvp
    
    ts = ts[mask]
    
    w = np.array((ts!=-3000)*1,dtype='double')

    # lrange shall be [-2, 4] in steps of 0.1 (TBC)
    lrange = array.array('d', np.linspace(-2, 4, 60))

    try:
        
        # apply whittaker filter with V-curve
        zv, loptv = ws2doptv(ts, w, lrange)
        zvp, loptvp = ws2doptvp(ts, w, lrange, p=0.9)
    except(IndexError):
        
        loptvp = -999
        zvp = 0

    except(SystemError):
        #loptvp = -999
        #zvp = 0
        print mask

    ndvi_smooth = pd.Series.to_frame(dates).merge(pd.DataFrame(np.column_stack((dates[mask], np.array(zvp))),
                                                 columns=['jday', 'sndvi']), 
                                    on='jday',
                                    how='left').drop(['jday'], axis=1).values.T[0]
    
    
    return tuple(chain(np.append([loptvp], ndvi_smooth)))  
