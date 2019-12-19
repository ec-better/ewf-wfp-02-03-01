from __future__ import absolute_import, division, print_function

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
from whittaker import ws2d, ws2doptv, ws2doptvp, lag1corr
from itertools import chain
import cioppy
import array
import geopandas as gpd

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
    
    fields = 'title,identifier,self,enclosure,cat,cc,wkt,updated,startdate,vs:"tileid"'
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

        print(entry.jday)

        src_ds = gdal.Open(get_vsi_url(entry.enclosure, 
                                       pipeline_parameters['username'], 
                                       pipeline_parameters['api_key']))


        step_x = src_ds.RasterXSize // tiling_factor
        step_y = src_ds.RasterYSize // tiling_factor

        for x in range(0, src_ds.RasterXSize // step_x):

                cols = step_x
                start_x = x * step_x 


                for y in range(0, src_ds.RasterYSize // step_y):

                    temp_dict = dict()

                    rows = step_y
                    start_y = y * step_y
                    
                    temp_dict['sub_tile'] = 'tile_{}_{}_{}'.format(x, y, entry.title[38:44])
                    temp_dict['start_x'] = start_x
                    temp_dict['start_y'] = start_y
                    temp_dict['cols'] = cols
                    temp_dict['rows'] = rows
                    temp_dict['self'] = entry.self
                    temp_dict['title'] = entry.title
                    temp_dict['day'] = entry.day
                    temp_dict['jday'] = entry.jday
                    temp_dict['enclosure'] = entry.enclosure

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

def analyse_merge_row(row, band_to_process):
    
    series = dict()

    if 's_' in row.enclosure:
        output_type = 's'
    
    elif 'original_{}'.format(band_to_process) in row.enclosure:
        output_type = 'original_{}'.format(band_to_process)
    
    elif 'lag1_{}'.format(band_to_process) in row.enclosure:
        output_type = 'lag1'
        
    else: 
        output_type = band_to_process

    series['output_type'] = output_type
    series['tile'] = os.path.basename(row.enclosure).split('_')[-1].split('.')[0]
    return pd.Series(series)    

def analyse_subtile(row, parameters, band_to_analyse):
    
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

    if band_to_analyse == 'NDVI':
        
        for band in ['B04', 'B08', 'SCL']:
            # read the data
            series[band] = np.array(ds_mem.GetRasterBand(bands[band]).ReadAsArray())

        series['MASK'] = ((series['SCL'] == 2) | (series['SCL'] == 4) | (series['SCL'] == 5) | (series['SCL'] == 6) |
                          (series['SCL'] == 7) | (series['SCL'] == 10) | (series['SCL'] == 11)) & (series['B08'] + series['B04'] != 0)

        series['NDVI'] = np.where(series['MASK'], (series['B08'] - series['B04'])/(series['B08'] + series['B04']), np.nan)

        ### Added to delete non-interpretable data of the NDVI. Check the source to understand the problem.
        series['NDVI'] = np.where(((series['NDVI'] > 1) | (series['NDVI'] < -1)), np.nan, series['NDVI'])

        # remove the no longer needed bands
        #for band in ['B04', 'B08', 'SCL']:
        #    series.pop(band, None)
        for band in ['B04', 'B08']:
            series.pop(band, None)

    else:
        for band in [band_to_analyse, 'SCL']:
            # read the data
            series[band] = np.array(ds_mem.GetRasterBand(bands[band]).ReadAsArray())
            
        series[band_to_analyse] = np.array(ds_mem.GetRasterBand(bands[band_to_analyse]).ReadAsArray())
        
        series['MASK'] = ((series['SCL'] == 2) | (series['SCL'] == 4) | (series['SCL'] == 5) | (series['SCL'] == 6) |
                          (series['SCL'] == 7) | (series['SCL'] == 10) | (series['SCL'] == 11))
        
        series[band_to_analyse] = np.where(series['MASK'], series[band_to_analyse], np.nan)
    
    
    series['SCL_mask'] = ((series['SCL'] == 2) | (series['SCL'] == 4) | (series['SCL'] == 5) | (series['SCL'] == 6) |
                          (series['SCL'] == 7) | (series['SCL'] == 10) | (series['SCL'] == 11))
        
    ds_mem.FlushCache()

    return pd.Series(series)

def fromjulian(x):
    """
    Parses julian date string to datetime object.

    Args:
        x: julian date as string YYYYJJJ

    Returns:
        datetime object parsed from julian date
    """

    return datetime.datetime.strptime(x, '%Y%j').date()
    
def generate_dates(startdate_string=None, enddate_string=None, delta=5):
    """
    Generates a list of dates from a start date to an end date.

    Args:
        startdate_string: julian date as string YYYYJJJ
        enddate_string: julian date as string YYYYJJJ
        delta: integer timedelta between each date

    Returns:
        list of string julian dates YYYYJJJ
    """

    
    startdate = datetime.datetime.strptime(startdate_string, '%Y%j').date()
    enddate = datetime.datetime.strptime(enddate_string, '%Y%j').date()
    
    date_generated = [startdate + datetime.timedelta(days=x) for x in range(0, (enddate-startdate).days+delta, delta)]
    
    datelist = ['{}{:03d}'.format(x.year, x.timetuple().tm_yday) for x in date_generated]

    return datelist

def whittaker(ts, date_mask):
    """
    Apply the whittaker smoothing to a 1d array of floating values.

    Args:
        ts: array of floating values
        date_mask: full list of julian dates as string YYYYJJJ

    Returns:
        list of floating values. The first value is the s smoothing parameter
    """
    
    #  mask is True when the value is np.nan
    mask = np.isnan(ts)
    
    # the output is an  array full of np.nan by default
    ndvi_smooth = np.array([np.nan]*len(date_mask))
    
    # check if all values are np.npn
    if not mask.all():
        
        # parameters needed for the first smoothing without interpolation
        ts_not_nan = ts[~mask]

        w = np.ones(len(ts_not_nan), dtype='double')

        lrange = array.array('d', np.linspace(-2, 4, 61))
        
        try: 
            # apply whittaker filter with V-curve
            zv, loptv = ws2doptvp(ts_not_nan, w, lrange, p=0.90)
            #parameters needed for the interpolation step
            dvec = np.zeros(len(date_mask))
            
            w = np.ones(len(ts), dtype='double')
            
            w[mask] = 0
            
            # adding new dates with no associated product to the weights
            for idx, el in enumerate(date_mask):
                if not el:
                    w = np.insert(w, idx, 0)

            dvec[w==1] = zv
            
            # apply whittaker filter with very low smoothing to interpolate
            ndvi_smooth = ws2d(dvec, 0.0001, w)
            
            # Calculates Lag-1 correlation
            #def test_lag1corr(self):
            #    """Test lag-1 correlation function."""
            #    self.assertAlmostEqual(lag1corr(self.y[:-1], self.y[1:], -3000.0), self.data['lag1corr'])
            lag1 = lag1corr(ts[:-1], ts[1:], -999)

        except Exception as e:
            loptv = -999
            lag1 = -999
            print(e)
            print(mask)

    else:
        loptv = -999
        lag1 = -999
        
    return tuple(np.append(np.append(loptv,lag1), ndvi_smooth))

def cog(input_tif, output_tif):
    
    translate_options = gdal.TranslateOptions(gdal.ParseCommandLine('-co TILED=YES ' \
                                                                    '-co COPY_SRC_OVERVIEWS=YES ' \
                                                                    ' -co COMPRESS=LZW'))

    ds = gdal.Open(input_tif, gdal.OF_READONLY)

    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'DEFLATE')
    ds.BuildOverviews('NEAREST', [2,4,8,16,32])
    
    ds = None

    ds = gdal.Open(input_tif)
    gdal.Translate(output_tif,
                   ds, 
                   options=translate_options)
    ds = None

    os.remove('{}.ovr'.format(input_tif))
    os.remove(input_tif)
