# QVP extraction code
# Extract QVP data using predefined functions and store the data in a netcdf file
# Written by David Dufton, Nov. 2016
# Updated by M. Lukach, May 2018
import sys
import os
from optparse import OptionParser

import numpy as np
from netCDF4 import Dataset
from netCDF4 import date2num
import qvp_functions_nu as qvp_functions
import sys
import glob
import datetime as dt
import calendar
import time

import pytz

def main():
    print("qvp_extraction.py is running")
    
    usage = """%prog [options] /input/directory/*.nc elevation /output/dyrectory/
    Use the option -h or --help to get all possible options
    """
    #determine basename of the application
    appbasename = os.path.basename(os.path.splitext(sys.argv[0])[0])
        
    #determine the basedir  of the application
    basedir = os.path.dirname(sys.argv[0]) + os.sep
    if(basedir == "/"):basedir = "./"
    
    #the default logfile is located in the same directory as the program and
    #has the program name with the extension '.log'
    dfltlogfile = basedir + appbasename + '.log'
    
    # default values, that will be changed by the parser
    tstart = '20170517T000000' 
    tstop = '20170517T235959'
    zoom_interval = ('20170517T000000','20170517T235959') #'KDP_VULPIANI_WRADLIB','KDP_VULPIANI',
    field_list = ['dBuZ','dBuZv','dBZv','dBZ','KDP_UKMO','KDP','ZDR','ZDRu','RhoHV','RhoHVu','uPhiDP','PhiDP','SQI','SNR','V','Vu','W','DOP','DOPu','temperature']
    
    # Other examples of the fields lists:
    #field_list = ['dBuZ','dBuZv','dBZv','dBZ','KDP_UKMO','KDP','ZDR','ZDRu','RhoHV','RhoHVu','uPhiDP','PhiDP','SQI','SNR','V','Vu','W','DOP','DOPu','temperature_4']
    #field_list = ['dBuZ','dBZ','dBZ_ac','ZDR_ac','KDP_UKMO','KDP','ZDR','ZDRu','RhoHV','RhoHVu','uPhiDP','PhiDP','SQI','SNR','V','Vu','W','DOP','DOPu','temperature_2']
    #field_list = ['dBuZ','dBZ','KDP_UKMO','KDP','ZDR','ZDRu','RhoHV','RhoHVu','uPhiDP','SQI','SNR','V','Vu','W','DOP','DOPu','temperature_2']
    
    # parsing options and arguments 'KDP_UKMO',
    parser = OptionParser(usage=usage)
    parser.add_option("-d", "--debug", action="store_true", dest="debug", default=False, help="output debug messages during program execution")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="print messages (debug,info,warning,...) about the program execution to the console (stderr)")
    parser.add_option("-s", "--start time (included)", dest="tstart", default = tstart, help="date and time in format: yyyymmdd'T'HHMMSS")
    parser.add_option("-e", "--end time (included)", dest="tstop", default = tstop, help="date and time in format: yyyymmdd'T'HHMMSS")
    parser.add_option("-f", "--fields for QVP generation", dest="field_list", default = field_list, help="list of variables for QVP calculation")
    parser.add_option("-z", "--zoom time-interval", dest="zoom_interval", default = zoom_interval, help="Time interval to zoom, Formatted as (yyyymmdd'T'HHMMSS,yyyymmdd'T'HHMMSS)")
    parser.add_option("-c", "--count threshold", dest="count_threshold", default = 0, help="Minimal number of points for the mean value calculation at each range")
    parser.add_option("-a", "--azimuth bounds to exclude", dest="azimuth_bounds_to_exclude", default = None, help="Azimuths, that should be excluded in the mean value calculation. Format [star1,end1;start2,end2]")
    	
    #parsing parameters
    (options, files) = parser.parse_args()
    debug = options.debug
    verbose = options.verbose
    tstart = options.tstart
    tstop = options.tstop
    field_list = options.field_list
    zoom_interval = options.zoom_interval
    count_threshold = int(options.count_threshold)
    
    if(options.azimuth_bounds_to_exclude != 'None'): 
        azimuth_bounds_to_exclude = options.azimuth_bounds_to_exclude.split(',')
        azimuth_exclude = np.arange(int(azimuth_bounds_to_exclude[0]),int(azimuth_bounds_to_exclude[1]))
    else:
        azimuth_exclude = None
    
    #Make aware datetime objects from the start and stop dates (with time zone info included in the object)
    utc = pytz.timezone("UTC")
    start_datetime = dt.datetime(int(tstart[0:4]),int(tstart[4:6]), 
                                int(tstart[6:8]), int(tstart[9:11]),
                                int(tstart[11:13]),int(tstart[13:15]), tzinfo=utc)
    stop_datetime = dt.datetime(int(tstop[0:4]),int(tstop[4:6]), 
                                int(tstop[6:8]), int(tstop[9:11]),
                                int(tstop[11:13]),int(tstop[13:15]), tzinfo=utc)
    if verbose:
        print("start date is")
        print(tstart)
        print("stop date is")
        print(tstop)
        print("both dates belong to the same day")
        print(stop_datetime - start_datetime < dt.timedelta(1))
    
    # Define event date TODO define list of event's dates
    if(stop_datetime - start_datetime < dt.timedelta(1)):event = start_datetime.strftime('%Y%m%d') 
    else: event = start_datetime.strftime('%Y%m%d') 
    
    if(len(zoom_interval)==2):
        #Make aware datetime objects from the zoom_interval dates (with time zone info included in the object)
        zoom_start = dt.datetime(int(zoom_interval[0][0:4]),int(zoom_interval[0][4:6]),
                                int(zoom_interval[0][6:8]), int(zoom_interval[0][9:11]),
                                int(zoom_interval[0][11:13]),int(zoom_interval[0][13:15]), tzinfo=utc)
        zoom_end = dt.datetime(int(zoom_interval[1][0:4]),int(zoom_interval[1][4:6]),
                               int(zoom_interval[1][6:8]), int(zoom_interval[1][9:11]),
                               int(zoom_interval[1][11:13]),int(zoom_interval[1][13:15]), tzinfo=utc)
    else:
        if verbose: print("Zoom interval should be a list with two dates. The event's start and end are used instead.")
        zoom_start = start_datetime
        zoom_end = stop_datetime
	
    folder_glob_spec = sys.argv[1]
    if(folder_glob_spec[-1]=="/"):folder_glob_spec = os.path.dirname(folder_glob_spec)
    if verbose: print(f"First argument is {folder_glob_spec}")
    elevation = float(sys.argv[2])
    if verbose: print(f"Second argument is {elevation}")
    #'./20170203/ver/*.nc'
    if (elevation == 90.0):folder_with_files = f'{folder_glob_spec}/{event}/ver/*.nc'
    else:folder_with_files = f'{folder_glob_spec}/{event}/*.nc'
    if verbose: print(f"Input is {folder_with_files}")
    
    
    output_glob_spec = sys.argv[3]
    if(output_glob_spec[-1]=="/"):output_glob_spec = os.path.dirname(output_glob_spec)
    output_dir = f"{output_glob_spec}/{event}_QVP"
    if not os.path.exists(output_dir):os.makedirs(output_dir)
    
    output_file = f'{output_dir}/{event}_QVP_{int(elevation)}deg.nc' #20170517_QVP_20deg
    if verbose: 
        print(f"Third argument is {output_glob_spec}")
        print(f"Output will be placed in {output_file}")
        print(f"count_threshold is {count_threshold}")
    
    # Another way to ignor the count threshold
    #count_threshold = 0
    
    file_list = glob.glob(folder_with_files)
    file_list.sort()
    if verbose:
        print("file_list")
        print(file_list)
        print('First file is:', file_list[0])
        print('Last file is:', file_list[-1])
        #print('There are %i files') % (len(file_list))
        print(f'There are {len(file_list)} files')
	
    result_dict, stddev_dict, count_dict, sweep_times = qvp_functions.time_height_qvp(file_list,
                                                                         elevation,
                                                                         field_list,
                                                                         count_threshold,azimuth_exclude=azimuth_exclude,verbose=verbose)
	
	
    qvp_output = Dataset(output_file, 'w', format='NETCDF4')
    # Create dimensions
    qvp_output.createDimension('Height', result_dict['alts'].shape[0])
    qvp_output.createDimension('Time', len(file_list))
    # Create coordinate variables for 2 dimensions
    times = qvp_output.createVariable('Time', np.float64,('Time',))
    heights = qvp_output.createVariable('Height', np.float32, ('Height',))
    
    time_units = f'seconds since {sweep_times[0].year}-01-01T00:00:00Z'
    number_times = [date2num(sweeptime,time_units,calendar='gregorian') for sweeptime in sweep_times]
    
    ## Add values to the dimension variables
    heights[:] = result_dict['alts']
    times[:] = np.array(number_times)
    times.units = time_units
    times.calendar = "gregorian"
    heights.units = "metres above sea level"
    
    unit_dictionary = qvp_functions.return_units(file_list[0], field_list)
    long_names, short_names = qvp_functions.return_names(file_list[0], field_list)
    
    # Create the field variables
    for field in field_list:
        group_data_construct_means = f'/{field}/Means'
        group_data_construct_deviations = f'/{field}/StdDevs'
        group_data_construct_counts = f'/{field}/Counts'
        temp_means = qvp_output.createVariable(group_data_construct_means, np.float32, ('Height','Time'))
        temp_stds = qvp_output.createVariable(group_data_construct_deviations, np.float32, ('Height', 'Time'))
        temp_counts = qvp_output.createVariable(group_data_construct_counts, np.float32, ('Height', 'Time'))
        temp_means[:] = result_dict[field]
        temp_stds[:] = stddev_dict[field]
        temp_counts[:] = count_dict[field]
        temp_means.long_name = f'Quasi-vertical mean {field} at an elevation angle of {round(elevation,1)} degrees'
        stdev_name = f'Quasi-vertical standard deviation of {field} at an elevation angle of {round(elevation,1)} degrees'
        count_name = f'Quasi-vertical number of observations of {field} at an elevation angle of {round(elevation,1)} degrees'
        temp_stds.long_name = stdev_name
        temp_counts.long_name = count_name
        if field in unit_dictionary.keys():
            qvp_output[field].units = unit_dictionary[field]
            qvp_output[field].long_name = long_names[field]
            qvp_output[field].standard_name = short_names[field]
        else:
            qvp_output[field].units = 'Default'
            qvp_output[field].long_name = 'Default'
            qvp_output[field].standard_name = 'Default'
    qvp_output.elevation = elevation
	
    qvp_output.close()
    
    if verbose:print(f'New netcdf file created at {output_file}')

	

	
#-----------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
	
