# QVP functional code
# Used to extract QVP data from cf-radial files using Pyart

import numpy as np
import My_pyart_functions
import pyart
import preprocessing_qvp
from netCDF4 import num2date
import kdp_functions_edit as kdpfun
import copy

def altitude_parameter_averaging_qvp_original(radar, elevation, field, azimuth_exclude, verbose=False):
    mask_field = field
    expected_ele = elevation
    
    try:
        sweep = np.where(radar.elevation['data'] == expected_ele)[0][0] / (radar.nrays / radar.nsweeps)

    except:
        print('Elevation not in sweep', radar.time['units'])
        zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
        zero_array[:] = np.nan
        timeofsweep = num2date(radar.time['data'][0],
                               radar.time['units'],
                               radar.time['calendar'])

        return zero_array, zero_array, zero_array, zero_array, timeofsweep
    else:

        if radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]] == expected_ele:
            #try:
            My_pyart_functions.field_fill_to_nan(radar, mask_field)
            #except:
            #    print('No fill value for:', mask_field)
           
            try:
                sweep_ind = (radar.sweep_start_ray_index['data'][sweep], radar.sweep_end_ray_index['data'][sweep])
                sweep_data = radar.fields[mask_field]['data'][sweep_ind[0]:sweep_ind[1]]
                timeofsweep = num2date(np.nanmean(radar.time['data'][sweep_ind[0]:sweep_ind[1]]),
                                       radar.time['units'],
                                       radar.time['calendar'])
                if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
                    sweep_data = np.power(10, sweep_data / 10.0)
                # unfold uPhiDP values
                # remove phi_dp wrap-around:
                #print 'unwrapping phidp ...'
                if mask_field in ['uPhiDP']:
                    rays = sweep_data.shape[0]
                    bins = sweep_data.shape[1]
					
					
					# METEO_THRESH=0.7 OK for meteo QVP should be lower for Biodar
                    METEO_THRESH=0.7
					# flags = np.where(rad2.fields['classification']['data']==1,0,1)
                    flags = np.zeros((rays, bins))
					# generate non-meteo mask:
					#print 'generating non-meteo mask ...'
                    rhohv = radar.fields['RhoHV']['data'][sweep_ind[0]:sweep_ind[1]]
                    (meteoMask) = kdpfun.generate_meteo_mask(rays, bins, flags, rhohv, METEO_THRESH)
                    sweep_data = kdpfun.unwrap_phidp(rays, bins, meteoMask, sweep_data)
                summed = np.zeros(sweep_data.shape)
                counted = np.zeros(sweep_data.shape)
                mask = np.ma.masked_invalid(sweep_data).mask
                inv_mask = np.where(mask, 0, 1)
                summed = np.nansum([summed, sweep_data], axis=0)
                counted = np.nansum([counted, inv_mask], axis=0)
                summed = np.where(counted == 0, np.nan, summed)
                counted = np.where(counted == 0, np.nan, counted)
                mean_values = np.nanmean((summed / counted), axis=0)
                std_values = np.nanstd((summed / counted), axis=0)
                if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:#mask_field == 'dBuZ' or mask_field == 'dBZ':
                    mean_values = 10 * np.log10(mean_values)
                    std_values = 10 * np.log10(std_values)
                observation_count = np.nansum(counted, axis=0)
                if expected_ele == 90.0:
                    altitudes = radar.range['data']
                else:
                    altitudes = radar.fields['scan_altitude']['data'][radar.sweep_start_ray_index['data'][sweep], :]
               
                # clean up the nearest beans influenced by sidelobs
                ranges=radar.range['data']
                mean_values[ranges<400] = np.nan
                
                return altitudes, observation_count, mean_values, std_values, timeofsweep
            except:
                print(radar.time, ' failed for unknown reason')
                zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
                zero_array[:] = np.nan
                timeofsweep = num2date(radar.time['data'][0],
                                       radar.time['units'],
                                       radar.time['calendar'])
                return zero_array, zero_array, zero_array, zero_array, timeofsweep

        else:
            print(radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]])
            zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
            zero_array[:] = np.nan
            timeofsweep = num2date(radar.time['data'][0],
                                   radar.time['units'],
                                   radar.time['calendar'])

            return zero_array, zero_array, zero_array, zero_array, timeofsweep

def altitude_parameter_averaging_qvp(radar, elevation, field, azimuth_exclude, verbose=False):
    
    if verbose:
	    print("in altitude_parameter_averaging_qvp")
		
	    #print "azimuth_exclude"
	    #print azimuth_exclude
	    
    #azimuth_exclude = [45,185]
    mask_field = field
    expected_ele = elevation

    try:
        sweep = np.where(radar.elevation['data'] == expected_ele)[0][0] / (radar.nrays / radar.nsweeps)

    except:
        print('Elevation not in sweep', radar.time['units'])
        zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
        zero_array[:] = np.nan
        timeofsweep = num2date(radar.time['data'][0],
                               radar.time['units'],
                               radar.time['calendar'])

        return zero_array, zero_array, zero_array, zero_array, timeofsweep
    else:
        if radar.elevation['data'][radar.sweep_start_ray_index['data'][int(sweep)]] == expected_ele:
            try:
                My_pyart_functions.field_fill_to_nan(radar, mask_field)
            except:print('No fill value for:', mask_field)
           
            #try:
            sweep_ind = (radar.sweep_start_ray_index['data'][int(sweep)], radar.sweep_end_ray_index['data'][int(sweep)])
            #print "sweep_ind.shape = {}".format(sweep_ind.shape)
            #print sweep_ind
            azimuth_data = radar.azimuth['data'][sweep_ind[0]:sweep_ind[1]]
            #print  " azimuth_data  np.where(azimuth_data>=45&azimuth_data<=185)"
            #print  azimuth_data
            #print np.where((azimuth_data==azimuth_exclude))
            sweep_data = radar.fields[mask_field]['data'][sweep_ind[0]:sweep_ind[1]]
            #sweep_data[np.where((azimuth_data>=45)&(azimuth_data<=185))[0]] = np.nan
            
            #print "sweep_data.shape = {}".format(sweep_data.shape)
            timeofsweep = num2date(np.nanmean(radar.time['data'][sweep_ind[0]:sweep_ind[1]]),
								   radar.time['units'],
								   radar.time['calendar'])
            if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
                sweep_data = np.power(10, sweep_data / 10.0)
			# unfold uPhiDP values
			# remove phi_dp wrap-around:
			#print 'unwrapping phidp ...'
            if mask_field in ['uPhiDP']:
                rays = sweep_data.shape[0]
                bins = sweep_data.shape[1]
                METEO_THRESH=0.7
				# flags = np.where(rad2.fields['classification']['data']==1,0,1)
                flags = np.zeros((rays, bins))
				# generate non-meteo mask:
				#print 'generating non-meteo mask ...'
                rhohv = radar.fields['RhoHV']['data'][sweep_ind[0]:sweep_ind[1]]
                (meteoMask) = kdpfun.generate_meteo_mask(rays, bins, flags, rhohv, METEO_THRESH)
                sweep_data = kdpfun.unwrap_phidp(rays, bins, meteoMask, sweep_data)
            summed = np.zeros(sweep_data.shape)
            counted = np.zeros(sweep_data.shape)
            #rays = sweep_data.shape[0]
            #bins = sweep_data.shape[1]
            #azimuth_mask = np.ones((rays, bins))
            #azimuth_mask[0:azimuth_exclude[0],:] = 0
            #azimuth_mask[azimuth_exclude[1]:360,:] = 0
            #print radar
            #sweep_data[np.where(azimuth_data[sweep_ind[0]:sweep_ind[1]] in azimuth_exclude)[0]] = np.nan
            
            
            # print "EXCLUDE azimuthes"
            # print "azimuth_data"
            # print azimuth_data.astype(int)
            # print len(azimuth_data)
            # print "np.isin(azimuth_data,azimuth_exclude)"
            # print np.isin(azimuth_data.astype(int),azimuth_exclude)
            # #print "azimuth_exclude"
            # #print azimuth_exclude
            
            # print "np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==True)"
            # print np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==True)
            # print len(np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==True))
            
            # print "sweep_data"
            # print sweep_data
            
            #print "np.where(sweep_data != np.nan)"
            #print np.where(sweep_data != np.nan)
            #print "len(np.where(sweep_data != np.nan)[0])"
            #print len(np.where(sweep_data != np.nan)[0])
            
            
            #print "sweep_data[np.where(sweep_data != np.nan)]"
            #print sweep_data[np.where(sweep_data != np.nan)]
            
           # print "sweep_data[0,0]"
            #print sweep_data[0,0]
            
            #print "np.isnan(sweep_data[0,0])"
            #print np.isnan(sweep_data[0,0])
            
            #print "sweep_data[0,0]==float('nan')"
            #print sweep_data[0,0]==float('NaN')
            
            
            
            
            #print "np.isnan(sweep_data)"
            #print np.isnan(sweep_data)
            
            #print "sweep_data[np.isfinite(sweep_data)]"
            #print sweep_data[np.isfinite(sweep_data)]
            #print len(sweep_data[np.isfinite(sweep_data)])
            
            #print "sweep_data[0,0] = np.nan"
            #sweep_data[0,0] = np.nan
            
            #print "np.isnan(sweep_data[0,0])"
            #print np.isnan(sweep_data[0,0])
            
            #print "sweep_data[0,0]==float('nan')"
            #print sweep_data[0,0]==float('nan')
            #print sweep_data[0,0]== np.nan
            
            
            
            #print "np.where(np.isnan(sweep_data))"
            #print np.where(sweep_data != NaN)
            #print "len(np.where(sweep_data != float('nan'))[0])"
            #print len(np.where(sweep_data != float('nan'))[0])
            
            
            #print "sweep_data[np.where(sweep_data != float('nan'))]"
            #print sweep_data[np.where(sweep_data != float('nan'))]
            
            
            #print "sweep_data[np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==True)]"
            
            in_excluded = sweep_data[np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==True)[0],:]
            #print in_excluded[np.isfinite(in_excluded)]
            #print len(in_excluded[np.isfinite(in_excluded)])
            #print "in_excluded.shape"
            #print in_excluded.shape
            #print "np.tile(in_excluded,(sweep_data.shape[1],1))"
            #print np.tile(np.isin(azimuth_data.astype(int),azimuth_exclude).reshape(len(azimuth_data),1),(1,sweep_data.shape[1]))
            #print np.tile(np.isin(azimuth_data.astype(int),azimuth_exclude).reshape(len(azimuth_data),1),(1,sweep_data.shape[1])).shape
            #print "sweep_data.shape[1]"
            #print sweep_data.shape[1]
            
            azimuth_mask = np.tile(np.isin(azimuth_data.astype(int),azimuth_exclude).reshape(len(azimuth_data),1),(1,sweep_data.shape[1]))
            
            
            #print "sweep_data[np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==False)]"
            
            #out_excluded = sweep_data[np.where(np.isin(azimuth_data.astype(int),azimuth_exclude)==False)[0],:]
            #print out_excluded[np.isfinite(out_excluded)]
            #print len(out_excluded[np.isfinite(out_excluded)])
            
            sweep_data = np.ma.masked_where(azimuth_mask,sweep_data)
            
            #print "np.where(azimuth_data[sweep_ind[0]:sweep_ind[1]] in azimuth_exclude)[0] = {}".format(np.where(azimuth_data[sweep_ind[0]:sweep_ind[1]] in azimuth_exclude))
            #print "np.where(azimuth_data[sweep_ind[0]:sweep_ind[1]] in azimuth_exclude)[0] len = {}".format(len(np.where(azimuth_data[sweep_ind[0]:sweep_ind[1]] in azimuth_exclude)[0]))
            
            #print "inv_mask = {}".format(inv_mask)
            
            mask = np.ma.masked_invalid(sweep_data).mask
            inv_mask = np.where(mask, 0, 1)
            #print "np.where(inv_mask==1) = {}".format(len(np.where(inv_mask==1)[1]))
            #print "np.where(inv_mask==0) = {}".format(len(np.where(inv_mask==0)[1]))
            
            summed = np.nansum([summed, sweep_data], axis=0)
            #print "summed"
            #print summed
            counted = np.nansum([counted, inv_mask], axis=0)
            #print "counted"
            #print counted
            summed = np.where(counted == 0, np.nan, summed)
            counted = np.where(counted == 0, np.nan, counted)
            mean_values = np.nanmean((summed / counted), axis=0)
            std_values = np.nanstd((summed / counted), axis=0)
            if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:#mask_field == 'dBuZ' or mask_field == 'dBZ':
                mean_values = 10 * np.log10(mean_values)
                std_values = 10 * np.log10(std_values)
            observation_count = np.nansum(counted, axis=0)
            if expected_ele == 90.0:
                altitudes = radar.range['data']
            else:
                altitudes = radar.fields['scan_altitude']['data'][radar.sweep_start_ray_index['data'][int(sweep)], :]
		   
			# clean up the nearest beans influenced by sidelobs
            ranges=radar.range['data']
            mean_values[ranges<400] = np.nan
			
            return altitudes, observation_count, mean_values, std_values, timeofsweep
		# except:
                # print(radar.time, ' failed for unknown reason')
                # zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
                # zero_array[:] = np.nan
                # timeofsweep = num2date(radar.time['data'][0],
                                       # radar.time['units'],
                                       # radar.time['calendar'])
                # return zero_array, zero_array, zero_array, zero_array, timeofsweep

        else:
            print(radar.elevation['data'][radar.sweep_start_ray_index['data'][int(sweep)]])
            zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
            zero_array[:] = np.nan
            timeofsweep = num2date(radar.time['data'][0],
                                   radar.time['units'],
                                   radar.time['calendar'])

            #return zero_array, zero_array, zero_array, zero_array, timeofsweep

def altitude_parameter_meteo_averaging_qvp(radar, elevation, field, azimuth_exclude, METEO_THRESH=0.7):
    mask_field = field
    expected_ele = elevation

    try:
        sweep = np.where(radar.elevation['data'] == expected_ele)[0][0] / (radar.nrays / radar.nsweeps)

    except:
        print('Elevation not in sweep', radar.time['units'])
        zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
        zero_array[:] = np.nan
        timeofsweep = num2date(radar.time['data'][0],
                               radar.time['units'],
                               radar.time['calendar'])

        return zero_array, zero_array, zero_array, zero_array, timeofsweep
    else:

        if radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]] == expected_ele:
            #try:
            My_pyart_functions.field_fill_to_nan(radar, mask_field)
            #except:
            #    print('No fill value for:', mask_field)
           
            try:
                sweep_ind = (radar.sweep_start_ray_index['data'][sweep], radar.sweep_end_ray_index['data'][sweep])
                this_data = radar.fields[mask_field]['data'][sweep_ind[0]:sweep_ind[1]]
                sweep_data = np.fill_like(this_data,np.nan)
                
                
                timeofsweep = num2date(np.nanmean(radar.time['data'][sweep_ind[0]:sweep_ind[1]]),
                                       radar.time['units'],
                                       radar.time['calendar'])  
                                                            
                rays = sweep_data.shape[0]
                bins = sweep_data.shape[1]
                flags = np.zeros((rays, bins))                      
                rhohv = radar.fields['RhoHV']['data'][sweep_ind[0]:sweep_ind[1]]
                (meteoMask) = kdpfun.generate_meteo_mask(rays, bins, flags, rhohv, METEO_THRESH)
                
                snrH = radar.fields['SNR']['data'][sweep_ind[0]:sweep_ind[1]]
                
                snrV = radar.fields['SNRv']['data'][sweep_ind[0]:sweep_ind[1]]
                
                snrMask = np.fill_like(meteoMask,0)
                snrMask[(snrH > 8) & (snrV > 8)] - 1
                
                sweep_data[(meteoMask == 1) & (snrMask == 1)] = this_data[(meteoMask == 1) & (snrMask == 1)]
				
				
                if mask_field in ['dBuZ', 'dBZ', 'dBuZv', 'dBZv']:
                    sweep_data[(sweep_data < -10) | (sweep_data > 60)] = np.nan
                    sweep_data = np.power(10, sweep_data / 10.0)
                # unfold uPhiDP values
                # remove phi_dp wrap-around:
                #print 'unwrapping phidp ...'
                elif mask_field in ['uPhiDP']:
                    sweep_data = kdpfun.unwrap_phidp(rays, bins, meteoMask, sweep_data)
					
                elif mask_field in ['ZDR']:
                    sweep_data[(sweep_data < -1.5) | (sweep_data > 5)] = np.nan	
                summed = np.zeros(sweep_data.shape)
                counted = np.zeros(sweep_data.shape)
                mask = np.ma.masked_invalid(sweep_data).mask
                inv_mask = np.where(mask, 0, 1)
                summed = np.nansum([summed, sweep_data], axis=0)
                counted = np.nansum([counted, inv_mask], axis=0)
                summed = np.where(counted == 0, np.nan, summed)
                counted = np.where(counted == 0, np.nan, counted)
                mean_values = np.nanmean((summed / counted), axis=0)
                std_values = np.nanstd((summed / counted), axis=0)
                if mask_field == 'dBuZ' or mask_field == 'dBZ':
                    mean_values = 10 * np.log10(mean_values)
                    std_values = 10 * np.log10(std_values)
                observation_count = np.nansum(counted, axis=0)
                if expected_ele == 90.0:
                    altitudes = radar.range['data']
                else:
                    altitudes = radar.fields['scan_altitude']['data'][radar.sweep_start_ray_index['data'][sweep], :]
                return altitudes, observation_count, mean_values, std_values, timeofsweep
            except:
                print(radar.time, ' failed for unknown reason')
                zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
                zero_array[:] = np.nan
                timeofsweep = num2date(radar.time['data'][0],
                                       radar.time['units'],
                                       radar.time['calendar'])
                return zero_array, zero_array, zero_array, zero_array, timeofsweep

        else:
            print(radar.elevation['data'][radar.sweep_start_ray_index['data'][sweep]])
            zero_array = np.zeros((radar.fields[mask_field]['data'].data.shape[1],))
            zero_array[:] = np.nan
            timeofsweep = num2date(radar.time['data'][0],
                                   radar.time['units'],
                                   radar.time['calendar'])

            return zero_array, zero_array, zero_array, zero_array, timeofsweep


def time_height_qvp(list_of_files, elevation, field_list, count_threshold=0, azimuth_exclude = [], verbose=False):
    #verbose = True
    if (verbose):
		#print "in time_height_qvp \n function input is:{} \nwith the number of files {} \n elevation is {} \n field list is {}".format(list_of_files,len(list_of_files),elevation,field_list)
        print(f"in time_height_qvp \n function input is:{list_of_files} \nwith the number of files {len(list_of_files)} \n elevation is {elevation} \n field list is {field_list}")
        print("count_threshold is ")
        print(count_threshold)
       
    result_dict = {}
    stddev_dict = {}
    count_dict = {}
    sweep_times = []
    
    for f, file_ in enumerate(list_of_files):

        radar = pyart.io.read(file_)
        try:
            if (verbose): print(f"{file_}\npreprocessing is running")
            preprocessing_qvp.preprocssing(radar)
        except:
            print('Preprocessing failed for', file_)


        sweep_time = []
        for field in field_list:
                #alts, counts, means, standard_deviations, timeofsweep = altitude_parameter_averaging_qvp(radar, elevation, field, azimuth_exclude = azimuth_exclude)
              
                try:
                    alts, counts, means, standard_deviations, timeofsweep = altitude_parameter_averaging_qvp(radar, elevation, field, azimuth_exclude = azimuth_exclude, verbose=verbose)
                except:
                    if f == 0:
                        print('''Starting file doesn't work, try again''')
                        raise
                    else:
                        #print('Error - need to setup logging on this code')
                        empty_array = result_dict[field]
                        stddev_array = stddev_dict[field]
                        counts_array = count_dict[field]
                        counts = np.zeros(counts_array.shape[0])
                        means = np.zeros(counts_array.shape[0])
                        standard_deviations = np.zeros(counts_array.shape[0])
                        counts[:] = np.nan
                        means[:] = np.nan
                        standard_deviations[:] = np.nan

                        empty_array[:, f] = np.where(counts > count_threshold, means, np.nan)
                        stddev_array[:, f] = np.where(counts > count_threshold, standard_deviations, np.nan)
                        counts_array[:, f] = counts
                        result_dict.update({field: empty_array})
                        stddev_dict.update({field: stddev_array})
                        count_dict.update({field: counts_array})
                else:
                    if f == 0:
                        empty_array = np.zeros((means.shape[0], len(list_of_files)))
                        empty_array[:] = np.nan
                        stddev_array = np.zeros((means.shape[0], len(list_of_files)))
                        counts_array = np.zeros((means.shape[0], len(list_of_files)))
                        stddev_array[:] = np.nan
                        result_dict.update({'alts': alts})
                    else:
                        empty_array = result_dict[field]
                        stddev_array = stddev_dict[field]
                        counts_array = count_dict[field]

                    empty_array[:, f] = np.where(counts > count_threshold, means, np.nan)
                    stddev_array[:, f] = np.where(counts > count_threshold, standard_deviations, np.nan)
                    counts_array[:, f] = counts
                    result_dict.update({field: empty_array})
                    stddev_dict.update({field: stddev_array})
                    count_dict.update({field: counts_array})
                    sweep_time.append(timeofsweep)
        sweep_times.append(max(sweep_time))
        del radar
    return result_dict, stddev_dict, count_dict, sweep_times


def return_units(filename, fields):

    radar = pyart.io.read(filename)
    unit_dictionary = {}
    for field in fields:
        if field in radar.fields.keys():
            unit_dictionary.update({field:radar.fields[field]['units']})
        else:
            print(field, 'not in original files, need to manually add units')
    return unit_dictionary


def return_names(filename, fields):

    radar = pyart.io.read(filename)
    long_dictionary = {}
    short_dictionary = {}
    for field in fields:
        if field in radar.fields.keys():
            print(field)
            long_dictionary.update({field:radar.fields[field]['long_name']})
            print(radar.fields[field])

            try:short_dictionary.update({field: radar.fields[field]['standard_name']})
            except:
                try:short_dictionary.update({field: radar.fields[field]['proposed_standard_name']})
                except:short_dictionary.update({field: radar.fields[field]['long_name']})
        else:
            print(field, 'not in original files, need to manually add names')
    return long_dictionary, short_dictionary

