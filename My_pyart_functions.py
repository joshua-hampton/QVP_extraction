import numpy as np


def field_fill_to_nan(radar, mask_field):
	
	
	if(np.ma.is_masked(radar.fields[mask_field]['data'])):
		mask = radar.fields[mask_field]['data'].mask
		radar.fields[mask_field]['data'] = radar.fields[mask_field]['data'].data
		radar.fields[mask_field]['data'][mask==True] = np.nan

	#radar[mask_field] = np.nan;
	return
