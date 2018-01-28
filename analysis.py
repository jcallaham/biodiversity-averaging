#!/usr/bin/python3

"""
Jared Callaham
6/28/17

Temperature averaging to account for coordinate uncertainty
    - Load NOAA data and csv of VertNet records
    - For each collection record, calculate an average temperature for the year of collection
    - Calculate correlation between temperature and mass for each species

"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset

records_path = 'VertNet.csv'            # File containing cleaned-up collection records       
save_filename = 'averaged_temps.csv'

R_earth = 6.371E6           # Radius of earth in meters

##############################################################
# Load temperature file and average most recent year of data
##############################################################
nc_file = 'air.mon.mean.v401.nc'

print("Loading temperature data....")
fh = Dataset(nc_file, mode='r')

# Load data into separate 1D numpy arrays
lons = fh.variables['lon'][:]
d_lon = lons[1]-lons[0]
lon_offset = min(lons) - d_lon/2

lats = fh.variables['lat'][:]
d_lat = lats[1]-lats[0]
lat_offset = min(lats) - d_lat/2

time = fh.variables['time'][:] 
time = time.astype(int)                # 'hours since 1900-1-1 0:0:0'

air = fh.variables['air'][:]           # Air temperature - index by [time, lat, lon]
air_units = fh.variables['air'].units

fh.close()

# Average temperature over past year
print("Averaging data for last year of collection")

# Find time value corresponding to the last year in measurements
one_year_ago = min( [t for t in time if (time[-1]-t)<(24*365)] )
start_index = np.nonzero(time==one_year_ago)
start_index = np.squeeze(start_index)      # Remove extra dimensions

# Keep average air temperature over last year in numpy MaskedArray
#   in order to track locations with missing data
temperature = np.mean(air[start_index:, :, :], 0) 

# Plot temperatures (looks cool)
plt.imshow(temperature, cmap='coolwarm', interpolation='nearest')
plt.show()

#############################################
# Accessing temperature data through pixels
################################################

# Width and height of pixels in meters
width = lambda lat_idx: R_earth*np.cos( lats[lat_idx]*(np.pi/180) )/len(lons)
height = R_earth/len(lons)

def coords_to_indices(this_lat, this_lon):
    """ Return index pair corresponding to pixel of lat, lon """
    lon_idx = int( (this_lon + lon_offset) /d_lon )
    lat_idx = int( (this_lat + lat_offset) /d_lat )
    return (lat_idx, lon_idx)


def get_temp( lat, lon ):
    """ Return temperature measurement for pixel corresponding to lat, lon """
    lat_idx, lon_idx = coords_to_indices(lat, lon)
    if temperature.mask[ lat_idx, lon_idx ]:
        return np.nan
    else:
        return temperature[ lat_idx, lon_idx ]



###########################################
# Weighted averaging for each record
###########################################

def dist(loc1, loc2):
    """ Distance in meters between the pixel centers corresponding to the pairs
    of (lat, lon) indices provided"""
    dx, dy = np.array(loc1) - np.array(loc2)    # Center-center distance in pixels
    return np.sqrt( ( dx*width(loc1[0]) )**2 + (dy*height)**2 )  # Scale by pixel width and height


def pdf(r, sigma, R_max=None):
    """ Estimated probability of collection at this distance from
            recorded coordinates, as a function of distance in meters.
        Note: Normalization depends on number of valid pixels, so
            this must be done externally to this function"""
    if R_max==None: R_max = sigma
    if r > R_max: return 0
    return 1        # Uniform distribution


def get_pixels(loc, R_max):
    """Search for indices of pixels within R_max of 'loc'.
    Probably kind of an inefficient search, since it is based on a square
    Assumes width of pixels is approximately constant near 'loc'
    I guess it also assumes there aren't any measurements sufficiently close
        to the poles to lead to problems with non-Euclidean geometry....

        Inputs:
            loc: (lat, lon) indices
            R_max: maximum distance from loc to return pixels
        Output:
            nearby: list of (lat, lon) indices corresponding to nearby pixels
    """

    # Number of pixels to search in vicinity of 'loc' (width will be less than or equal to height)
    to_search = int( R_max//width(loc[0]) ) + 1

    nearby = []         # Will be filled with pixels within R_max of 'loc'
    for d_lon in range(-to_search, to_search):     # Longitude steps
        lon_idx = (loc[1]+d_lon)%len(lons)   # Calculate neighbor index, including possibility of wrapping around
        for d_lat in range(-to_search, to_search): # Latitude steps
            lat_idx = (loc[0]+d_lat)%len(lats)
            if dist( loc, (lat_idx, lon_idx) ) < R_max:
                nearby.append( (lat_idx, lon_idx) )

    return np.array(nearby)


def avg_temp(this_lat, this_lon, sigma, R_max=None):
    """ Find average temperature at this latitude and longitude according to the
    specified uncertainty sigma and maximum distance R_max """
    if R_max==None: R_max = sigma

    loc = coords_to_indices(this_lat, this_lon)      # (lat, lon) indices
    nearby = get_pixels(loc , R_max )               # Array of (lat, lon) indices

    # Find temperatures of nearby pixels
    nearby_temps = [temperature[ nearby[i][0], nearby[i][1] ] for i in range(len(nearby))]
    # Convert to masked array to account for lack of measurements (e.g. near coast)
    masked_temps = np.ma.masked_equal(nearby_temps, value=temperature.fill_value ) 

    # If using a uniform distribution, could now just return np.mean(masked_temps)
    """
    nearby_temps = masked_temps[ ~masked_temps.mask ]      # Now contains only valid measurements
    nearby = nearby[ ~masked_temps.mask ]      # Forget about nearby pixels without temperature measurements

    if len(nearby)==0:
        return get_temp(this_lat, this_lon)  # If no valid neighbors left, forget about averaging

    # Now calculate normalization
    N = sum( [pdf( dist(loc, nearby[i]), sigma, R_max) for i in range(len(nearby)) ] )
    T_avg = (1/N)*sum( [pdf( dist(loc, nearby[i]), sigma, R_max)*
            temperature[ nearby[i][0], nearby[i][1] ] for i in range(len(nearby)) ] )
    return T_avg
    """

    return np.mean(masked_temps)
 

# Read lat, lon, coordinate uncertainty in from database and calculate avg temp
#   Add this as column to dataframe and save dataframe, then  plot results

# Load VertNet collection records into pandas DataFrame
records = pd.read_csv(records_path)

# Add columns to store temperature data
records['coord_temp'] = np.zeros( records.shape[0] )   # Temperature with no location uncertainty (for comparison)
records['avg_temp'] = np.zeros( records.shape[0] )     # Averaged temperature

records = records.rename( columns={'coordinateuncertaintyinmeters': 'uncertainty'} )  # For code convenience
records.loc[ records['uncertainty'] > R_earth, 'uncertainty' ] = np.nan  # If uncertainty greater than Earth's radius, make NaN

records = records.dropna( subset=['uncertainty'] )  # For purposes of seeing effect of averaging

counter = 0
# Should be able to do this with DataFrame.apply()
for index, row in records.iterrows():
    lat, lon, sigma = records.loc[ index, ['latitude', 'longitude', 'uncertainty'] ]
    #print("{0:0.2f}\t{1:0.2f}\t{2}".format(lat, lon, sigma) ) 
    records.loc[ index, 'coord_temp' ] = get_temp( lat, lon )
    records.loc[ index, 'avg_temp' ] = avg_temp( lat, lon, sigma )  # Have already dropped NaN sigmas

    counter += 1
    if (counter%1000 == 0): print("Entry {0} of {1}".format(counter, records.shape[0]) )


# Save DataFrame including average and coordinate temperature
records.to_csv(save_filename)
 

#####################################################
# Plot histogram of differences due to averaging
#######################################################

records = pd.read_csv(save_filename) 

records = records.dropna( subset=['coord_temp', 'avg_temp'] )
diff = records.loc[:, 'coord_temp'] - records.loc[:, 'avg_temp']
print(len(diff))

plt.hist(diff, bins=np.linspace(-1, 1, num=30))
plt.ylabel('Frequency')
plt.xlabel('Shift in temperature due to spatial average')
plt.show()  
