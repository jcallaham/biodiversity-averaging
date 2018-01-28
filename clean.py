#!/usr/bin/python3

"""
Loads multiple VertNet csv files into database, extracts desired information,
    and saves results into one single csv file

Jared Callaham
6/4/17
"""

import pandas as pd
import sqlite3
import glob
import os
import numpy as np
import pytaxize
from netCDF4 import Dataset
import time


#####################################################
# Load raw VertNet files into sqlite database
#   (at this stage there are approximately 15 million entries)
#####################################################

print("Loading VertNet database")
data_dir = 'vertnet-temp/'
csv_filenames = glob.glob(data_dir+'*.csv')       # Assume that all .csv files in this directory are raw VertNet databases
table_names = [name.split('_')[2][:-4] for name in csv_filenames]    # Filename format is 'vertnet_latest_{table_name}.csv'

    
db_filename = data_dir+'vertnet.db'          # Temporary sqlite database to extract desired columns
out_filename = 'all_vertnet.csv'    # Single output filename

"""
chunksize=5000                      # Number of rows to load at a time from input files (to avoid running out of memory)
         
# Build database out of csv files
with sqlite3.connect(db_filename) as cnx :
    for (csv_name, table_name) in zip(csv_filenames, table_names):
        row = 0     # Counter to display progress
        print("Processing file {0}".format(csv_name))
        df = pd.read_csv(csv_name, chunksize=chunksize, dtype=str)
        for chunk in df:
            chunk.to_sql(table_name, cnx, if_exists='append')
            row += chunksize
            print( 'Row {0}'.format(row) )


#################################################
# Load and query the sqlite database
################################################

with sqlite3.connect(db_filename) as cnx:
    c = cnx.cursor()

    # Build query to extract desired information
    query = ""
    table_names = iter(table_names)
    for name in table_names:
        query += "SELECT scientificname, class, family, year, decimallongitude, decimallatitude, coordinateuncertaintyinmeters, massing, citation, license, isfossil, underivedlifestage FROM {0} WHERE massing IS NOT NULL UNION ALL ".format(name)
    query = query[:-len(" UNION ALL ")]    # Get rid of last " UNION ALL "

    print("Executing query...")
    # Execute query and convert to pandas DataFrame
    data = [row for row in c.execute(query)]

    print("Formatting and saving DataFrame...")
    columns = [c.description[i][0] for i in range(len(c.description))]
    data = pd.DataFrame(data=data, columns=columns)


# Save DataFrame to csv
data.to_csv(out_filename)

# Remove temporary database file
#os.remove(db_filename)


"""

######################################################################
# Clean taxonomy records
#       (could be done before saving original DataFrame, but saving and loading doesn't take long)
#######################################################################

print("Cleaning taxonomy records....")

tax_filename = 'clean_taxonomy_py.csv'
cleaned_filename = 'VertNet.csv'

data = pd.read_csv(out_filename, dtype=str) 
        
orig_names = data['scientificname'].as_matrix()   # Contains NaN values
print(type(orig_names))
orig_names[pd.isnull(orig_names)] = "nan"          # Will be converted to np.nan in loop

# Get list of all reported names and the indices to reconstruct the array
orig_names, orig_idx = np.unique(orig_names, return_inverse=True)
print(orig_names[:10])
end = len(orig_names)

cleaned_names = np.empty(end, dtype=object)           # Taxize matches for the reported names

print("Loaded data")

# Clean taxonomy names (so that species identifiers match)  Start 1:22
for i in range(end):
    words = orig_names[i].split()
    named = False
    if len(words) >= 2:
        name = words[0]+' '+words[1]       # Remove subspecies, if exists 

        # Search for match with pytaxize
        matched = False
        while not matched:
            try:
                tax_match = pytaxize.gnr_resolve(name)[0]
                matched = True
            except ConnectionResetError:
                print("No connection. Waiting...")
                time.sleep(60)      # If internet fails, wait and try again
        
        if (len(tax_match)!=0) and len(tax_match[0]['canonical_form'].split()) >= 2:
            matches = np.array([tax_match[j]['canonical_form'] for j in range(min([len(tax_match), 5]))])
            best_match = matches[0].split()
            best_match = best_match[0]+' '+best_match[1]   # Take only genus and species
            if name == best_match:
                print("Submitted name accepted")
                cleaned_names[i] = name   # Submitted name equals match
                named = True
            elif all( matches[1:] == matches[0] ):   # Slight error here: should split matches to first two words... may not matter
                print("Submitted name corrected")
                cleaned_names[i] = best_match   # Resolved names are consistent but don't match submission
                named = True

        #cleaned_names.append(name)  # Delete when using pytaxize
                
    if not named:       # Less than two words, no match, etc.
        print("{0} rejected".format(orig_names[i]))
        cleaned_names[i] = np.nan
    print("Name {0} of {1}: {2}".format(i, end, cleaned_names[i]))
    
    # Save in case taxize interrupted
    if (i%100 == 0):
        data['genus_species'] = cleaned_names[orig_idx]
        data.to_csv(tax_filename)
        
# Remove 'Environmental Halophage' listings
cleaned_names[ np.nonzero( cleaned_names=="Environmental Halophage" ) ] = np.nan
    
# Replace VertNet entries with corresponding cleaned-up taxonomy
data['genus_species'] = cleaned_names[orig_idx]  # Reverse index "cleaned_names"

print(data.shape)
data = data.loc[data['genus_species'] != 'nan', :]      # Remove entries with no name

"""
(721486, 14)
sys:1: DtypeWarning: Columns (4) have mixed types. Specify dtype option on import or set low_memory=False.
(721486, 14)
Subsetting by range and valid coordinates...
(476254, 14)
(371911, 14)
Subsetting by first adult...
(337234, 14)
(333421, 14)
"""


del data['scientificname']

data.to_csv(tax_filename)

####################
# Subset Data
#   Result should be 275,000 individuals (as per paper)
#   TO DO... CLEAN THIS UP!
#######################
data = pd.read_csv(tax_filename)
print(data.shape)
print("Subsetting by range and valid coordinates...")

# Convert some columns to numbers
data = data.rename( columns={'decimallatitude': 'latitude', 'decimallongitude': 'longitude',
        'Unnamed: 0': 'row_index', 'massing': 'mass'} )
data.loc[:, ['latitude','longitude', 'mass', 'year']] = \
        data.loc[:, ['latitude','longitude', 'mass', 'year']].apply(pd.to_numeric, errors='coerce')

# Remove coordinates outside of range and transform longitudes to expected range
data = data[ (data['latitude'] < 90) & (data['latitude'] > -90) ]
data = data[ (data['longitude'] < 180) & (data['longitude'] > -180) ]
data.loc[data['longitude'] < 0, 'longitude'] = data.loc[data['longitude'] < 0, 'longitude'] + 360

# Remove invalid 'year' entries and subset by 1900 to 2010 collection dates
data = data[ (data['year'] >= 1900) & (data['year'] <= 2010) ]      # Collection data between 1900 and 2010

# Rename columns
data = data.dropna( subset=['mass', 'latitude', 'longitude', 'year'] )  # Clear invalid entries

print( data.shape )

# Remove species that don't have at least 30 individuals, 20 years of collection, and 5 degrees of latitude
data = data.groupby('genus_species').filter( lambda species: ( (len(species['row_index']) >= 30) and
                                             ( (max(species['year']) - min(species['year']) ) >= 20) and
                                             ( (max(species['latitude']) - min(species['latitude']) ) >= 5) ))

print( data.shape )

# For each species ordered by mass, remove all individuals before the first 
# individual with an adult lifestage designation
print("Subsetting by first adult...")

def first_adult_mass(group):
    """ Helper function for groupby filter or apply
    Given a pandas DataFrameGroupBy object, find the mass of the
    lightest recorded 'adult' measurement.
    """
    adults = group[ np.in1d(group['underivedlifestage'], ['adult', 'ad', 'U-Ad.', 'U-Ad', 'Adult', 'Ad.']) ]
    if len(adults) == 0:
        return 0        # No adults found - delete all individuals
    else:
        return min(adults['mass'])   # Return lightest mass of any individual identified as an adult
    

data = data.groupby('genus_species').apply( lambda species: species[ species['mass'] >= first_adult_mass(species) ] )
print(data.shape)

# Subset again by number of individuals and ranges
data = data.groupby('genus_species').filter( lambda species: ( (len(species['row_index']) >= 30) and
                                             ( (max(species['year']) - min(species['year']) ) >= 20) and
                                             ( (max(species['latitude']) - min(species['latitude']) ) >= 5) ))

print( data.shape )

# Save results
data.to_csv(cleaned_filename)
