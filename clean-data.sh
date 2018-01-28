#!/bin/bash

# Download and process data files
# Result should be one csv file each for VertNet and NOAA data
#   NOTE: download and extraction of VertNet files requires at least 70 Gb free (temporarily)

if [ ! -f VertNet.csv ]; then
    echo "VertNet data not found. Downloading..."

    mkdir vertnet-temp
    cd vertnet-temp

    # Download VertNet datasets
    wget https://de.iplantcollaborative.org/anon-files//iplant/home/shared/commons_repo/curated/Vertnet_Amphibia_Sep2016/VertNet_Amphibia_Sept2016.zip
    wget https://de.iplantcollaborative.org/anon-files//iplant/home/shared/commons_repo/curated/Vertnet_Aves_Sep2016/VertNet_Aves_Sept2016.zip
    wget https://de.iplantcollaborative.org/anon-files//iplant/home/shared/commons_repo/curated/Vertnet_Mammalia_Sep2016/VertNet_Mammalia_Sept2016.zip
    wget https://de.iplantcollaborative.org/anon-files//iplant/home/shared/commons_repo/curated/Vertnet_Reptilia_Sep2016/VertNet_Reptilia_Sept2016.zip

    unzip *.zip
    rm -r _MACOSX *.xml
    rm *.zip
    cd ..

    # Combine datasets and extract desired columns
    echo "Combining and cleaning VertNet data..."
    python3 clean.py

    # Do you need these files?
    mv all_vertnet.csv ..   
    mv clean_taxonomy.csv ..

    mv VertNet.csv ..
    cd ..

    rm -r vertnet-temp
fi


if [ ! -f NOAA.csv ]; then
    echo "Temperature data not found. Creating...."

    if [ ! -f air.mon.mean.v301.nc ]; then
        echo "NOAA dataset not found. Downloading..."
        wget ftp://ftp.cdc.noaa.gov/Datasets/udel.airt.precip/air.mon.mean.v301.nc
    fi
    
fi  
END
