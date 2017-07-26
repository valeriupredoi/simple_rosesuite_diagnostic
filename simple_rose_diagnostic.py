import iris
import matplotlib
matplotlib.use('Agg')   # directs output to file rather than screen
import matplotlib.pyplot as plt
import subprocess
import time
import os
import numpy as np

"""
Simple time-series plotting script that takes in PP files
from a ROSE suite simulation and plots a variable of interest
defined by STASH code versus time

Valeriu Predoi, University of Reading, 2017
valeriu.predoi@ncas.ac.uk
"""
__author__ = "Valeriu Predoi <valeriu.predoi@ncas.ac.uk>"

############
# Functions
############
# get the area average
# climate model monitor function
def area_average(cube, region=None):
    """
    Calculates the area-weighted spatial mean of the passed-in cube. By default
    the global mean is calculated. If a region is defined then the mean is
    calculated over the corresponding lat/long subset of the input cube.

    :param iris.cube.Cube cube: The cube for which to calculate the spatial mean.
    :param list region: If specified, defines the geographical region over which
        to calculate the area-weighted mean. The region should be a list of
        lat/long coordinates in the order: [min-long, min-lat, max-long, max-lat]
    :returns: A scalar cube containing the spatial mean.
    """
    # Ensure that lat/long coordinates have bounds.
    for coord_name in ['latitude', 'longitude']:
        coord = cube.coord(coord_name)
        if not coord.has_bounds(): coord.guess_bounds()

    if region:
        minlon, minlat, maxlon, maxlat = region[:]
        cube = cube.intersection(latitude=(minlat, maxlat), longitude=(minlon, maxlon))

    # Calculate grid weights.
    grid_weights = iris.analysis.cartography.area_weights(cube)

    # Calculate the spatial mean.
    mean_cube = cube.collapsed(['latitude', 'longitude'], iris.analysis.MEAN,
        weights=grid_weights)

    return mean_cube

# remove unwanted coordinates
# clime model monitor function
def remove_coords(cube):
    """
    Remove non-essential coordinates from ``cube``. This should prevent a number
    of subsequent cube concatenation issues, plus reduce output file size. Note,
    however, that scalar coordinates, such as level number in the case of a
    vertical subset, are also removed.
    """
    removed_coord_names = []

    # First remove any auxiliary coordinates.
    for coord in cube.coords(dim_coords=False):
        removed_coord_names.append(coord.name())
        cube.remove_coord(coord)

    # Then remove any dimension coordinates not in the following list.
    dim_coords_to_keep = ['time', 'latitude', 'longitude']
    for coord in cube.coords(dim_coords=True):
        if coord.standard_name not in dim_coords_to_keep:
            removed_coord_names.append(coord.name())
            cube.remove_coord(coord)

    return removed_coord_names

######
# RUN
######
# Fixed and file-read parameters
months = ['jan','feb','mar','apr','may','jun','jul','aug','se','oct','nov','dec']#fixed
nr_months = ['01','02','03','04','05','06','07','08','09','10','11','12']#fixed
tz_root = '01T0000Z/'#fixed
ppfile = open('simple_rose_diagnostic_pars.txt','r')
# initiate empty lists to hold parameters
suite_names = []
rootdir = []
stash_codes = []
for line in ppfile:
    if line.strip().split('=')[0] == 'suite_names':
        sn = line.strip().split('=')[1].split(',')
        for snl in sn:
            suite_names.append(snl)
    if line.strip().split('=')[0] == 'rootdir':
        rootdir = line.strip().split('=')[1]
    if line.strip().split('=')[0] == 'stash_codes':
        stc = line.strip().split('=')[1].split(',')
        for stcl in stc:
            stash_codes.append(stcl)
    if line.strip().split('=')[0] == 'years':
        yrs = line.strip().split('=')[1]
        minyrs = int(yrs.split(',')[0])
        maxyrs = int(yrs.split(',')[1])
        years = [str(a) for a in range(minyrs,maxyrs+1)]

# loop over suites
for suite_name in suite_names:
    print('Analyzing suite %s' % suite_name)
    # atmos = a
    # TODO extend to other media
    fileDescriptor = '*' + suite_name + 'a.pm*'
    suite_rootdir = rootdir + '/' + 'u-' + suite_name + '/'
    # loop over stash codes
    for stash_code in stash_codes:
        # create NC files dir
        print('Analyzing STASH %s' % stash_code)
        ncFilesDir = 'NC_' + suite_name + '_' + stash_code +\
                      '_' + years[0] + '-' + years[-1] + '/'
        if os.path.isdir(ncFilesDir) is False:
            mkc = 'mkdir -p ' + ncFilesDir
            proc = subprocess.Popen(mkc, stdout=subprocess.PIPE, shell=True)
            (out, err) = proc.communicate()
        
        # Find PP files and load/convert to netCDF
        for yr in years:
            for mo in nr_months:
                basedir = suite_rootdir + yr + mo + tz_root
                lf = 'ls  ' + basedir + fileDescriptor
                proc = subprocess.Popen(lf, stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
                fp = out.strip()
                print('We will extract variable from file: %s' % fp)
                if os.path.isfile(fp):
                    stash_cons = iris.AttributeConstraint(STASH=stash_code)
                    t1 = time.time()
                    cb = iris.load(fp, constraints=stash_cons)[0]
                    frt = fp.split('/')[-1].strip('.pp')
                    locf = ncFilesDir + frt + '.nc'
                    iris.save(cb,locf)
                    t2 = time.time()
                    dt = t2-t1
                    print('netCDF file: %s' % locf)
                    print('Time it took to load variable and convert to netCDF: %.1f seconds' % dt)
        
        # Read each of the netCDF files and get global yearly mean
        cubeListAll = []
        skip_years = []
        for yr in years:
            cubeList_y = []
            for mo in months:
                ff = suite_name + 'a.pm' + yr + mo + '.nc'
                fp = ncFilesDir + ff
                print('Reading netCDF file: %s' % fp)
                if os.path.isfile(fp):
                    cb = iris.load_cube(fp)
                    cb_a = area_average(cb)
                    cubeList_y.append(np.mean(cb_a.data))
            # check if data is from exactly 12 months
            if len(cubeList_y) == 12:
                cubeListAll.append(np.mean(np.array(cubeList_y)))
            else:
                print('Data from year %s does not contain 12 months exactly, skipping year' % yr)
                skip_years.append(yr)

        #######
        # Plot
        #######
        # get rid of years we dont want to see
        new_years = [a for a in years if a not in skip_years]
        iyears = [int(b) for b in new_years]
        # time series
        fig = plt.figure()
        ax1 = fig.add_subplot(211)
        plotTitle = suite_name + ' ' + stash_code + ' time series and histogram'
        plt.suptitle(plotTitle)
        ax1.plot(iyears,cubeListAll)
        ax1.scatter(iyears,cubeListAll,marker='v')
        yLabel = stash_code + ' variable'
        plt.ylabel(yLabel)
        xl1 = min(iyears) - 1
        xl2 = max(iyears) + 1
        plt.xlim(xl1,xl2)
        plt.grid()
        # histogram
        ax2 = fig.add_subplot(212)
        ax2.hist(cubeListAll, cumulative=True)
        plt.grid()
        plt.ylabel('#')
        xLabel = stash_code + ' variable'
        plt.xlabel(xLabel)
        figName = suite_name + '_' + stash_code + '_time_hist.png'
        plt.savefig(figName)
        plt.close()
########
# ze end
########
