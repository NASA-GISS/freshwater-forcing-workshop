"""

This script is a guideline to redistribute the melt of each basin to the front of the ice shelves of the basins (instead of redistributing uniformly in front of the whole basin)
Step 1: Define the ice-shelf fronts (by regridding BedMachine to your native grid)
Step 2: Redistribute the melt at the ice-shelf fronts of given basins

BONUS AT THE END: Guideline on how to assign the minimum and maximum injection depth to points at the front and cut the maximum injection depth when local point at the front is shallower than it

Author: Clara Burgard

"""

import xarray as xr
import subprocess


inputpath = 'where_you_store_your_files_regarding_subshelfmelt'
inputpath_grid = 'where_your_gridfile_is' 
inputpath_BedMachine='where_your_BedMachine_file_has_been_downloaded' # download BedMachine v3 from here: https://nsidc.org/data/nsidc-0756/versions/3
outputpath = 'where_you_want_to_store_the_melt_forcing_file'

#### PREPARE MASKS FOR YOUR MODEL GRID

your_model_gridfile = inputpath_grid + 'gridfile.nc' # or gridfile.txt

# Regrid BedMachine ice shelf/grounded ice/ocean mask to your model grid
subprocess.run('cdo remapnn,'+your_model_gridfile+' -selname,mask '+inputpath_BedMachine+'BedMachineAntarctica-v3.nc '+inputpath+'BedMachine3_mask_onmodelgrid.nc', shell = True, executable="/bin/bash")
# Regrid the IMBIE basin mask to your model grid
subprocess.run('cdo remapnn,'+your_model_gridfile+' -selname,mask_region_extrapolated_to_ocean_2D '+inputpath+'IMBIE_mask_basins_extrapolated_to_ocean.nc '+inputpath+'IMBIE_mask_basins_extrapolated_to_ocean_onmodelgrid.nc',shell = True, executable="/bin/bash")

#### READ IN DATA
IMBIE_mask = xr.open_dataset(inputpath  + 'IMBIE_mask_basins_extrapolated_to_ocean_onmodelgrid.nc')
BedMachine_mask = xr.open_dataset(inputpath  + 'BedMachine3_mask_onmodelgrid.nc')
# set Lake Vostok and non-ice covered ground point to grounded ice
BedMachine_mask_clean = BedMachine_mask['mask'].where(BedMachine_mask['mask'] != 4,2).where(BedMachine_mask['mask']!= 1,2) 
# prescribed melt per basin, from zenodo
melt_file = xr.open_dataset(inputpath  + 'AQ_subshelf_melt.nc') 
# cell area of your model grid
cell_area = xr.open_dataset(inputpath + 'yourfile_withcellarea.nc')['cell_area'] 

#### IDENTIFY THE ICE SHELF FRONTS
# Create a mask discriminating between land, ocean and ice shelf
mask_0_1_2 = BedMachine_mask_clean.copy()
mask_0_1_2 = mask_0_1_2.where(BedMachine_mask_clean != 2, 400) # land
mask_0_1_2 = mask_0_1_2.where(BedMachine_mask_clean != 0, 0) # ocean
mask_0_1_2 = mask_0_1_2.where(BedMachine_mask_clean != 3, 200) # ice shelf

# set all ice shelves to 300
mask_front0 = mask_0_1_2.where((mask_0_1_2 == 0) | (mask_0_1_2 == 400), 300).copy()

mask_front = mask_front0.copy()
# check all directions and set points at border between ocean and ice shelf (300-0) to 500
# ATTENTION: replace lon and lat by the name of your horizontal coordinates 
mask_front = mask_front.where((mask_front0.shift(lon=-1)-mask_front0)!=300,500)
mask_front = mask_front.where((mask_front0.shift(lon=1)-mask_front0)!=300,500)
mask_front = mask_front.where((mask_front0.shift(lat=-1)-mask_front0)!=300,500)
mask_front = mask_front.where((mask_front0.shift(lat=1)-mask_front0)!=300,500)
# cut out all front points
mask_front = mask_front.where(mask_front==500)
# set the ice shelf number
mask_front = mask_front.where(mask_front!=500,IMBIE_mask['mask_region_extrapolated_to_ocean_2D'])


### DISTRIBUTE THE MELT
yearinsec = 86400.0 * 365.2422 # number of seconds in a year

melt_all = melt_file['baseline'] + melt_file['anomaly'] # melt in Gt/yr
melt_kg_per_s = melt_all / yearinsec *10**(12) # melt in kg per s

melt_front = mask_front.copy()
# redistribute the melt uniformly over the front of the ice shelves in the given basin
for id_basin in range(1,19):  
    cell_area_bb = cell_area.where(mask_front == id_basin)
    # melt in kg per s per m^2 
    melt_front = melt_front.where(mask_front != id_basin, melt_kg_per_s.sel(region=id_basin) / cell_area_bb.sum(['lon','lat'])).drop('region') 

melt_ds = xr.Dataset()
melt_ds['melt_to_be_injected'] = melt_front
melt_ds['melt_to_be_injected'].attrs['standard_name'] = 'Sub-shelf melt to be injected between minimal and maximal injection depth given in AQ_subshelf_melt.nc'
melt_ds['melt_to_be_injected'].attrs['units'] = 'kg s^-1 m^-2'
melt_ds.to_netcdf(outputpath + 'subshelfmelt_distributed_iceshelf_fronts.nc')

### BONUS: GUIDELINE HOW TO ASSIGN THE MINIMUM AND MAXIMUM INJECTION DEPTH TO POINTS AT THE FRONT AND CUT THE MAXIMUM INJECTION DEPTH WHEN LOCAL POINT AT THE FRONT IS SHALLOWER THAN IT
z_file = mask_front.copy()

zmin_all = melt_file['min_injection_depth']
zmax_all = melt_file['max_injection_depth']
bathy = xr.open_dataset('bathymetry_from_your_model.nc')

zmin_front = mask_front.copy()
zmax_front = mask_front.copy()
for id_basin in range(1,19): 
    zmin_front = zmin_front.where(mask_front != id_basin, zmin_all.sel(region=id_basin)).drop('region')
    zmax_front = zmax_front.where(mask_front != id_basin, zmax_all.sel(region=id_basin)).drop('region')
# cut where bathymetry is shallower
zmax_front = zmax_front.where(bathy >= zmax_front, bathy).where(np.isfinite(mask_front))

melt_ds['max_injection_depth_map'] = zmax_front.where(np.isfinite(zmax_front), 0)
melt_ds['min_injection_depth_map'] = zmin_front.where(np.isfinite(zmin_front), 0)
melt_ds.to_netcdf(outputpath + 'subshelfmelt_distributed_iceshelf_fronts_with_injectiondepthlimits.nc')

