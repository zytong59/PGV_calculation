#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install geopandas')
get_ipython().system('pip install pandas')


# In[ ]:


import geopandas as gpd
import rasterio
from rasterio.mask import mask
from shapely.geometry import mapping
import numpy as np
import os
import pandas as pd



## Input files
# Set workspace
input_workspace = r'E:/20241202GE_calculation_ThuCloud/inputfile_GE/'
output_workspace = r'E:/20241202GE_calculation_ThuCloud/outputfile_GE/'

# Set AHP Results for effective GE
ahp_results = {
    'ahp_1': 0.11143422,
    'ahp_2': 0.22316116,
    'ahp_3': 0.265157492,
    'ahp_4': 0.084136101,
    'ahp_5': 0.043952652,
    'ahp_6': 0.162014177,
    'ahp_7': 0.110144197,
    # Set non-working ratio for GE_2
    'non_working': 0.3150685
}



## Zonal for PK with RCRS_buffer
# Read shpfile
shp_file_path = input_workspace + "RCRS_900m_buffer.shp"
gdf = gpd.read_file(shp_file_path)
# Read tiffile
tif_file_path = input_workspace + "PK_BeijingVeg2020.tif"
with rasterio.open(tif_file_path) as src:
    for i in range(1, 8):
        for index, row in gdf.iterrows():
            # Get Polygon from the shp file
            geom = row['geometry']
            geoms = [mapping(geom)]
            # Read the number of pixels from the TIF file
            out_image, out_transform = mask(src, geoms, crop=True)
            # Calculate the number of pixels with the value of i
            veg_count = np.count_nonzero(out_image == i)
            # Write in shp
            gdf.at[index, 'PK_' + str(i)] = veg_count
        print("     Complete the processing PK classification of " + str(i))
# Save new shpfile
output_shp_file_path = output_workspace + "zonal_PK.shp"
gdf.to_file(output_shp_file_path)
print("New SHPfile has been saved to：" + output_shp_file_path)
#Save the attribute table as a CSV file
gdf2 = gpd.read_file(output_shp_file_path)
csv_file_path = output_workspace + "zonal_PK.csv"
gdf2.to_csv(csv_file_path,index=False)
print("New CSVfile has been saved to：" + csv_file_path)


## Zonal for CI with RCRS_buffer
# Read shpfile
shp_file_path = input_workspace + "RCRS_900m_buffer.shp"
gdf = gpd.read_file(shp_file_path)
# Read tiffile
tif_file_path = input_workspace + "CI_BeijingVeg2020.tif"
with rasterio.open(tif_file_path) as src:
    for i in range(1, 8):
        for index, row in gdf.iterrows():
            # Get Polygon from the shp file
            geom = row['geometry']
            geoms = [mapping(geom)]
            # Read the number of pixels from the TIF file
            out_image, out_transform = mask(src, geoms, crop=True)
            # Calculate the number of pixels with the value of i
            veg_count = np.count_nonzero(out_image == i)
            # Write in shp
            gdf.at[index, 'CI_' + str(i)] = veg_count
        print("     Complete the processing CI classification of " + str(i))
# Save new shpfile
output_shp_file_path = output_workspace + "zonal_CI.shp"
gdf.to_file(output_shp_file_path)
print("New SHPfile has been saved to：" + output_shp_file_path)
#Save the attribute table as a CSV file
gdf2 = gpd.read_file(output_shp_file_path)
csv_file_path = output_workspace + "zonal_CI.csv"
gdf2.to_csv(csv_file_path,index=False)
print("New CSVfile has been saved to：" + csv_file_path)


## Zonal for RC with RC_boundary
# Read shpfile
shp_file_path = input_workspace + "RC_boundary.shp"
gdf = gpd.read_file(shp_file_path)
# Read tiffile
tif_file_path = input_workspace + "RC_BeijingVeg2020.tif"
with rasterio.open(tif_file_path) as src:
    for i in range(1, 8):
        for index, row in gdf.iterrows():
            # Get Polygon from the shp file
            geom = row['geometry']
            geoms = [mapping(geom)]
            # Read the number of pixels from the TIF file
            out_image, out_transform = mask(src, geoms, crop=True)
            # Calculate the number of pixels with the value of i
            veg_count = np.count_nonzero(out_image == i)
            # Write in shp
            gdf.at[index, 'RC_' + str(i)] = veg_count
        print("     Complete the processing RC classification of " + str(i))
# Save new shpfile
output_shp_file_path = output_workspace + "zonal_RC.shp"
gdf.to_file(output_shp_file_path)
print("New SHPfile has been saved to：" + output_shp_file_path)
#Save the attribute table as a CSV file
gdf2 = gpd.read_file(output_shp_file_path)
csv_file_path = output_workspace + "zonal_RC.csv"
gdf2.to_csv(csv_file_path,index=False)
print("New CSVfile has been saved to：" + csv_file_path)


## Zonal for GE0 with RCRS_buffer
# Read shpfile
shp_file_path = input_workspace + "RCRS_900m_buffer.shp"
gdf = gpd.read_file(shp_file_path)
# Read tiffile
tif_file_path = input_workspace + "BeijingVeg2020_1.tif"
with rasterio.open(tif_file_path) as src:
    for i in range(1, 2):
        for index, row in gdf.iterrows():
            # Get Polygon from the shp file
            geom = row['geometry']
            geoms = [mapping(geom)]
            # Read the number of pixels from the TIF file
            out_image, out_transform = mask(src, geoms, crop=True)
            # Calculate the number of pixels with the value of i
            veg_count = np.count_nonzero(out_image == i)
            # Write in shp
            gdf.at[index, 'All_' + str(i)] = veg_count
        print("     Complete the processing all greensapce classification of " + str(i))
# Save new shpfile
output_shp_file_path = output_workspace + "zonal_allbuffer.shp"
gdf.to_file(output_shp_file_path)
print("New SHPfile has been saved to：" + output_shp_file_path)
#Save the attribute table as a CSV file
gdf2 = gpd.read_file(output_shp_file_path)
csv_file_path = output_workspace + "zonal_allbuffer.csv"
gdf2.to_csv(csv_file_path,index=False)
print("New CSVfile has been saved to：" + csv_file_path)



## Calculate for GE
# Read csvfile
rcrs_origin = pd.read_csv(input_workspace + 'RCRS_origin.csv')
zonal_allbuffer = pd.read_csv(output_workspace + 'zonal_allbuffer.csv')
zonal_pk = pd.read_csv(output_workspace + 'zonal_PK.csv')
zonal_rc = pd.read_csv(output_workspace + 'zonal_RC.csv')
zonal_ci = pd.read_csv(output_workspace + 'zonal_CI.csv')
print("Read csv files already")

# Calculate GE_0
# Connect GE_calculation_output and zonal_allbuffer according to RCRS_id
ge_calculation_output = pd.merge(rcrs_origin, zonal_allbuffer[['RCRS_id', 'All_1']], on='RCRS_id', how='left')
ge_calculation_output['GE_0'] = (ge_calculation_output['All_1']) * 0.000005163162
print("     GE_0 calculation finished")

# Calculate GE_1
# Connect GE_calculation_output and zonal_rc &zonal_rc according to RCRS_id
ge_calculation_output = pd.merge(ge_calculation_output, zonal_pk[['RCRS_id', 'PK_1', 'PK_2', 'PK_3', 'PK_4', 'PK_5', 'PK_6', 'PK_7']], on='RCRS_id', how='left')
ge_calculation_output = pd.merge(ge_calculation_output, zonal_rc[['OBJECTID', 'RC_1', 'RC_2', 'RC_3', 'RC_4', 'RC_5', 'RC_6', 'RC_7']], on='OBJECTID', how='left')
# Fill 0
ge_calculation_output.fillna(0, inplace=True)
# GE_1 formula
ge_calculation_output['GE_1'] = (ge_calculation_output['PK_1']
                                 + ge_calculation_output['PK_2']
                                 + ge_calculation_output['PK_3']
                                 + ge_calculation_output['PK_4']
                                 + ge_calculation_output['PK_5']
                                 + ge_calculation_output['PK_6']
                                 + ge_calculation_output['PK_7']
                                 + ge_calculation_output['RC_1']
                                 + ge_calculation_output['RC_2']
                                 + ge_calculation_output['RC_3']
                                 + ge_calculation_output['RC_4']
                                 + ge_calculation_output['RC_5']
                                 + ge_calculation_output['RC_6']
                                 + ge_calculation_output['RC_7']) * 0.000005163162
print("     GE_1 calculation finished")

# Calculate EGE_1
# EGE_1 formula
ge_calculation_output['EGE_1'] = (ge_calculation_output['PK_1'] * ahp_results['ahp_1']
                                 + ge_calculation_output['PK_2'] * ahp_results['ahp_2']
                                 + ge_calculation_output['PK_3'] * ahp_results['ahp_3']
                                 + ge_calculation_output['PK_4'] * ahp_results['ahp_4']
                                 + ge_calculation_output['PK_5'] * ahp_results['ahp_5']
                                 + ge_calculation_output['PK_6'] * ahp_results['ahp_6']
                                 + ge_calculation_output['PK_7'] * ahp_results['ahp_7']
                                 + ge_calculation_output['RC_1'] * ahp_results['ahp_1']
                                 + ge_calculation_output['RC_2'] * ahp_results['ahp_2']
                                 + ge_calculation_output['RC_3'] * ahp_results['ahp_3']
                                 + ge_calculation_output['RC_4'] * ahp_results['ahp_4']
                                 + ge_calculation_output['RC_5'] * ahp_results['ahp_5']
                                 + ge_calculation_output['RC_6'] * ahp_results['ahp_6']
                                 + ge_calculation_output['RC_7'] * ahp_results['ahp_7']) * 0.000005163162
print("     EGE_1 calculation finished")

# Calculate GE_2
# Connect GE_calculation_output and zonal_ci according to RCRS_id
ge_calculation_output = pd.merge(ge_calculation_output, zonal_ci[['RCRS_id', 'CI_1', 'CI_2', 'CI_3', 'CI_4', 'CI_5', 'CI_6', 'CI_7']], on='RCRS_id', how='left')
# Fill 0
ge_calculation_output.fillna(0, inplace=True)
# GE_2 formula
ge_calculation_output['GE_2'] = (ge_calculation_output['CI_1'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_2'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_3'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_4'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_5'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_6'] * ahp_results['non_working']
                                 + ge_calculation_output['CI_7'] * ahp_results['non_working']) * 0.000005163162 + ge_calculation_output['GE_1']
print("     GE_2 calculation finished")

# Calculate EGE_2
# EGE_2 formula
ge_calculation_output['EGE_2'] = (ge_calculation_output['CI_1'] * ahp_results['non_working'] * ahp_results['ahp_1']
                                 + ge_calculation_output['CI_2'] * ahp_results['non_working'] * ahp_results['ahp_2']
                                 + ge_calculation_output['CI_3'] * ahp_results['non_working'] * ahp_results['ahp_3']
                                 + ge_calculation_output['CI_4'] * ahp_results['non_working'] * ahp_results['ahp_4']
                                 + ge_calculation_output['CI_5'] * ahp_results['non_working'] * ahp_results['ahp_5']
                                 + ge_calculation_output['CI_6'] * ahp_results['non_working'] * ahp_results['ahp_6']
                                 + ge_calculation_output['CI_7'] * ahp_results['non_working'] * ahp_results['ahp_7']) * 0.000005163162 + ge_calculation_output['EGE_1']
print("     EGE_2 calculation finished")

# Save results to GE_calculation.csv
ge_calculation_output.to_csv(output_workspace + 'GE_calculation_output.csv', index=False)
print("All calculations finished")

# Summary of GE numerical results
# Select the specified column
selected_columns = ['RCRS_id', 'GE_0', 'GE_1', 'GE_2', 'EGE_1', 'EGE_2']
# Create a new DataFrame that contains only the specified columns
ge_value_final = ge_calculation_output[selected_columns]
# Save results to GE_FINAL.csv
ge_value_final.to_csv(output_workspace + 'GE_value_RESULT.csv', index=False)
print("Output files finished")

