from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from scripts.utils import create_tif_dict, generate_key_vals_from_covariate_list_file, read_list_file


def raster_points(tifs) -> pd.DataFrame:
    """
    Function to convert a list of raster files into a pandas dataframe
    Args:
        tifs (list): List of raster files
    Returns:
        covariates (pd.DataFrame): Dataframe containing the coordinates and values of each pixel in the raster
    """
    geotifs_dict = create_tif_dict(tifs)
    with rio.Env():
        # Get the coordinates of the top left and bottom right corners of the raster
        with rio.open(tifs[0]) as src:
            crs = src.crs
            xmin, ymax = np.around(src.xy(0.00, 0.00), 9)
            xmax, ymin = np.around(src.xy(src.height - 1, src.width - 1), 9)
            # Create a grid of x and y coordinates
            x = np.linspace(xmin, xmax, src.width)
            y = np.linspace(ymax, ymin, src.height)
            xs, ys = np.meshgrid(x, y)
        # Create a dictionary to store the coordinates and values of each pixel in the raster
        data = {
            "X_REF": pd.Series(xs.ravel()),
            "Y_REF": pd.Series(ys.ravel())
        }
        covariate_cols = []
        # Iterate through each raster file
        for tif in tifs:
            with rio.open(tif) as src:
                # Read the data from the raster
                zs = src.read(1, masked=True)
                # Apply NoData mask
                zs_data = zs.data
                zs_data[zs.mask] = np.nan
                # xs, ys, zs = xs[mask], ys[mask], zs[mask]
                data[geotifs_dict[tif]] = pd.Series(zs_data.flatten())
                # Add the column name to the list of covariates
                covariate_cols.append(geotifs_dict[tif])

    # the clhs R code cannot handle nans, so we filter them out
    # Create a pandas dataframe from the dictionary
    covariates = pd.DataFrame(data=data)
    print(covariates.head())
    num_covaraites = covariates[covariate_cols].shape[0]
    return covariates


#
# rasters = [
#     # "/home/sudipta/repos/uncover-ml/configs/data/LATITUDE_GRID1.tif",
#     # "/home/sudipta/repos/uncover-ml/configs/data/LONGITUDE_GRID1.tif",
#     # "/home/sudipta/repos/uncover-ml/configs/data/sirsam/k_15v5.tif",
#     # "/home/sudipta/repos/uncover-ml/configs/data/sirsam/outcrop_dis2.tif",
#     # "/home/sudipta/repos/uncover-ml/configs/data/sirsam/PM_Gravity.tif",
#     # "/home/sudipta/repos/uncover-ml/configs/data/sirsam/Th_v1.tif",
#     "/home/sudipta/repos/clhc_sampling/LHC/crop_G_Dose.tif",
#     "/home/sudipta/repos/clhc_sampling/LHC/crop_G_K.tif",
#     "/home/sudipta/repos/clhc_sampling/LHC/crop_T_DEM_S.tif",
#     "/home/sudipta/repos/clhc_sampling/LHC/crop_T_Saga_Wet_H.tif"
# ]

# covariates_list = "/home/sudipta/repos/uncover-ml/configs/data/sirsam/covariates_list.txt"
covariates_list = "/home/sudipta/repos/clhc_sampling/LHC/lhc_covs.txt"
rasters = read_list_file(covariates_list)
print(rasters)

geotifs_dict = generate_key_vals_from_covariate_list_file(covariates_list)
print(geotifs_dict)
# import sys; sys.exit()
raster_df = raster_points(rasters)
raster_df.iloc[:2000, :].to_csv("LHC/T1_covs_part.txt", header=True, index=False)
raster_df.iloc[:, :].to_csv("LHC/T1_covs_full.txt", header=True, index=False)
