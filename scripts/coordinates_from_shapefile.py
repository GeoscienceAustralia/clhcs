import geopandas as gpd
import pandas as pd

# shape_file = "/g/data/ge3/data/LHC/LHC_labels_K_sites.shp"

# shape_file = "additional/HV_dat_shape.shp"
# shape_file = "/home/sudipta/repos/uncover-ml/configs/data/geochem_sites.shp"
shape_file = "/home/sudipta/repos/clhc_sampling/LHC/LHC_labels_K_sites.shp"
shp_gdf = gpd.read_file(shape_file)

df = pd.DataFrame({'FID': shp_gdf.index, 'X_REF':shp_gdf.geometry.x, 'Y_REF': shp_gdf.geometry.y})
df.to_csv('T1_dat_K_WA.txt', index=False)
