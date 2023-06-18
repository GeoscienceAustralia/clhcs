from typing import Dict
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from scripts.utils import generate_key_vals_from_covariate_list_file


def intersect_and_sample_shp(shp: Path, geotifs: Dict[str, str], intersect_parallel: bool = False):
    print("====================================\n", f"intersecting {shp.as_posix()}")
    pts = gpd.read_file(shp)
    coords = np.array([(p.x, p.y) for p in pts.geometry])
    geom = pd.DataFrame(coords, columns=geom_cols, index=pts.index)
    pts = pts.merge(geom, left_index=True, right_index=True)

    src = rasterio.open(list(geotifs.keys())[0])
    extent = src.bounds
    cell_size = src.res[0]

    # Get the cell number from the pixel coordinates
    pts["cell_no"] = [int(((c[0] - extent[0]) + (c[1] - extent[1]) * src.width)/cell_size) for c in coords]

    for i, (k, v) in enumerate(geotifs.items()):
        print(f"adding {i}/{len(geotifs)}: {k} to output dataframe")
        pts[v] = __extract_raster_points(k, coords)

    return pts


def __extract_raster_points(raster: str, coordinates: np.ndarray):
    with rasterio.open(Path(raster)) as src:
        print(f"---- intersecting {raster}--------")
        return [x[0] for x in src.sample(coordinates)]


def intersect_sample_and_clean(shp, write_dropped: bool = False, intersect_parallel: bool = False):
    print("=========here here ==========")
    df2 = intersect_and_sample_shp(shp, geotifs, intersect_parallel=intersect_parallel)
    output_dir.mkdir(exist_ok=True, parents=True)
    output_dir.joinpath('clean')
    out_shp = output_dir.joinpath(shp.name)
    df3 = df2[list(geotifs.values())]
    finite_idices = ((np.isfinite(df3)).sum(axis=1) == len(geotifs)) & \
                    ((np.abs(df3) < 1e10).sum(axis=1) == len(geotifs))
    df4 = df2.loc[finite_idices, :]
    df5 = df2.loc[~finite_idices, :]
    try:
        if df5.shape[0]:
            df4.to_file(out_shp.parent.joinpath(out_shp.stem + '_cleaned.shp'))
            print(f"Wrote clean shapefile {out_shp.parent.joinpath(out_shp.stem + '_cleaned.shp')}")
            if write_dropped:
                df5.to_file(out_shp.parent.joinpath(out_shp.stem + '_cleaned_dropped.shp'))
        else:
            d = {'X_REF': df2.geometry.x, 'Y_REF': df2.geometry.y}
            d.update({v: df2[v] for v in geotifs.values()})
            df = pd.DataFrame(d)
            df.iloc[:, :].to_csv("intersected_covs.txt", header=True, index=False)
            df2.to_file(out_shp.as_posix())
            print(f"saved intersected shapefile at {out_shp.as_posix()}")
            print(f"No points dropped and there for _cleaned.shp file is not createed'")
            print(f"No points dropped and there for _cleaned_dropped.shp file is not created")
    except Exception as e:
        print(e)
        print(f"Check this shapefile {shp}")


if __name__ == '__main__':
    # local
    output_dir = Path('out_resampled')
    shape = Path("/home/sudipta/repos/uncover-ml/configs/data/geochem_sites.shp")
    geom_cols = ['POINT_X', 'POINT_Y']

    covariastes_list = "/home/sudipta/repos/uncover-ml/configs/data/sirsam/covariates_list.txt"
    geotifs = generate_key_vals_from_covariate_list_file(covariastes_list)

    # check all required files are available on disc
    for k, v in geotifs.items():
        print(f"checking if {k} exists")
        assert Path(k).exists()

    print("Add under the 'intersected_features:' section in yaml")
    print('='*100)
    for k, v in geotifs.items():
        print(f"\"{Path(k).name}\": \"{v}\"")

    print('='*100)
    intersect_sample_and_clean(shape, write_dropped=True, intersect_parallel=True)
