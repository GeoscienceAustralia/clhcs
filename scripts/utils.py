import csv
from collections import defaultdict
from os import path
from pathlib import Path
from typing import Dict, List


def generate_key_val(gl, shorts):
    eight_char_name = Path(gl).stem[:8]
    if eight_char_name not in shorts:
        shorts[eight_char_name] += 1
    else:
        shorts[eight_char_name] += 1
    return gl, eight_char_name + str(shorts[eight_char_name])


def generate_key_vals_from_covariate_list_file(cov_list: str) -> Dict:
    geotifs_list = read_list_file(list_path=cov_list)
    geotifs = create_tif_dict(geotifs_list)
    return geotifs


def create_tif_dict(geotifs_list: List[str]) -> Dict:
    shorts = defaultdict(int)
    geotifs = {}
    for gl in geotifs_list:
        k, v = generate_key_val(gl, shorts)
        geotifs[k] = v
    return geotifs


def read_list_file(list_path: str) -> List[str]:
    files = []
    csvfile = path.abspath(list_path)
    with open(csvfile, 'r') as f:
        reader = csv.reader(f)
        tifs = list(reader)
        tifs = [f[0].strip() for f in tifs
                if (len(f) > 0 and f[0].strip() and
                    f[0].strip()[0] != '#')]
    for f in tifs:
        files.append(path.abspath(f))
    return files
