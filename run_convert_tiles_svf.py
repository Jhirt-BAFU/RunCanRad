"""
--------------------------------------------------------------------------------
Script:        run_convert_tiles_svf.py
Author:        Sandro Bischof
Date:          2025-05-26
Description:   Converts sky view fraction tiles from NetCDF format to GeoTIFF.
               Processes specified batch folders and writes georeferenced raster
               outputs using rasterio and xarray.
Usage:         python run_convert_tiles_svf.py [OPTION]... [BATCH_TILE]

Options:
  -h, --help      Show this help message and exit.

Example:
  python run_convert_tiles_svf.py batch_01
--------------------------------------------------------------------------------
"""

import glob
import os
import sys
import getopt
from datetime import datetime

import numpy as np
import rasterio as rio
import xarray as xr
from rasterio.transform import from_origin

# --- Constants: Input and Output base directories ---
INPUT_PATH = r'W:\CanRad_Euler\02_collated_tiles'
OUTPUT_PATH = r'W:\CanRad_Euler\03_converted_tiles'


def convert_tiles(input_file: str, output_folder: str) -> None:
    """
    Converts a NetCDF file to a GeoTIFF.
    Extracts the sky view fraction layer, applies the necessary transformations,
    and writes out a georeferenced raster.

    Parameters:
        input_file (str): Path to the input NetCDF file.
        output_folder (str): Directory where output GeoTIFF will be saved.
    """
    SCALE = 1
    FILL_VALUE = -99

    # Prepare the output filename from input filename parts
    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    parts = base_filename.split('_')
    if len(parts) < 3:
        raise ValueError(f"Filename '{base_filename}' does not follow the expected format with at least 3 parts.")
    filename_out = f"{parts[1]}_{parts[2]}"

    # Open NetCDF dataset with xarray
    with xr.open_dataset(input_file) as ds:
        svf_s = ds.svf_hemi_summer.to_numpy()
        svf_w = ds.svf_hemi_winter.to_numpy()
        loc_x = ds.Easting.to_numpy()
        loc_y = ds.Northing.to_numpy()

    # Flip and rotate array to match GeoTIFF orientation
    svf_s_rotated = np.rot90(np.flip(svf_s, axis=1), k=1)
    svf_w_rotated = np.rot90(np.flip(svf_w, axis=1), k=1)

    svf_s_scaled = np.nan_to_num(svf_s_rotated * SCALE, nan=FILL_VALUE)
    svf_w_scaled = np.nan_to_num(svf_w_rotated * SCALE, nan=FILL_VALUE)

    svf_s_clip = np.clip(svf_s_scaled, -32768, 32767).astype(np.int8)
    svf_w_clip = np.clip(svf_w_scaled, -32768, 32767).astype(np.int8)

    # Calculate raster metadata
    res_x = loc_x[1] - loc_x[0]
    res_y = loc_y[1] - loc_y[0]  # Usually negative, invert below
    x0 = loc_x[0] - (res_x / 2.0)
    y0 = loc_y[0] - (res_y / 2.0)

    # Construct output filepath
    out_file = os.path.join(output_folder, f"{filename_out}.tif")

    # Define affine transform (top-left origin)
    transform = from_origin(x0, y0, res_x, -res_y)

    # Set coordinate reference system (EPSG:2056)
    crs = rio.crs.CRS({"init": "epsg:2056"})

    # Write GeoTIFF raster
    with rio.open(
            out_file, 'w', driver='GTiff',
            height=svf_s_clip.shape[0],
            width=svf_s_clip.shape[1],
            count=2,
            dtype="int8",
            crs=crs,
            transform=transform,
            nodata=FILL_VALUE,
            compress='lzw',
            tiled=True
    ) as dst:
        dst.write(svf_s_clip, 1)
        dst.set_band_description(1, 'svf_hemi_summer')
        dst.write(svf_w_clip, 2)
        dst.set_band_description(2, 'svf_hemi_winter')


def format_datetime(dt: datetime) -> str:
    """Formats datetime object to a readable string 'YYYY-MM-DD HH:MM:SS'."""
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def time_difference(start_time: datetime, end_time: datetime) -> str:
    """Computes elapsed time between two datetime formatted as 'd-hh:mm:ss'."""
    delta = end_time - start_time
    total_seconds = int(delta.total_seconds())

    days = total_seconds // 86400
    hours = (total_seconds % 86400) // 3600
    minutes = (total_seconds % 3600) // 60
    seconds = total_seconds % 60

    return f"{days}-{hours:02}:{minutes:02}:{seconds:02}"


def print_usage() -> None:
    """Prints usage information."""
    print("Usage: run_convert_tiles_svf.py [OPTION]... [BATCH_TILE]")
    print("Options:")
    print("-h, --help                           Show this help message and exit.")


def start_task(batch_tile: str) -> None:
    """
    Starts the process of converting NetCDF tile files to GeoTIFF format for a specific batch of tiles.
    This function organises the directories for input and output, then processes each file within the
    specified batch directory and logs progress and execution time.

    :param batch_tile: A unique identifier for the batch of tiles to be processed. This ID determines
        the input directory, output directory, and the files to process.
    :return: None
    """
    start_time = datetime.now()

    output_folder = os.path.join(OUTPUT_PATH, batch_tile)
    os.makedirs(output_folder, exist_ok=True)

    input_folder = os.path.join(INPUT_PATH, f"output_{batch_tile}")

    print("CanRad - Convert tiles from NetCDF to GeoTIFF")
    print(f"{'Batch tile:':<30}{batch_tile}")
    print(f"{'Input path:':<30}{input_folder}")
    print(f"{'Output path:':<30}{output_folder}")
    print(f"{'Started:':<30}{format_datetime(start_time)}")

    files = glob.glob(os.path.join(input_folder, '*.nc'))
    i_tiles = len(files)
    print(f"{'Number of input files:':<30}{i_tiles}")

    # --- Key control ---
    for i, file in enumerate(files, start=1):
        i_percent = round(i / i_tiles * 100, 2)
        print(f"{'Tile in progress:':<30}{os.path.basename(file)}  ({i}/{i_tiles} = {i_percent}%)")
        convert_tiles(file, output_folder)

    end_time = datetime.now()
    elapsed = time_difference(start_time, end_time)

    print(f"{'Ended:':<30}{format_datetime(end_time)}")
    print(f"{'Elapsed time:':<30}{elapsed}")
    print(f"{'Summary:':<30}start={format_datetime(start_time)} end={format_datetime(end_time)} duration={elapsed}")


def main() -> None:
    """Parses command line arguments and launches the tile conversion process."""
    argument_list = sys.argv[1:]
    short_options = "h"
    long_options = ["help"]

    try:
        arguments, values = getopt.getopt(argument_list, short_options, long_options)

        batch = values[0] if values else None

        for current_arg, _ in arguments:
            if current_arg in ("-h", "--help"):
                print_usage()
                sys.exit(0)

        if batch is None:
            print("Error: No batch name specified.", file=sys.stderr)
            print_usage()
            sys.exit(1)

        start_task(batch)

    except getopt.GetoptError as err:
        print(f"Error: {err}", file=sys.stderr)
        print_usage()
        sys.exit(2)


# When Python executes a file, it sets a special variable, __name__, to "__main__" if the file is being run directly.
# If the file is imported as a module in another script, __name__ is set to the module's name instead.
# The if __name__ == "__main__": construct checks if the script is being run directly and only executes the main code in that case
if __name__ == "__main__":
    main()
