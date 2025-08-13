"""
--------------------------------------------------------------------------------
Script:        run_convert_tiles_svf.py
Author:        Sandro Bischof
Date:          2025-05-26
Description:   Merges sky view fraction GeoTIFF tiles into a seamless mosaic.
               Processes tiles from a specified batch directory and generates
               a georeferenced output using rasterio.
Usage:         python run_convert_tiles_svf.py [OPTIONS] BATCH_TILE

Options:
  -h, --help              Show this help message and exit.
  -i, --input=DIR         Specify custom input directory.
  -o, --output=DIR        Specify custom output directory.

Example:
  python run_convert_tiles_svf.py -i W:/input_path -o W:/output_path batch_01
--------------------------------------------------------------------------------
"""

import getopt
import glob
import os
import sys
from datetime import datetime
from typing import List

import rasterio as rio
from rasterio.io import DatasetReader
from rasterio.merge import merge

# --- Constants: Input and Output base directories ---
DEFAULT_INPUT_PATH = r'W:\CanRad_Euler\03_converted_tiles'
DEFAULT_OUTPUT_PATH = r'W:\CanRad_Euler\04_stuck_tiles'


def stick_tiles_together(all_tiles: List[DatasetReader], meta: dict, descriptions, batch_tile: str, output_folder: str) -> None:
    """
    Merges a list of raster tiles into a single mosaic and writes it to disk as a GeoTIFF.

    Parameters:
    - all_tiles: List of open rasterio DatasetReader objects representing tiles to be merged.
    - meta: Metadata dictionary (copied from one of the input tiles) for configuring output.
    - batch_tile: Identifier used to name the output file.
    - output_folder: Path to the folder where the output GeoTIFF should be saved.
    """
    # Merge the tiles into a mosaic
    mosaic, out_transform = merge(all_tiles)

    # Set target CRS (EPSG:2056)
    crs = rio.crs.CRS({"init": "epsg:2056"})

    # Update metadata to match the new mosaic
    meta.update({
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "transform": out_transform,
        "crs": crs,
        "compress": "lzw"
    })

    # Construct an output file path
    output_file = os.path.join(output_folder, f"map_{batch_tile}.tif")

    # Write mosaic to a file using a context manager
    # with rio.open(output_file, "w", **meta) as dest:
    #     dest.write(mosaic)
    with rio.open(output_file, "w", **meta) as dst:
        dst.write(mosaic)
        if descriptions and len(descriptions) == mosaic.shape[0]:
            for i, desc in enumerate(descriptions, start=1):
                if desc:
                    dst.set_band_description(i, desc)

def format_datetime(dt: datetime) -> str:
    """Formats datetime object to a readable string 'YYYY-MM-DD HH:MM:SS'."""
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def time_difference(start_time: datetime, end_time: datetime) -> str:
    """Computes elapsed time between two datetimes formatted as 'd-hh:mm:ss'."""
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
    print("-i, --input                          Set input directory.")
    print("-o, --output                         Set output directory.")


def start_task(batch_tile: str, input_path: str, output_path: str) -> None:
    input_path = input_path or DEFAULT_INPUT_PATH
    output_path = output_path or DEFAULT_OUTPUT_PATH
    start_time = datetime.now()

    input_folder = os.path.join(input_path, batch_tile)
    output_folder = output_path
    os.makedirs(output_folder, exist_ok=True)

    print("CanRad - Stick GeoTIFFs together")
    print(f"{'Batch tile:':<30}{batch_tile}")
    print(f"{'Input path:':<30}{input_folder}")
    print(f"{'Output path:':<30}{output_folder}")
    print(f"{'Started:':<30}{format_datetime(start_time)}")

    files = glob.glob(os.path.join(input_folder, "*.tif"))
    print(f"{'Number of input files:':<30}{len(files)}")

    if not files:
        print("Error: No input files found. Aborting task.", file=sys.stderr)
        return

    all_tiles = []
    for i, file in enumerate(files, start=1):
        print(f"{'Tile in progress:':<30}{os.path.basename(file)}  ({i}/{len(files)})")
        tile = rio.open(file)
        all_tiles.append(tile)

    meta = all_tiles[-1].meta.copy()
    descriptions = all_tiles[-1].descriptions
    stick_tiles_together(all_tiles, meta, descriptions, batch_tile, output_folder)

    for tile in all_tiles:
        tile.close()

    end_time = datetime.now()
    elapsed = time_difference(start_time, end_time)
    print(f"{'Ended:':<30}{format_datetime(end_time)}")
    print(f"{'Elapsed time:':<30}{elapsed}")
    print(f"{'Summary:':<30}start={format_datetime(start_time)} end={format_datetime(end_time)} duration={elapsed}")


def main() -> None:
    """Parses command line arguments and launches the tile conversion process."""
    short_options = "hi:o:"
    long_options = ["help", "input=", "output="]

    input_path = None
    output_path = None

    try:
        arguments, remaining_args = getopt.getopt(sys.argv[1:], short_options, long_options)

        for current_arg, current_value in arguments:
            if current_arg in ("-h", "--help"):
                print_usage()
                sys.exit(0)
            elif current_arg in ("-i", "--input"):
                input_path = current_value
            elif current_arg in ("-o", "--output"):
                output_path = current_value

        if not remaining_args:
            print("Error: No batch tile specified.", file=sys.stderr)
            print_usage()
            sys.exit(1)

        batch_tile = remaining_args[0]
        start_task(batch_tile, input_path, output_path)

    except getopt.GetoptError as err:
        print(f"Error: {err}", file=sys.stderr)
        print_usage()
        sys.exit(2)


# When Python executes a file, it sets a special variable, __name__, to "__main__" if the file is being run directly.
# If the file is imported as a module in another script, __name__ is set to the module's name instead.
# The if __name__ == "__main__": construct checks if the script is being run directly and only executes the main code in that case
if __name__ == "__main__":
    main()
