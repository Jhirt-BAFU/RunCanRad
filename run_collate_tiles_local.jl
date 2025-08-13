################################################################################
# Program:    CanRad Tile Collation
# Purpose:    Collates individual radiation model sub-tiles into full 1km x 1km
#             tiles using spatial and model data, for further climate analysis.
#
# Usage:      julia this_script.jl <batch_tile>
#
# Arguments:  batch_tile - Batch name that corresponds to a tile list
#
# Author:       Clare Webster
# Maintainer:   Sandro Bischof
# Institution:  Federal Institute for Forest, Snow and Landscape Research WSL
#
# Adapted:
# 2025-05-14, Sandro — Modularized and documented for clarity and reuse.
#
# Dependencies: CanRad.jl, SpatialFileIO.jl, NCDatasets.jl, DelimitedFiles.jl,
#               collate_tiles_functions.jl, user-defined settings files.
################################################################################

using DelimitedFiles, Dates, Printf
using CanRad, SpatialFileIO, Formatting, NCDatasets

# Include custom functions
include("collate_tiles_functions_svf.jl")

# ---------------------- Utility Functions ----------------------
include("utilitiy_functions.jl")

"Parse the tile name string into spatial limits based on TILESIZE."
function extract_tile_limits(tile_name::String, tilesize::Int)
    xllcorner = parse(Int, split(tile_name, "_")[1])
    yllcorner = parse(Int, split(tile_name, "_")[2])
    return hcat(xllcorner, xllcorner + tilesize, yllcorner, yllcorner + tilesize)
end
 
"Process a single tile: verify input, collect data, and generate output file."
function process_tile(tile_name::String, input_folder::String, output_folder::String,
                      pts_file::String, par_in, pt_spacing::Int, tilesize::Int)

    indir  = joinpath(input_folder, tile_name)
    outdir = indir

    limits = extract_tile_limits(tile_name, tilesize)

    ptsx, ptsy, ptdx, _ = read_griddata_window(pts_file, limits, true, true)

    if isdir(indir)
        subtiles = readdir(indir)
        output_file = joinpath(output_folder, "Output_" * tile_name * ".nc")
        if !isfile(output_file) #&& length(subtiles) == 100
            collate2tilefile_fulltile(outdir, limits, tile_name, ptsx, ptsy,
                                      indir, output_folder, par_in, ptdx,
                                      pt_spacing, tilesize)
        end
    end
end

# ---------------------- Main Execution ----------------------

# Get batch argument from command line
batch_tile = "GRENZ" #ARGS[1]

# Constants
const TILESIZE = 1000
const PT_SPACING = 5

# Paths
input_path    = "E:/canrad_output/local/01_calculated_tiles"
model_path    = "E:/canrad_data/model"
settings_path = "E:/canrad_data/settings/local"
output_path   = "E:/canrad_output/local/02_collated_tiles"

# File and folder names
input_folder_name   = "output_" * batch_tile
input_folder        = joinpath(input_path, input_folder_name)
output_folder       = joinpath(output_path, input_folder_name)

pts_file   = joinpath(model_path, "ForestMask_NFI_TLM_Buffer25m_BorderClip_5m.tif")

# Load model settings
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
_, par_in = C2R_Settings(model_path)

include(joinpath(settings_path, "S2R_Settings_CHCW_local.jl"))
par_in_shi = S2R_Settings()

# Prepare execution
ensure_directory(output_folder)

start_time = now()

# Log start
println("CanRad - Collate tiles from ha to 1 km²")
println("Input path:                  $input_path")
println("Output path:                 $output_path")
println("Batch tile:                  $batch_tile")

# Get tile names directly from subdirectories
tile_dirs = filter(entry -> isdir(joinpath(input_folder, entry)), readdir(input_folder))
n_tiles = length(tile_dirs)

println("Number of tiles:             $n_tiles")
println("Started:                     " * format_datetime(start_time))

# Process tiles
for (n, tile_name) in enumerate(tile_dirs)
    n_percent = round(n/n_tiles * 100, digits=2)
    print("\rTile in progress:            $tile_name ($n/$n_tiles = $n_percent%)")
    process_tile(tile_name, input_folder, output_folder,
                 pts_file, par_in, PT_SPACING, TILESIZE)
end

# Final reporting
end_time = now()
elapsed = time_difference(start_time, end_time)
println("")
println("Ended:                       " * format_datetime(end_time))
println("Elapsed time:                $elapsed")
println("Summary:                     start=" * format_datetime(start_time) *
        " end=" * format_datetime(end_time) * " duration=" * elapsed)