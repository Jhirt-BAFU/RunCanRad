################################################################################
# CanRad Batch Tile Processor
#
# Description:
# This script processes spatial tiles for the CanRad radiation model pipeline.
# It reads tile definitions, loads corresponding forest mask data, and computes
# radiative parameters (terrain and canopy-based) on a defined subtile grid.
# The process includes failure handling and progress tracking to support 
# batch execution in high-performance or distributed environments.
#
# Input:
#   ARGS[1] : Batch name (used to identify the tile list and output folder)
#   ARGS[2] : Start index of tile list
#   ARGS[3] : End index of tile list
#
# Output:
#   - Radiation model output files for each subtile
#   - Progress and failure logs for monitoring and reprocessing
#
# Author:       Clare Webster
# Maintainer:   Sandro Bischof
# Institution:  Federal Institution for Forest, Snow and Landscape Research WSL
#
# Adapted:    
# 2025-05-14, Sandro — Restructured and added English documentation and improved
#                      clarity.
#
# Requirements:
#   Julia packages: CanRad, SpatialFileIO, Formatting, NCDatasets, DelimitedFiles
#   Spatial inputs: Forest mask raster, tile list file, settings scripts
#
################################################################################
using DelimitedFiles
using CanRad, SpatialFileIO, Formatting, NCDatasets
using Profile, ProfileView, StatProfilerHTML

Profile.clear()

Profile.init(n = 10^7)
@profilehtml begin
# === Input arguments from command line ===
batch = "ZH" #ARGS[1]
job_start = 2    # parse(Int, ARGS[2])  # Index of the first job to process
job_end = 2 #parse(Int, ARGS[3])    # Index of the last job to process

# === Configuration ===
tile_size = 100        # Tile size in meters (must match the value in prep_cluster_input.jl)
sub_tile_size = 100     # Subtile size in meters (must evenly divide tile_size; used to reduce RAM usage)

# === Directory paths ===
#output_path = "C:/Users/Z70AJHI/project/RunCanRad/output"
#input_path = "W:/GIS/Projekte/Wald/Projekte/2025_LFI_CanRad/input"
#settings_path = "C:/Users/Z70AJHI/project/RunCanRad"

output_path = "C:/Users/joshu/Documents/BAFU/RunCanRad/output"
input_path = "C:/Users/joshu/Documents/BAFU/input"
settings_path = "C:/Users/joshu/Documents/BAFU/RunCanRad"

output_folder_name = "output_" * batch
output_folder = joinpath(output_path, output_folder_name)
output_folder_june = joinpath(output_path, output_folder_name * "_june")

progress_temp_folder = joinpath(output_path, "temp_progress_" * output_folder_name)
progress_folder = joinpath(output_path, "progress_" * output_folder_name)
fail_folder = joinpath(output_path, "failed_" * output_folder_name)

tiles_file = joinpath(settings_path, "tiles_" * batch * ".txt")
pts_file = joinpath(input_path, "ForestMask_NFI_TLM_Buffer25m_BorderClip_5m.tif")

# === Load user-defined settings ===
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
include(joinpath(settings_path, "S2R_Settings_CHCW_local.jl"))

dat_in, par_in = C2R_Settings(input_path)
par_in_shi = S2R_Settings()

@info dat_in["chmf"]


# === Main processing loop for each tile ===
for job_idx in job_start:job_end

    tile_name = string(readdlm(tiles_file)[job_idx])

    try
        # Skip tile if it has already been processed
        if isfile(joinpath(progress_folder, tile_name * ".txt")) || isfile(joinpath(progress_temp_folder, tile_name))
            println("Already completed: " * tile_name)
            continue
        end

        # === Determine spatial extent of current tile ===
        xllcorner = parse(Int, split(tile_name, "_")[1])
        yllcorner = parse(Int, split(tile_name, "_")[2])
        limits = hcat(xllcorner, xllcorner + tile_size, yllcorner, yllcorner + tile_size)

        # === Load forest mask data within current tile limits ===
        ptsx, ptsy, ptdx, _ = read_griddata_window(pts_file, limits, true, true)

        if sum(ptdx) > 0  # If the tile contains forest or terrain within Switzerland
            pts_all = hcat(ptsx[ptdx .> 0], ptsy[ptdx .> 0], ptdx[ptdx .> 0])

            # === Define subtiles within current tile ===
            limx = Int(floor(limits[1] / sub_tile_size)) * sub_tile_size : sub_tile_size : Int(floor(limits[2] / sub_tile_size)) * sub_tile_size
            limy = Int(floor(limits[3] / sub_tile_size)) * sub_tile_size : sub_tile_size : Int(floor(limits[4] / sub_tile_size)) * sub_tile_size

            outdir = joinpath(output_folder, tile_name)
            outdir_june = joinpath(output_folder_june, tile_name)

            dim = Int(tile_size / sub_tile_size)

            # === Process all subtiles within the current tile ===
            for x in 1:dim, y in 1:dim
                idx = (limx[x] .<= pts_all[:, 1] .< limx[x+1]) .& (limy[y] .<= pts_all[:, 2] .< limy[y+1])

                if sum(idx) > 0
                    pts = pts_all[idx, :]
                    tilestr = string(limx[x]) * "_" * string(limy[y])
                    taskID = string(job_idx) * "_" * tilestr

                    # Check if output already exists for this subtile
                    if !check_output(outdir, pts[1:10, :], true, taskID)

                        # If all points are terrain (value 2), use ter2rad!
                        if sum(pts[:, 3]) == size(pts, 1) * 2
                            @info "ter2rad! " * taskID
                            ter2rad!(pts, dat_in, par_in, outdir, taskID)
                        else
                            @info "chm2rad! " * taskID
                            chm2rad!(pts, dat_in, par_in, outdir, taskID)
                        end

                        # Remove temporary radiation file after computation
                        shi_file = joinpath(outdir, tilestr, "SHIs_" * tilestr * ".nc")
                        rm(shi_file, force = true)

                    else
                        println("Already completed: " * taskID)
                    end
                end
            end

            # === Mark tile as temporarily completed ===
            mkpath(progress_temp_folder)
            writedlm(joinpath(progress_temp_folder, tile_name), NaN)

            # Remove failure log if it exists
            if isfile(joinpath(fail_folder, tile_name * ".txt"))
                rm(joinpath(fail_folder, tile_name * ".txt"), force = true)
            end

            println("Tile completed: " * tile_name)

        else
            # Tile is entirely outside Switzerland – no processing needed
            mkpath(progress_folder)
            writedlm(joinpath(progress_folder, tile_name * ".txt"), NaN)
        end

    catch e
        # === Error handling for failed tile ===
        mkpath(fail_folder)
        writedlm(joinpath(fail_folder, tile_name * ".txt"), NaN)
        println("Failed to process tile: " * tile_name)
        println("Caught exception: ", e)
        println("Stack trace:")
        for (i, frame) in enumerate(catch_backtrace())
            println("[$i] $frame")
        end
    end
end
end #profile end

ProfileView.view()
statprofilehtml()