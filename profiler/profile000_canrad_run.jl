### profile canrad at one camera location
using Profile, PProf, BenchmarkTools
###

using DelimitedFiles
using CanRad, SpatialFileIO, Formatting, NCDatasets
using Base.Threads

# === Directory paths ===
output_path = "/mnt/output"
input_path = "/mnt/input"
settings_path = ".."

output_folder_name = "profiling"
output_folder = joinpath(output_path, output_folder_name)

progress_temp_folder = joinpath(output_path, "temp_progress_" * output_folder_name)
progress_folder = joinpath(output_path, "progress_" * output_folder_name)
fail_folder = joinpath(output_path, "failed_" * output_folder_name)

#tiles_file = joinpath(settings_path, "tiles_" * batch * ".txt")
pts_file = joinpath(input_path, "ForestMask_NFI_TLM_Buffer25m_BorderClip_5m.tif")

# === Load user-defined settings ===
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
include(joinpath(settings_path, "S2R_Settings_CHCW_local.jl"))

dat_in, par_in = C2R_Settings(input_path)
par_in_shi = S2R_Settings()

# === Input arguments from command line ===
#batch = "ZH" #ARGS[1]
#job_start = 16 # parse(Int, ARGS[2])  # Index of the first job to process
#job_end = 16 #parse(Int, ARGS[3])    # Index of the last job to process

# === Configuration ===

#par_in["make_geotiff"] = true
#tile_size = 10        # Tile size in meters
#sub_tile_size = 10     # Subtile size in meters (must evenly divide tile_size; used to reduce RAM usage)

par_in["make_geotiff"] = false
tile_size = 5        # Tile size in meters
sub_tile_size = 5     # Subtile size in meters (must evenly divide tile_size; used to reduce RAM usage)

par_in["calc_trans"] = false

#tile_name = "2670000_1261000"
tile_name = "2669420_1250090"

@info dat_in["chmf"]

#par_in["save_images"] = true
#par_in["make_pngs"] = true

# === Main processing loop for each tile ===
#for job_idx in job_start:job_end

#    tile_name = string(readdlm(tiles_file)[job_idx])

#    try
        # Skip tile if it has already been processed
#        if isfile(joinpath(progress_folder, tile_name * ".txt")) || isfile(joinpath(progress_temp_folder, tile_name))
#            println("Already completed: " * tile_name)
#            continue
#        end

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

            dim = Int(tile_size / sub_tile_size)

            # === Process all subtiles within the current tile ===
            subtile_indices = [(x, y) for y in 1:dim, x in 1:dim]
            #Threads.@threads for x in 1:dim, y in 1:dim
##            Threads.@threads for i in 1:length(subtile_indices)
            for i in 1:length(subtile_indices)
                x, y = subtile_indices[i]
                idx = (limx[x] .<= pts_all[:, 1] .< limx[x+1]) .& (limy[y] .<= pts_all[:, 2] .< limy[y+1])

                if sum(idx) > 0
                    pts = pts_all[idx, :]
                    tilestr = string(limx[x]) * "_" * string(limy[y])
                    taskID = "ZH"*"_"*tilestr

                    # Check if output already exists for this subtile
                    if true #!check_output(outdir, pts[1:10, :], true, taskID)

                        # If all points are terrain (value 2), use ter2rad!
                        if sum(pts[:, 3]) == size(pts, 1) * 2
                            @info "ter2rad! " * taskID
			    func! = ter2rad!
                        else
                            @info "chm2rad! " * taskID
			    func! = chm2rad!
			end

# generate shi png
par_in["save_images"] = true
par_in["make_pngs"] = true
b1 = @benchmark			    $func!($pts, $dat_in, $par_in, $outdir, $taskID)
display(b1)

# disable shi png
par_in["save_images"] = false
par_in["make_pngs"] = false

# profile allocs
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1                            func!(pts, dat_in, par_in, outdir, taskID)
#PProf.Allocs.pprof(from_c=false)
PProf.Allocs.pprof()
@info "PROFILE ALLOCS - press enter to continue"
readline()

# profile runtime
Profile.clear()
Profile.@profile                            func!(pts, dat_in, par_in, outdir, taskID)
PProf.pprof()
@info "PROFILE RUNTIME - press enter to continue"
readline()

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

#    catch e
#        # === Error handling for failed tile ===
#        mkpath(fail_folder)
#        writedlm(joinpath(fail_folder, tile_name * ".txt"), NaN)
#        println("Failed to process tile: " * tile_name)
#        println("Caught exception: ", e)
#        println("Stack trace:")
#        for (i, frame) in enumerate(catch_backtrace())
#            println("[$i] $frame")
#        end
#    end
#end
