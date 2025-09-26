### profile canrad at one camera location
using Profile, PProf, BenchmarkTools
###
using CanRad

# === Directory paths ===
outdir = "/mnt/output/profiling"
input_path = "/mnt/input"
settings_path = "/mnt/input"

# === Load user-defined settings ===
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
dat_in, par_in = C2R_Settings(input_path)

par_in["batch"] = false ; taskID = "task"
par_in["make_geotiff"] = false
pts = [2669422.5 1250092.5 2.0]

# If all points are terrain (value 2), use ter2rad!
if sum(pts[:, 3]) == size(pts, 1) * 2
    func! = ter2rad!
    @info "ter2rad!"
else
    func! = chm2rad!
    @info "chm2rad!"
end
print("pts: "); println(pts)

# generate shi png
par_in["save_images"] = true
par_in["make_pngs"] = true
# benchmark
b1 = @benchmark func!($pts, $dat_in, $par_in, $outdir, $taskID)
io = IOContext(stdout)
show(io, MIME("text/plain"), b1)

# disable shi png
par_in["save_images"] = false
par_in["make_pngs"] = false

# profile allocs
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 func!(pts, dat_in, par_in, outdir, taskID)
#PProf.Allocs.pprof(from_c=false)
PProf.Allocs.pprof()
@info "PROFILE ALLOCS - press enter to continue" ; readline()

# profile runtime
Profile.clear()
Profile.@profile func!(pts, dat_in, par_in, outdir, taskID)
PProf.pprof()
@info "PROFILE RUNTIME - press enter to continue" ; readline()

println("completed")
