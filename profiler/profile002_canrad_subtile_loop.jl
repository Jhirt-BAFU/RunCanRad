### profile pcd2pol2cart
using Profile, PProf, BenchmarkTools
###
using CanRad, Parameters, NearestNeighbors
using Base.Threads

# === Directory paths ===
outdir = "/mnt/output/profiling"
input_path = "/mnt/input"
settings_path = ".."

# === Load user-defined settings ===
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
dat_in, par_in = C2R_Settings(input_path)

par_in["calc_trans"] = false

#par_in["terrain_highres"] = false
par_in["lowres_peri"] = 3000
#par_in["terrain_lowres"] = false
pts = [
  2669722.5 1249652.5 2.0;
  2669722.5 1249657.5 2.0;
  2669722.5 1249662.5 2.0;
  2669722.5 1249667.5 2.0;
  2669727.5 1249652.5 2.0;
  2669727.5 1249657.5 2.0;
  2669727.5 1249662.5 2.0;
  2669727.5 1249667.5 2.0;
  2669732.5 1249652.5 2.0;
  2669732.5 1249657.5 2.0;
  2669732.5 1249662.5 2.0;
  2669732.5 1249667.5 2.0;
  2669737.5 1249652.5 2.0;
  2669737.5 1249657.5 2.0;
  2669737.5 1249662.5 2.0;
  2669737.5 1249667.5 2.0;
]
print("pts: "); println(pts)

# generate shi png
par_in["save_images"] = true
par_in["make_pngs"] = true

###################################

    ################################################################################
    # Initialise

    # run compatability check then extract settings
    compatability_check!(par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    # separate the points to vectors
    pts_x, pts_y = [pts[:,ptdx] for ptdx in 1:size(pts,2)]

    canrad  = CANRAD()
    ter2rad = TER2RAD(pts_sz = size(pts_x,1))

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    outstr = splitpath(outdir)[end-1]

    dataset = CanRad.createfiles_terrain(outdir,outstr,pts,calc_trans,calc_swr)

    save_images && (images = CanRad.create_exmat_terrain(outdir,outstr,pts,canrad.mat2ev))

    save_horizon && (hlm = CanRad.create_exhlm(outdir,outstr,pts,ter2rad))

    ################################################################################
    # > Import surface data

    if terrain_highres

        limits_highres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,highres_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = CanRad.read_griddata_window(dtmf,limits_highres,true, false);

        if !isempty(dtm_x) || !isnan(sum(dtm_z))
            pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,Vector{Float64}(undef,size(pts_x)))
            rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:highres_peri)
            rows = findall(isnan,dtm_z); deleteat!(dtm_x,rows); deleteat!(dtm_y,rows); deleteat!(dtm_z,rows)
        # silly workaround to dealing with terrain data that doesn't have the same boundaries as grids.
        # this workaround loads the data including nans so that the above doesn't fail, and a horizon line
        #   can be calculated where high-res data exists. however, it will fill pts_e with nans, so this is
        #   fixed below where the low-res terrain information is used.
        # This fix is intended for the swissalti3d data which is clipped to the swiss border.
        end

    else

        dtm_x = [] # define an empty variable to run with terrain_highres only (fix in the future to something less blunt)

    end

    if terrain_lowres

        limits_lowres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,lowres_peri)
        dem_x, dem_y, dem_z, dem_cellsize = CanRad.read_griddata_window(demf,limits_lowres,true,true)
        pts_e_dem = findelev!(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,limits_lowres,50,Vector{Float64}(undef,size(pts_x)))
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:lowres_peri)

    end

    if isempty(dtm_x) || isnan(sum(pts_e))
        pts_e = copy(pts_e_dem)
    end

    ###############################################################################
    # > Image matrix  preparation

    # create an empty image matrix
    @unpack radius, mat2ev, g_rad, g_coorcrt = canrad
    g_rad[g_rad .> radius] .= NaN
    outside_img = isnan.(g_rad)
    g_coorcrt .= ((g_coorcrt .- radius) ./ radius) .* 90

    # make g_coorcrt a KDtree for easy look up
    kdtree = KDTree(g_coorcrt')

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    ###############################################################################
    # > Loop through the points

#b1 = @benchmark     for crx = 1:size(pts_x,1)
#b1 = @benchmark     Threads.@threads for crx = 1:size(pts_x,1)
b1 = @benchmark     @simd for crx = 1:size(pts_x,1)

        # get the high-res local terrain
        if !isempty(dtm_x) && terrain_highres

            pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],highres_peri);
            pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)

            if save_horizon
                dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
            end

        end

        # get the low-res regional terrain
        if terrain_lowres

            pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],lowres_peri);
            pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);

            if save_horizon
                if !isempty(dtm_x) && terrain_highres
                    dem_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                    hlm["tht"][:,crx] = Int.(round.(minimum(hcat(dtm_mintht,dem_mintht),dims=2) .*100))
                else
                    hlm["tht"][:,crx] = Int.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1]) .* 100))
                end
            end

        end

        # combine the datasets and occupy the image matrix
        fill!(mat2ev,1);
        if !isempty(dtm_x) && terrain_highres
            prepterdat!(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),13,mat2ev);
        else
            prepterdat!(pt_dem_x,pt_dem_y);
            fillmat!(canrad,kdtree,hcat(pt_dem_x,pt_dem_y),13,mat2ev);
        end


        # create the image matrix
        mat2ev[outside_img] .= 1;

        save_images && (images["SHI_terrain"][:,:,crx] = mat2ev;)

        # calculate svf and transmissivity
        svf_p, svf_h = calc_svf(canrad,mat2ev)

        dataset["svf_planar_t"][crx] = Int8(round(svf_p*100));
        dataset["svf_hemi_t"][crx]   = Int8(round(svf_h*100));

    end


###################################
# clean up

    close(dataset)
    save_images && close(images)
    save_horizon && close(hlm)

    (save_images && make_pngs) && make_SHIs(outdir,"none","none",true)

###################################

println("completed")

display(b1)
