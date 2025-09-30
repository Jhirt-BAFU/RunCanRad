### profile pcd2pol2cart
using Profile, PProf, BenchmarkTools
###
using CanRad, Parameters

# === Directory paths ===
outdir = "/mnt/output/profiling"
input_path = "/mnt/input"
settings_path = ".."

# === Load user-defined settings ===
include(joinpath(settings_path, "C2R_Settings_CHCW_local.jl"))
dat_in, par_in = C2R_Settings(input_path)

par_in["calc_trans"] = false

#par_in["terrain_highres"] = false
#par_in["lowres_peri"] = 3000
par_in["terrain_lowres"] = false
pts = [
  2669722.5 1249652.5;
]
print("pts: "); println(pts)

# generate shi png
#par_in["save_horizon"] = true
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

#    dataset = CanRad.createfiles_terrain(outdir,outstr,pts,calc_trans,calc_swr)

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
    kdtree = scipyspat.cKDTree(g_coorcrt)

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    ###############################################################################
    # process the first point only - instead of looping through all the pts points
    #    @simd for crx = 1:size(pts_x,1)
    crx = 1 # unloop


        # get the high-res local terrain
        if !isempty(dtm_x) && terrain_highres

###################################
## benchmark highres
@info "highres: radius "*string(highres_peri)*"m"
@info string(length(dtm_x))*" points within bbox"
@info minimum(dtm_x), maximum(dtm_x), minimum(dtm_y), maximum(dtm_y), minimum(dtm_z), maximum(dtm_z)
            pt_dtm_x0, pt_dtm_y0, pt_dtm_z0 = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],highres_peri);
@info string(length(pt_dtm_x0))*" points within 3d radius (<pi/4)"
@info minimum(pt_dtm_x0), maximum(pt_dtm_x0), minimum(pt_dtm_y0), maximum(pt_dtm_y0), minimum(pt_dtm_z0), maximum(pt_dtm_z0)
            pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,copy(pt_dtm_x0), copy(pt_dtm_y0), copy(pt_dtm_z0),pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)
b1 = @benchmark            pt_dtm_x, pt_dtm_y = pcd2pol2cart!($ter2rad,copy($pt_dtm_x0), copy($pt_dtm_y0), copy($pt_dtm_z0),$pts_x[$crx],$pts_y[$crx],$pts_e[$crx],"terrain",$rbins_dtm,$image_height)
@info string(length(pt_dtm_x))*" resulting points in azimuthal equidistant projection"
@info minimum(pt_dtm_x), maximum(pt_dtm_x), minimum(pt_dtm_y), maximum(pt_dtm_y)
dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
@info string(length(dtm_mintht))*" results on horizon line"
@info minimum(dtm_mintht), maximum(dtm_mintht)
display(b1); println()
###################################

            if save_horizon
                dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
            end

        end

        # get the low-res regional terrain
        if terrain_lowres

###################################
## benchmark lowres
@info "lowres: radius "*string(lowres_peri)*"m"
@info string(length(dem_x))*" points within bbox"
@info minimum(dem_x), maximum(dem_x), minimum(dem_y), maximum(dem_y), minimum(dem_z), maximum(dem_z)
            pt_dem_x0, pt_dem_y0, pt_dem_z0 = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],lowres_peri)
@info string(length(pt_dem_x0))*" points within 3d radius (<pi/4)"
@info minimum(pt_dem_x0), maximum(pt_dem_x0), minimum(pt_dem_y0), maximum(pt_dem_y0), minimum(pt_dem_z0), maximum(pt_dem_z0)
            pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,copy(pt_dem_x0), copy(pt_dem_y0), copy(pt_dem_z0),pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height)
b2 = @benchmark            pt_dem_x, pt_dem_y = pcd2pol2cart!($ter2rad,copy($pt_dem_x0), copy($pt_dem_y0), copy($pt_dem_z0),$pts_x[$crx],$pts_y[$crx],$pts_e_dem[$crx],"terrain",$rbins_dem,$image_height)
@info string(length(pt_dem_x))*" resulting points in azimuthal equidistant projection"
@info minimum(pt_dem_x), maximum(pt_dem_x), minimum(pt_dem_y), maximum(pt_dem_y)
dem_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
@info string(length(dem_mintht))*" results on horizon line"
@info minimum(dem_mintht), maximum(dem_mintht)
display(b2); println()

###################################

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
        if !isempty(dtm_x) && terrain_highres && terrain_lowres
            prepterdat!(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),13,mat2ev);
        elseif terrain_lowres
            prepterdat!(pt_dem_x,pt_dem_y);
            fillmat!(canrad,kdtree,hcat(pt_dem_x,pt_dem_y),13,mat2ev);
	else
	    prepterdat!(pt_dtm_x,pt_dtm_y);
	    fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),13,mat2ev);
        end


        # create the image matrix
        mat2ev[outside_img] .= 1;

        save_images && (images["SHI_terrain"][:,:,crx] = mat2ev;)

        # calculate svf and transmissivity
        svf_p, svf_h = calc_svf(canrad,mat2ev)
@info svf_p, svf_h

#        dataset["svf_planar_t"][crx] = Int8(round(svf_p*100));
#        dataset["svf_hemi_t"][crx]   = Int8(round(svf_h*100));

    #    end


###################################
## clean up
#    close(dataset)
    save_images && close(images)
    save_horizon && close(hlm)

    (save_images && make_pngs) && make_SHIs(outdir,"none","none",true)


###################################
## profiling

# profile allocs
Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=1 pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,copy(pt_dtm_x0), copy(pt_dtm_y0), copy(pt_dtm_z0),pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)
#PProf.Allocs.pprof(from_c=false)
PProf.Allocs.pprof()
@info "PROFILE ALLOCS - press enter to continue" ; readline()

# profile runtime
Profile.clear()
Profile.@profile pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,copy(pt_dtm_x0), copy(pt_dtm_y0), copy(pt_dtm_z0),pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)
PProf.pprof()
@info "PROFILE RUNTIME - press enter to continue" ; readline()
###################################

println("completed")
