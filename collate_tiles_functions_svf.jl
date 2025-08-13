using Printf


function collate2tilefile_fulltile(outdir::String, limits::Matrix{Int64}, input::String, ptsx::Vector{Float64},
	ptsy::Vector{Float64}, inputdir::String, exdir::String, par_in::Dict{String, Any},
	ptdx::Vector{Float64}, pt_spacing::Int64, fulltilesize::Int64)

	if (par_in["special_implementation"] !== "swissrad") && (par_in["special_implementation"] !== "oshd") && (par_in["special_implementation"] !== "oshd-alps")
		@warn "collate2tilefile() only written and tested for swissrad and oshd implementation of CanRad"
	end

	numptsgrid = div(fulltilesize, pt_spacing)
	numptsvec = numptsgrid^2

	# combine tile output to one file
	tiles = readdir(inputdir)

	xlim = collect(limits[1]:pt_spacing:(limits[2]-pt_spacing))
	ylim = collect(limits[3]:pt_spacing:(limits[4]-pt_spacing))

	# create the Radiation output dataset
	ds = NCDataset(joinpath(exdir, "Output_" * input * ".nc"), "c", format = :netcdf4_classic)
	defDim(ds, "locX", size(xlim, 1))
	defDim(ds, "locY", size(ylim, 1))
	defVar(ds, "Easting", Float64.(unique(ptsx)), ("locX",))
	defVar(ds, "Northing", Float64.(reverse(unique(ptsy))), ("locY",))

	#loc_time = NCDataset(joinpath(inputdir, tiles[1], "Output_" * tiles[1] * ".nc"))["datetime"][:]
	#defDim(ds, "DateTime", size(loc_time, 1))
	#time_zone = par_in["time_zone"]
	#if time_zone >= 0
	#	dt_comment = "time zone UTC+" * string(time_zone)
	#elseif time_zone < 0
	#	dt_comment = "time zone UTC-" * string(time_zone) * "; timestamp is for beginning of averaged period"
	#end
	#defVar(ds, "datetime", loc_time, ("datetime",), attrib = ["comments" => dt_comment])

	#numptstime = size(loc_time)[1]

	# create collate flags for terrain, summer and winter
	if (haskey(par_in, "calc_terrain") && par_in["calc_terrain"]) || (haskey(par_in, "SHI_terrain") && par_in["SHI_terrain"])
		terrain = true
	else
		terrain = false
	end

	if (haskey(par_in, "season") && ((par_in["season"] == "both") || (par_in["season"] == "summer"))) ||
	   (haskey(par_in, "SHI_summer") && par_in["SHI_summer"])
		summer = true
	else
		summer = false
	end

	if (haskey(par_in, "season") && ((par_in["season"] == "both") || (par_in["season"] == "winter"))) ||
	   (haskey(par_in, "SHI_winter") && par_in["SHI_winter"])
		winter = true
	else
		winter = false
	end

	for tx in tiles

		tds = NCDataset(joinpath(inputdir, tx, "Output_" * tx * ".nc"))

		# FIRST SUB-TILE
		# First sub-tile has the equal lower left coordinates as the entire 1km2 tile 
		if tx == tiles[1]
			global east  = tds["easting"][:]
			global north = tds["northing"][:]

			terrain && (global svf_p_t = tds["svf_planar_t"][:])
			terrain && (global svf_h_t = tds["svf_hemi_t"][:])
#			terrain && (global for_tau_t = tds["trans_t"][:, :])

			if sum(ptdx) != (size(ptdx, 1) * 2) # check for if the sub-tile is terrain only


				if haskey(tds, "svf_planar_w")
					winter && (global svf_p_w = tds["svf_planar_w"][:])
				else
					winter && (global svf_p_w = fill(-1, size(svf_p_t)))
				end

				if haskey(tds, "svf_planar_s")
					summer && (global svf_p_s = tds["svf_planar_s"][:])
				else
					summer && (global svf_p_s = fill(-1, size(svf_p_t)))
				end

				if haskey(tds, "svf_hemi_w")
					winter && (global svf_h_w = tds["svf_hemi_w"][:])
				else
					# winter && (global svf_h_w = fill(-1, size(svf_p_t)))
					winter && (global svf_h_w = tds["svf_hemi_t"][:])
				end

				if haskey(tds, "svf_hemi_s")
					summer && (global svf_h_s = tds["svf_hemi_s"][:])
				else
					# summer && (global svf_h_s = fill(-1, size(svf_p_t)))
					summer && (global svf_h_s = tds["svf_hemi_t"][:])
				end

				# if haskey(tds, "for_trans_w")
				# 	winter && (global for_tau_w = tds["for_trans_w"][:, :])
				# else
				# 	winter && (global for_tau_w = fill(-1, size(for_tau_t)))
				# end

				# if haskey(tds, "for_trans_s")
				# 	summer && (global for_tau_s = tds["for_trans_s"][:, :])
				# else
				# 	summer && (global for_tau_s = fill(-1, size(for_tau_t)))
				# end

			else
				winter && (global svf_p_w = fill(-1, size(svf_p_t)))
				summer && (global svf_p_s = fill(-1, size(svf_p_t)))

				winter && (global svf_h_w = fill(-1, size(svf_p_t)))
				summer && (global svf_h_s = fill(-1, size(svf_p_t)))

#				winter && (global for_tau_w = fill(-1, size(for_tau_t)))
#				summer && (global for_tau_s = fill(-1, size(for_tau_t)))
			end

		# ALL FOLLOWING SUB-TILES
		# iterates throuh this section (exept the first) 
		else
			append!(east, tds["easting"][:])
			append!(north, tds["northing"][:])

			terrain && (append!(svf_p_t, tds["svf_planar_t"][:]))
			terrain && (append!(svf_h_t, tds["svf_hemi_t"][:]))
#			terrain && (for_tau_t = hcat(for_tau_t, tds["trans_t"][:, :]))
	
			if sum(ptdx) != (size(ptdx, 1) * 2)

				dummvf = fill(-1, size(tds["svf_planar_t"][:]))

				if haskey(tds, "svf_planar_w")
					winter && (append!(svf_p_w, tds["svf_planar_w"][:]))
				else
					winter && (append!(svf_p_w, dummvf))
				end

				if haskey(tds, "svf_planar_s")
					summer && (append!(svf_p_s, tds["svf_planar_s"][:]))
				else
					summer && (append!(svf_p_s, dummvf))
				end

				if haskey(tds, "svf_hemi_w")
					winter && (append!(svf_h_w, tds["svf_hemi_w"][:]))
				else
					#winter && (append!(svf_h_w, dummvf))
					winter && (append!(svf_h_w, tds["svf_hemi_t"][:]))
				end

				if haskey(tds, "svf_hemi_s")
					summer && (append!(svf_h_s, tds["svf_hemi_s"][:]))
				else
					# summer && (append!(svf_h_s, dummvf))
					summer && (append!(svf_h_s, tds["svf_hemi_t"][:]))
				end

				# dummytau = fill(-1, size(tds["trans_t"][:, :]))
				# if haskey(tds, "for_trans_w")
				# 	winter && (for_tau_w = hcat(for_tau_w, tds["for_trans_w"][:, :]))
				# else
				# 	winter && (for_tau_w = hcat(for_tau_w, dummytau))
				# end

				# if haskey(tds, "for_trans_s")
				# 	summer && (for_tau_s = hcat(for_tau_s, tds["for_trans_s"][:, :]))
				# else
				# 	summer && (for_tau_s = hcat(for_tau_s, dummytau))
				# end

			else

				dummvf = fill(-1, size(tds["svf_planar_t"][:]))
				winter && (append!(svf_p_w, dummvf))
				summer && (append!(svf_p_s, dummvf))

				winter && (append!(svf_h_w, dummvf))
				summer && (append!(svf_h_s, dummvf))

#				dummytau = fill(-1, size(tds["trans_t"][:, :]))
#				winter && (for_tau_w = hcat(for_tau_w, dummytau))
#				summer && (for_tau_s = hcat(for_tau_s, dummytau))

			end

		end

		close(tds)

	end

	B = hcat(east, north)
	p = sortperm(view.(Ref(B), 1:size(B, 1), :))

	locs = (ptdx .> 0)


	# if winter
	# 	svf_p_all_w, svf_h_all_w, ft_newdim_w = getnewdims(svf_p_w, svf_h_w, for_tau_w, p, locs, numptsvec, numptsgrid, numptstime)
	# end

	# if summer
	# 	svf_p_all_s, svf_h_all_s, ft_newdim_s = getnewdims(svf_p_s, svf_h_s, for_tau_s, p, locs, numptsvec, numptsgrid, numptstime)
	# end

	# if terrain
	# 	svf_p_all_t, svf_h_all_t, ft_newdim_t = getnewdims(svf_p_t, svf_h_t, for_tau_t, p, locs, numptsvec, numptsgrid, numptstime)
	# end

	
	if winter
		svf_p_all_w, svf_h_all_w, ft_newdim_w = getnewdims1(svf_p_w, svf_h_w, missing, p, locs, numptsvec, numptsgrid, missing)
	end

	if summer
		svf_p_all_s, svf_h_all_s, ft_newdim_s = getnewdims1(svf_p_s, svf_h_s, missing, p, locs, numptsvec, numptsgrid, missing)
	end

	if terrain
		svf_p_all_t, svf_h_all_t, ft_newdim_t = getnewdims1(svf_p_t, svf_h_t, missing, p, locs, numptsvec, numptsgrid, missing)
	end

	if winter
		defVar(ds, "svf_planar_winter", Int8.(reverse(reshape(svf_p_all_w, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for winter (leaf-off) canopy conditions"])
		defVar(ds, "svf_hemi_winter", Int8.(reverse(reshape(svf_h_all_w, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for winter (leaf-off) canopy conditions"])
		# defVar(ds, "Forest_Transmissivity_winter", Int8.(ft_newdim_w), ("datetime", "locY", "locX"), fillvalue = Int8(-1),
		#	deflatelevel = 5, attrib = ["comments" => "calcualated for winter (leaf-off) canopy conditions"])
	end

	if summer
		defVar(ds, "svf_planar_summer", Int8.(reverse(reshape(svf_p_all_s, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for summer (leaf-on) canopy conditions"])
		defVar(ds, "svf_hemi_summer", Int8.(reverse(reshape(svf_h_all_s, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for summer (leaf-on) canopy conditions"])
		#defVar(ds, "Forest_Transmissivity_summer", Int8.(ft_newdim_s), ("datetime", "locY", "locX"), fillvalue = Int8(-1),
		#	deflatelevel = 5, attrib = ["comments" => "calcualated for summer (leaf-on) canopy conditions"])
	end

	if terrain
		defVar(ds, "svf_planar_terrain", Int8.(reverse(reshape(svf_p_all_t, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for terrain only"])
		defVar(ds, "svf_hemi_terrain", Int8.(reverse(reshape(svf_h_all_t, numptsgrid, numptsgrid), dims = 1)), ("locY", "locX"), deflatelevel = 5, fillvalue = Int8(-1),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for terrain only"])

		#defVar(ds, "Forest_Transmissivity_terrain", Int8.(ft_newdim_t), ("datetime", "locY", "locX"), fillvalue = Int8(-1),
		#	deflatelevel = 5, attrib = ["comments" => "calcualated for terrain only"])
	end

	close(ds)

end


function writeNetCDF(folder, input, pts, par_in)

	fid = joinpath(folder, input, input, "Output_" * input * ".nc")
	ds_in = NCDataset(fid)
	xcoords = ds_in["easting"][:]
	ycoords = ds_in["northing"][:]

	num_smpls = length(xcoords)
	gridsize = 1

	# define whether winter data is in netCDF
	has_winter = any(occursin.("svf_planar_w", keys(ds_in)))

	# create the Radiation output dataset
	ds = NCDataset(joinpath(folder, "Output_" * input * ".nc"), "c", format = :netcdf4_classic)
	defDim(ds, "locX", length(xcoords))
	defDim(ds, "locY", length(ycoords))
	defVar(ds, "Easting", Float64.(xcoords), ("locX",))
	defVar(ds, "Northing", Float64.(ycoords), ("locY",))

		loc_time = ds_in["datetime"][:]
		defDim(ds, "datetime", size(loc_time, 1))
		time_zone = par_in["time_zone"]
		if time_zone >= 0
			dt_comment = "time zone UTC+" * string(time_zone)
		elseif time_zone < 0
			dt_comment = "time zone UTC-" * string(time_zone) * "; timestamp is for beginning of averaged period"
		end
		defVar(ds, "datetime", loc_time, ("datetime",), attrib = ["comments" => dt_comment])

		numptstime = size(loc_time)[1]

	# create collate flags for terrain, summer and winter
	if (haskey(par_in, "calc_terrain") && par_in["calc_terrain"]) || (haskey(par_in, "SHI_terrain") && par_in["SHI_terrain"])
		terrain = true
	else
		terrain = false
	end

	if (haskey(par_in, "season") && ((par_in["season"] == "both") || (par_in["season"] == "summer"))) ||
	   (haskey(par_in, "SHI_summer") && par_in["SHI_summer"])
		summer = true
	else
		summer = false
	end

	if has_winter
		if (haskey(par_in, "season") && ((par_in["season"] == "both") || (par_in["season"] == "winter"))) ||
		   (haskey(par_in, "SHI_winter") && par_in["SHI_winter"])
			winter = true
		else
			winter = false
		end
	else
		winter = false
	end

	global east  = ds_in["easting"][:]
	global north = ds_in["northing"][:]

	terrain && (global svf_p_t = ds_in["svf_planar_t"][:])
	terrain && (global svf_h_t = ds_in["svf_hemi_t"][:])
	terrain && (global for_tau_t = ds_in["trans_t"][:, :])

	ptdx = pts[:, 3]
	if sum(ptdx) != (size(ptdx, 1) * 2) # check for if the sub-tile is terrain only
		if has_winter
			winter && (global svf_p_w = ds_in["svf_planar_w"][:])
		end
		summer && (global svf_p_s = ds_in["svf_planar_s"][:])

		if has_winter
			winter && (global svf_h_w = ds_in["svf_hemi_w"][:])
		end

		summer && (global svf_h_s = ds_in["svf_hemi_s"][:])

		if has_winter
			winter && (global for_tau_w = ds_in["for_trans_w"][:, :])
		end
		summer && (global for_tau_s = ds_in["for_trans_s"][:, :])

	else

		if has_winter
			winter && (global svf_p_w = fill(-1, size(svf_p_t)))
		end
		summer && (global svf_p_s = fill(-1, size(svf_p_t)))

		if has_winter
			winter && (global svf_h_w = fill(-1, size(svf_p_t)))
		end
		summer && (global svf_h_s = fill(-1, size(svf_p_t)))

		if has_winter
			winter && (global for_tau_w = fill(-1, size(for_tau_t)))
		end
		summer && (global for_tau_s = fill(-1, size(for_tau_t)))

	end



	if winter & has_winter
		svf_p_all_w, svf_h_all_w, ft_newdim_w = getnewdims_smplspts(svf_p_w, svf_h_w, for_tau_w, gridsize, num_smpls, numptstime)
	end

	if summer
		svf_p_all_s, svf_h_all_s, ft_newdim_s = getnewdims_smplspts(svf_p_s, svf_h_s, for_tau_s, gridsize, num_smpls, numptstime)
	end

	if terrain
		svf_p_all_t, svf_h_all_t, ft_newdim_t = getnewdims_smplspts(svf_p_t, svf_h_t, for_tau_t, gridsize, num_smpls, numptstime)
	end

	if winter & has_winter
		defVar(ds, "svf_planar_winter", Int8.(reshape(svf_p_all_w, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for winter (leaf-off) canopy conditions"])

		defVar(ds, "svf_hemi_winter", Int8.(reshape(svf_h_all_w, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for winter (leaf-off) canopy conditions"])

		defVar(ds, "Forest_Transmissivity_winter", Int8.(ft_newdim_w), ("datetime", "locX"), attrib = ["comments" => "calcualated for winter (leaf-off) canopy conditions"])
	end

	if summer
		defVar(ds, "svf_planar_summer", Int8.(reshape(svf_p_all_s, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for summer (leaf-on) canopy conditions"])
		defVar(ds, "svf_hemi_summer", Int8.(reshape(svf_h_all_s, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for summer (leaf-on) canopy conditions"])
		defVar(ds, "Forest_Transmissivity_summer", Int8.(ft_newdim_s), ("datetime", "locX"), attrib = ["comments" => "calcualated for summer (leaf-on) canopy conditions"])
	end

	if terrain
		defVar(ds, "svf_planar_terrain", Int8.(reshape(svf_p_all_t, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of a horizontal flat uplooking surface;
									zenith rings weighted by surface area projected onto a horizontal flat surface;
									calcualated for terrain only"])
		defVar(ds, "svf_hemi_terrain", Int8.(reshape(svf_h_all_t, num_smpls, gridsize)), ("locX",),
			attrib = ["comments" => "values are represented as percentage;
									perspective of hemipherically shaped surface or plant;
									zenith rings weighted by surface area on the hemisphere;
									calcualated for terrain only"])

		defVar(ds, "Forest_Transmissivity_terrain", Int8.(ft_newdim_t), ("datetime", "locX"), attrib = ["comments" => "calcualated for terrain only"])
	end

	close(ds)
	close(ds_in)

end

function getnewdims_smplspts(svf_p, svf_h, for_tau, gridsize, num_smpls, numptstime)

	svf_p[ismissing.(svf_p)] .= -1
	svf_p_all = fill(-1, num_smpls, 1)
	svf_p_all[:] = svf_p[:]
	# svf_p_all[locs] .= (svf_p[p])

	svf_h[ismissing.(svf_h)] .= -1
	svf_h_all = fill(-1, num_smpls, 1)
	svf_h_all = svf_h
	# svf_h_all[locs] .= (svf_h[p])

	for_tau[ismissing.(for_tau)] .= -1
	ft_all = fill(-1, numptstime, num_smpls)
	# ft_all[:,locs] .= for_tau[:,p]
	ft_all = for_tau
	ft_newdim = fill(-1, numptstime, num_smpls, gridsize)
	for x in 1:1:size(ft_all, 1)
		ft_newdim[x, :, :] = reshape(ft_all[x, :], num_smpls, gridsize)
	end

	return svf_p_all, svf_h_all, ft_newdim

end

function getnewdims1(svf_p,svf_h,for_tau,p,locs,numptsvec,numptsgrid,numptstime)

    svf_p[ismissing.(svf_p)] .= -1
    svf_p_all = fill(-1,numptsvec,1)
    svf_p_all[locs] .= (svf_p[p])

    svf_h[ismissing.(svf_h)] .= -1
    svf_h_all = fill(-1,numptsvec,1)
    svf_h_all[locs] .= (svf_h[p])

    # for_tau[ismissing.(for_tau)] .= -1
    # ft_all = fill(-1,numptstime,numptsvec)
    # ft_all[:,locs] .= for_tau[:,p]
    # ft_newdim = fill(-1,numptstime,numptsgrid,numptsgrid)
    # for x in 1:1:size(ft_all,1)
    #     ft_newdim[x,:,:] = reverse(reshape(ft_all[x,:],numptsgrid,numptsgrid),dims=1)
    # end
    
    return svf_p_all, svf_h_all, missing #ft_newdim

end