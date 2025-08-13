function C2R_Settings(datfolder)

    par_in = Dict(
    
    "special_implementation" => "swissrad", # "none" or "swissrad" 

    # terrain settings
    "terrain_highres"     => true,  # include high high resolution local terrain (1-5 m)
    "terrain_lowres"      => true,  # include low resolution regional terrain (> 25 m)
    "lowres_peri"         => 25000, # radius for low-res terrain (calculating regional horizon line)
    "terrainmask_precalc" => false, # if the terrain mask has been pre-calculated by CanRad
    "calc_terrain"        => true,  # calculate variables for terrain-only as well as forest
    "buildings"           => false, # include buildings in terrain calculation - requires blmf file
    "hlm_precalc"         => false,  # use pre-calculated horizon line matrix if true - make sure "hlmf" is specified below, terrain_lowres = false and terrain_tile = true

    # forest settings
    "forest_type"  => "mixed", # evergreen, deciduous, mixed
    "tree_species" => "both",  # needleleaf, broadleaf, both
    "season"       => "both",  # summer, winter, both
    "surf_peri"    => 100,     # radius within which forest is included in SHI calculation
    "trunks"       => false,
    
    # calculation settings
    "calc_trans" => true,
    "calc_swr"   => 0, # 0 = off; 1 = potential swr (atm_trans = 1)

    "t1"    => "01.01.2020 00:00:00", # "dd.mm.yyyy HH:MM:SS"
    "t2"    => "31.12.2020 23:00:00",
    "tstep" => 60,

    "time_zone"   => 1,
    "coor_system" => "CH1903+",

    # camera settings
    "image_height" => 0.5, # enter float value
    "tilt"         => false,

    # run settings
    "batch"        => true, # running in parallel or single process
    "save_images"  => true,
    "make_pngs"    => false,
    "save_horizon" => false, # save the calculated terrain horizon line  
    "progress"     => false # save progress of each calculation step
    )

    dat_in = Dict(
    # Vegetation Hight Modle with 1m resolution
    "chmf" => joinpath(datfolder, "20250217_VHM_ALS_CH_SWISS_1m_2056.tif"),
    # Hight Model: Digital terain model from swisstopo with 5 m resolution
     "dtmf" => joinpath(datfolder,"swissALTI3D_5M_CHLV95_LN02_2020.tif"),
    # Hight Model ??
    "demf" => joinpath(datfolder,"BAFU_COP_DEM_merged_50m_2056.tif"),
    # Mixed Ratio Data
     "mrdf" => joinpath(datfolder,"WMG2020_RF_Geobasisdatensatz_nomask_LV95.tif"),
    # Forest/non-forest classified data
     "fcdf" => joinpath(datfolder,"ForestClass.tif"),
    # Forest Type Data
    "ftdf" => joinpath(datfolder,"Copernicus_ForestType_UInt32.tif")
    )

    return dat_in, par_in


end
