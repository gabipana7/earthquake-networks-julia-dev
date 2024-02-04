using CSV, DataFrames, GMT
using Graphs, GraphIO, FileIO
using Geodesy

include("./cubes.jl")
include("./network.jl")

function network_plot_gmt(region::String, cell_size::Float64, minimum_magnitude::Int64, marker_size::Float64, ratio_for_z::Int64, legend_info::Vector{Any}, perspective::Tuple{Int64,Int64}, quality::Int64, colorbar_info::Vector{Any})
    """
    Description:
    Takes in the region, cell_size and minimum_magnitude and plots with GMT
        a 3D visualization of the earthquake network with the above specifications, 
        with a semi-3D relief map of the region.

    Parameters:
    - `region`: a string (capital letter) the seismic region under scrutiny.
    - `cell_size`: discretization length of the cubes being split into.
    - `minimum_magnitude`: minimum magnitude to filter the database
    - `persepctive`: tuple (x,y) to set the view of the 3D plot
    - `quality`: the quality of the relief map, higher is better (>1000), and computationally demanding

    Returns:
        Plots the GMT visualization.
    """
    
    ##################################### READ DATA #####################################
    path = "./data/"
    filepath = path * region * ".csv"
    mkpath("./gmt/$region")
    df = CSV.read(filepath, DataFrame);
    df_filtered = df[df.Magnitude .> minimum_magnitude,:] 
    #########################################################################################

    ##################################### MAP OF REGION #####################################
    #########################################################################################
    # Split into cubes
    df_filtered, df_filtered_cubes = region_cube_split(df_filtered,cell_size=cell_size,energyRelease=true);
    # Region coords based on cubes
    min_lon = minimum(df_filtered_cubes.cubeLongitude) -0.1
    max_lon = maximum(df_filtered_cubes.cubeLongitude) +0.1
    min_lat = minimum(df_filtered_cubes.cubeLatitude) -0.1
    max_lat = maximum(df_filtered_cubes.cubeLatitude) +0.1
    min_dep = minimum(df_filtered_cubes.cubeDepth);
    max_dep = maximum(df_filtered_cubes.cubeDepth);
    # Create the map coordinates
    map_coords = [min_lon,max_lon,min_lat,max_lat]
    map_coords_depth = [min_lon,max_lon,min_lat,max_lat,-max_dep,-min_dep]
    # Colormap for the region topography
    C_map = makecpt(cmap=:geo, range=(-8000,8000), continuous=true);
    # Relief map of the region
    relief_map = grdcut("@earth_relief_30s", region=map_coords);

    # z_ratio (ratio between longitude and depth, exaggerated by a factor of 2 for better vis)
    x0_lla = LLA(min_lat,min_lon,-min_dep)
    xf_lla = LLA(min_lat,max_lon,-min_dep)
    lon_dist_in_km = Geodesy.euclidean_distance(xf_lla,x0_lla) / 1000
    z_ratio = ratio_for_z * max_dep / lon_dist_in_km

    #########################################################################################

    ####################################### GMT PLOT ########################################
    #########################################################################################
    #Basemap to define the axes, annot and labels
    basemap(limits=map_coords_depth, proj=:merc, zsize=z_ratio, frame="SEnwZ1+b xafg yafg zaf+lDepth(km)", 
            par=(FONT_LABEL="20p,Palatino-Roman",MAP_LABEL_OFFSET=0.8, FONT_ANNOT="18p,Palatino-Roman,black"), view=perspective)
    
    
    # Edges, plotted manually
    # Create network
    MG = create_network(df_filtered, df_filtered_cubes)
    # edgelist_array = Matrix(edgelist); for edges 
    edgelist = collect(edges(MG)) |> DataFrame;
    for i in range(1,nrow(edgelist))
        line_coords = DataFrame(lats = [df_filtered_cubes.cubeLatitude[edgelist.src[i]],df_filtered_cubes.cubeLatitude[edgelist.dst[i]]],
                        lons =[df_filtered_cubes.cubeLongitude[edgelist.src[i]],df_filtered_cubes.cubeLongitude[edgelist.dst[i]]],
                        deps= [df_filtered_cubes.cubeDepth[edgelist.src[i]],df_filtered_cubes.cubeDepth[edgelist.dst[i]]])
    
        plot3d!(line_coords.lons, line_coords.lats, -line_coords.deps, JZ=string(z_ratio) * "c", proj=:merc, pen=(0.3,:gray), alpha=45, view=perspective)
    end
    
    # Nodes
    # connectivity for nodes degree
    connectivity = degree(MG);
    # C_markers = makecpt(cmap=:berlin, range=(minimum(connectivity),maximum(connectivity)), inverse=true);
    # C_markers = makecpt(cmap=:berlin, range=(minimum(connectivity),maximum(connectivity)));
    C_markers = makecpt(cmap=:split, range=(minimum(connectivity),maximum(connectivity)));

    scatter3!(df_filtered_cubes.cubeLongitude, df_filtered_cubes.cubeLatitude, -df_filtered_cubes.cubeDepth,
                limits=map_coords_depth,frame="SEnwZ1+b",proj=:merc, marker=:cube, markersize=marker_size,
                cmap=C_markers, zcolor=connectivity, alpha=60, view=perspective)
    
    # Colorbar
    colorbar!(limits=map_coords, pos=(paper=colorbar_info[1], size=colorbar_info[2]), shade=0.4, xaxis=(annot=colorbar_info[3]), frame=(xlabel="Degree",),
                par=(FONT_LABEL="22p,Palatino-Roman,black",MAP_LABEL_OFFSET=0.6,FONT_ANNOT="18p,Palatino-Roman,black"),view=(180,90))
    
    # Relief map
    grdview!(relief_map, proj=:merc, surftype=(image=quality,), 
                cmap=C_map, zsize=1.0, alpha=45 ,yshift=z_ratio-0.3, view=perspective,
                # savefig="./gmt/$region/$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).png",
                )



    legend!((label=(txt="L=$(cell_size) , M"*subscript(min)*"=$(minimum_magnitude) ", justify=:L, font=(14, "Palatino-Roman")),),
                limits=map_coords, pos=(paper=legend_info[1], width=legend_info[2], justify=:BL, spacing=2.4),
                clearance=(0.1,0.1), box=(pen=0.1, fill=:azure1),
                figsize=10, proj=:Mercator, view = (180,90),
                savefig="./gmt/$region/$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).png",)
    #########################################################################################
end








#############################################################################################################################################
function network_plot_gmt_Vrancea(region::String, cell_size::Float64, minimum_magnitude::Int64, marker_size::Float64, ratio_for_z::Int64, legend_info::Vector{Any}, perspective::Tuple{Int64,Int64}, quality::Int64, colorbar_info::Vector{Any})
    """
    Description:
    Takes in the region, cell_size and minimum_magnitude and plots with GMT
        a 3D visualization of the earthquake network with the above specifications, 
        with a semi-3D relief map of the region.

    Parameters:
    - `region`: a string (capital letter) the seismic region under scrutiny.
    - `cell_size`: discretization length of the cubes being split into.
    - `minimum_magnitude`: minimum magnitude to filter the database
    - `persepctive`: tuple (x,y) to set the view of the 3D plot
    - `quality`: the quality of the relief map, higher is better (>1000), and computationally demanding

    Returns:
        Plots the GMT visualization.
    """
    
    ##################################### READ DATA #####################################
    path = "./data/"
    filepath = path * "Romania" * ".csv"
    mkpath("./gmt/$region")
    df = CSV.read(filepath, DataFrame);
    
    # Read Vrancea
    # Read data
    path = "./data/"
    filepath_Vrancea = path * region * ".csv"

    df_Vrancea = CSV.read(filepath_Vrancea, DataFrame);
    minlon_Vrancea = minimum(df_Vrancea.Longitude)
    maxlon_Vrancea = maximum(df_Vrancea.Longitude)
    minlat_Vrancea = minimum(df_Vrancea.Latitude)
    maxlat_Vrancea = maximum(df_Vrancea.Latitude)

    df = df[(df.Longitude .> minlon_Vrancea) .& (df.Longitude .< maxlon_Vrancea) .& 
        (df.Latitude .> minlat_Vrancea) .& (df.Latitude .< maxlat_Vrancea), :]

    df_filtered = df[df.Magnitude .> minimum_magnitude,:] 
    #########################################################################################

    ##################################### MAP OF REGION #####################################
    #########################################################################################
    # Split into cubes
    df_filtered, df_filtered_cubes = region_cube_split(df_filtered,cell_size=cell_size,energyRelease=true);
    # Region coords based on cubes
    min_lon = minimum(df_filtered_cubes.cubeLongitude) -0.2
    max_lon = maximum(df_filtered_cubes.cubeLongitude) +0.2
    min_lat = minimum(df_filtered_cubes.cubeLatitude) -0.2
    max_lat = maximum(df_filtered_cubes.cubeLatitude) +0.2
    min_dep = minimum(df_filtered_cubes.cubeDepth);
    max_dep = maximum(df_filtered_cubes.cubeDepth);
    # Create the map coordinates
    map_coords = [min_lon,max_lon,min_lat,max_lat]
    map_coords_depth = [min_lon,max_lon,min_lat,max_lat,-max_dep,-min_dep]
    # Colormap for the region topography
    C_map = makecpt(cmap=:geo, range=(-8000,8000), continuous=true);
    # Relief map of the region
    relief_map = grdcut("@earth_relief_30s", region=map_coords);

    # z_ratio (ratio between longitude and depth, exaggerated by a factor of 2 for better vis)
    x0_lla = LLA(min_lat,min_lon,-min_dep)
    xf_lla = LLA(min_lat,max_lon,-min_dep)
    lon_dist_in_km = Geodesy.euclidean_distance(xf_lla,x0_lla) / 1000
    z_ratio = ratio_for_z * max_dep / lon_dist_in_km

    #########################################################################################

    ####################################### GMT PLOT ########################################
    #########################################################################################
    #Basemap to define the axes, annot and labels
    basemap(limits=map_coords_depth, proj=:merc, zsize=z_ratio, frame="SEnwZ1+b xafg yafg zaf+lDepth(km)", 
            par=(FONT_LABEL="20p,Palatino-Roman",MAP_LABEL_OFFSET=0.8, FONT_ANNOT="18p,Palatino-Roman,black"), view=perspective)
    
    
    # Edges, plotted manually
    # Create network
    MG = create_network(df_filtered, df_filtered_cubes)
    # edgelist_array = Matrix(edgelist); for edges 
    edgelist = collect(edges(MG)) |> DataFrame;
    for i in range(1,nrow(edgelist))
        line_coords = DataFrame(lats = [df_filtered_cubes.cubeLatitude[edgelist.src[i]],df_filtered_cubes.cubeLatitude[edgelist.dst[i]]],
                        lons =[df_filtered_cubes.cubeLongitude[edgelist.src[i]],df_filtered_cubes.cubeLongitude[edgelist.dst[i]]],
                        deps= [df_filtered_cubes.cubeDepth[edgelist.src[i]],df_filtered_cubes.cubeDepth[edgelist.dst[i]]])
    
        plot3d!(line_coords.lons, line_coords.lats, -line_coords.deps, JZ=string(z_ratio) * "c", proj=:merc, pen=(0.3,:black), alpha=65, view=perspective)
    end
    
    # Nodes
    # connectivity for nodes degree
    connectivity = degree(MG);
    # C_markers = makecpt(cmap=:berlin, range=(minimum(connectivity),maximum(connectivity)), inverse=true);
    # C_markers = makecpt(cmap=:berlin, range=(minimum(connectivity),maximum(connectivity)));
    C_markers = makecpt(cmap=:split, range=(minimum(connectivity),maximum(connectivity)));

    scatter3!(df_filtered_cubes.cubeLongitude, df_filtered_cubes.cubeLatitude, -df_filtered_cubes.cubeDepth,
                limits=map_coords_depth,frame="SEnwZ1+b",proj=:merc, marker=:cube, markersize=marker_size,
                cmap=C_markers, zcolor=connectivity, alpha=60, view=perspective)
    
    # Colorbar
    colorbar!(limits=map_coords, pos=(paper=colorbar_info[1], size=colorbar_info[2]), shade=0.4, xaxis=(annot=colorbar_info[3]), frame=(xlabel="Degree",),
                par=(FONT_LABEL="22p,Palatino-Roman,black",MAP_LABEL_OFFSET=0.6,FONT_ANNOT="18p,Palatino-Roman,black"),view=(180,90))
    
    # Relief map
    grdview!(relief_map, proj=:merc, surftype=(image=quality,), 
                cmap=C_map, zsize=1.0, alpha=25 ,yshift=z_ratio-0.6, view=perspective,
                # savefig="./gmt/$region/$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).png",
                )


    legend!((label=(txt="L=$(cell_size) , M"*subscript(min)*"=$(minimum_magnitude) ", justify=:L, font=(20, "Palatino-Roman")),),
                limits=map_coords, pos=(paper=legend_info[1], width=legend_info[2], justify=:BL, spacing=2.4),
                clearance=(0.1,0.1), box=(pen=0.1, fill=:azure1),
                figsize=10, proj=:Mercator, view = (180,90),
                savefig="./gmt/$region/$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).png",)
    #########################################################################################
end
#############################################################################################################################################