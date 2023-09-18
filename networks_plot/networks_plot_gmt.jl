using Graphs, GraphIO
using CSV, DataFrames, FileIO

using GMT
using DelimitedFiles

include("./src/cubes.jl")
include("./src/network.jl")


function networks_plot_gmt(region,cell_size,minimum_magnitude,motif)

    # Read data
    path = "./data/"
    filepath = path * region * ".csv"
    mkpath("./maps/$region")

    df = CSV.read(filepath, DataFrame);
    df_filtered = df[df.Magnitude .> minimum_magnitude,:] 

    # Split into cubes
    df_filtered, df_filtered_cubes = region_cube_split(df_filtered,cell_size=cell_size,energyRelease=true);

    # Create network
    MG = create_network(df_filtered, df_filtered_cubes)
    connectivity = degree(MG);

    # edgelist_array = Matrix(edgelist);
    edgelist = collect(edges(MG)) |> DataFrame;


    # Get region's coordinates
    min_lon = minimum(df_filtered_cubes.cubeLongitude)
    max_lon = maximum(df_filtered_cubes.cubeLongitude)
    min_lat = minimum(df_filtered_cubes.cubeLatitude)
    max_lat = maximum(df_filtered_cubes.cubeLatitude)
    min_dep = minimum(df_filtered_cubes.cubeDepth)
    max_dep = maximum(df_filtered_cubes.cubeDepth);

    # Create the map coordinates
    map_coords = [min_lon,max_lon,min_lat,max_lat]
    map_coords_depth = [min_lon,max_lon,min_lat,max_lat,-max_dep,-min_dep]

    # Colormap for the region topography
    C_map = makecpt(cmap=:geo, range=(-8000,8000), continuous=true);
    # Relief map of the region
    relief_map = grdcut("@earth_relief_30s", region=map_coords);

    # control marker size based on degree
    marker_size = connectivity ./10

    # control marker color by connectivity
    C_markers = makecpt(cmap=:seis, range=(minimum(connectivity),maximum(connectivity)), inverse=true);

    # View of the plot
    perspective = (140,15)

    # Basemap to define the axes
    basemap(limits=map_coords_depth, proj=:merc, zsize=6, frame="SEnwZ1+b xafg yafg zaf+lDepth(km)", view=perspective)

    # Edges, plotted manually
    for i in range(1,nrow(edgelist))
        line_coords = DataFrame(lats = [df_filtered_cubes.cubeLatitude[edgelist.src[i]],df_filtered_cubes.cubeLatitude[edgelist.dst[i]]],
                        lons =[df_filtered_cubes.cubeLongitude[edgelist.src[i]],df_filtered_cubes.cubeLongitude[edgelist.dst[i]]],
                        deps= [df_filtered_cubes.cubeDepth[edgelist.src[i]],df_filtered_cubes.cubeDepth[edgelist.dst[i]]])

        plot3d!(line_coords.lons, line_coords.lats, -line_coords.deps, JZ="6c", proj=:merc, pen=(:thinner,:black), alpha=50, view=perspective)
    end

    # Nodes
    scatter3!(df_filtered_cubes.cubeLongitude, df_filtered_cubes.cubeLatitude, -df_filtered_cubes.cubeDepth,
    limits=map_coords_depth,frame="SEnwZ1+b xafg yafg zaf",proj=:merc, marker=:cube,markersize=0.1, #markersize=marker_size,
    cmap=C_markers, zcolor=connectivity, 
    alpha=50, view=perspective)

    if motif=="Triangle"
        network_target_path ="./networks/$region/cell_size_$(cell_size)km/"
        motif_filename = "motif$(motif)_$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).csv"

        motifs = readdlm(network_target_path * motif_filename, ',', Int64);

        # Motifs
        for i in range(1,size(motifs,1))
            motif_coords = DataFrame(lats = [df_filtered_cubes.cubeLatitude[motifs[i,1]], df_filtered_cubes.cubeLatitude[motifs[i,2]], df_filtered_cubes.cubeLatitude[motifs[i,3]]],
                            lons =[df_filtered_cubes.cubeLongitude[motifs[i,1]], df_filtered_cubes.cubeLongitude[motifs[i,2]], df_filtered_cubes.cubeLongitude[motifs[i,3]]],
                            deps= [df_filtered_cubes.cubeDepth[motifs[i,1]], df_filtered_cubes.cubeDepth[motifs[i,2]], df_filtered_cubes.cubeDepth[motifs[i,3]]])
            plot3d!(motif_coords.lons, motif_coords.lats, -motif_coords.deps, JZ="6c", proj=:merc, L=true, G=:red, alpha=50, view=perspective)
        end
    end

    # Colorbar
    colorbar!(limits=map_coords, pos=(paper=(20,2.0), size=(8,0.5)), shade=0.4, xaxis=(annot=2,), frame=(xlabel="Degree",),par=(MAP_LABEL_OFFSET=0.6,),view=(180,90))

    # Relief map
    grdview!(relief_map, proj=:merc, surftype=(image=1000,), 
    cmap=C_map, zsize=0.5, alpha=10 ,yshift=5.6, view=perspective,
    savefig="./networks/$region/cell_size_$(cell_size)km/$(region)_cell_size_$(cell_size)km_minmag_$(minimum_magnitude).png", 
    show=true)


end


# Romania
# pos=(paper=(16,2.0), size=(8,0.5))
# perspective = (135,20)
# networks_plot_gmt("romania", 5, 4, "")


# California
# pos=(paper=(20,2.0), size=(8,0.5))
# perspective = (140,15)
# networks_plot_gmt("california", 2.0, 4, "")


# Japan
# pos=(paper=(21,2.0), size=(8,0.5))
# perspective = (140,15)
# networks_plot_gmt("japan", 5, 5, "")


# Italy
# pos=(paper=(24,2.0), size=(8,0.5))
# perspective = (140,15)
# networks_plot_gmt("italy", 5.0, 4, "")
