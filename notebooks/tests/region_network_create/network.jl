using Graphs, MetaGraphs

# function that takes dataframe containing seismic data 
# and dataframe containing cube index and cube info
# and creates a graph with Vertices being the cubes
# and Edges being the subsequent earthquakes in the cubes
# seismic database is parsed chronologically
# Seismic event x has cubeIndex XXXX and becomes node
# Seisnic event x+1 has cubeIndex YYYY and becomes node
# Between XXXX and YYYY an edge is created
# With edgeWeight you control if you take into account Edge Weights

function create_network(df, df_cubes; edgeWeight=false)

    # Initialize empty Graph 
    G = Graph(length(df_cubes.cubeIndex))
    # And metagraph (for adding properties to vertices and edges)
    MG = MetaGraph(G)

    # Use cubeIndex to add index property to all graph nodes
    for i in 1:nv(MG)
        set_prop!(MG, i, :cubeIndex, df_cubes.cubeIndex[i])
    end
    # Used to easily access the information, based on cubeIndex and not graph index
    set_indexing_prop!(MG, :cubeIndex)

    # Create Seismic Network
    i=0
    # Length of database-1
    while i < length(df.cubeIndex)-1
        # increment i first
        i+=1
        # Get the two nodes using cubeIndex indexing
        current_node = MG[df.cubeIndex[i],:cubeIndex]
        target_node = MG[df.cubeIndex[i+1],:cubeIndex]

        # Control for edgeWeight
        if edgeWeight == false
            add_edge!(MG,current_node,target_node)
        # If edgeWeight == true
        else
            # Check if there is an edge already
            if has_edge(MG, current_node, target_node)
                # If there is, increment by 1
                set_prop!(MG,current_node,target_node,
                        :weight, get_prop(MG,Edge(current_node,target_node),:weight)+1)
                # Then return to beginning (I do not want to reset the edge below)
                continue
            else
                # If there is no edge, create one with weight 1
                add_edge!(MG,current_node,target_node,:weight,1)
            end
        end
    end


    return(MG)

end


# Might delete this
# Not really necessary if you already have the df_cubes
# Keep them separately and access df_cubes when needed ?

# function that adds cube properties from df_cubes to the network
function add_properties(MG, df_cubes)
    for i in 1:nv(MG)
        set_prop!(MG, i, :xLatitude, df_cubes.xLatitude[i])
        set_prop!(MG, i, :yLongitude, df_cubes.yLongitude[i])
        set_prop!(MG, i, :zDepth, df_cubes.zDepth[i])
    end

    return MG
end