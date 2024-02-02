using Geodesy

#-------------------------MOTIFS FILE PROCESSING--------------------------------#

#-------------------------TOTAL / MEAN  ENERGY--------------------------------#

# Works for any legth motifs
function total_mean_energy(motifs,df,df_cubes)
    #=
    motifs = Matrix{Int64}
    df = DataFrame
    df_cubes = DataFrame
    =#

    # STEP 1
    # Calculate total energy and mean energy only for cubes that are part of motifs
    # This is done to reduce calculations. I do not need to calculate these for each motif
    # (because, nodes are repeated in motifs.)
    # (I calculate separately, then I just add them up for each motif)

    # Get all the cubes, uniquely
    vec_motifs = unique!(vec(motifs));

    # Will apend to a dictionary: cubeIndex => [totalEnergy, meanEnergy]
    cube_energy = Dict()
    for value in vec_motifs
        # Get the cube from df_cubes
        respective_cube = df_cubes.cubeIndex[value]
        # Use the cube to get all the quakes in that cube from df
        quakesInNode = df[df.cubeIndex .== respective_cube,:];

        # Calculate total energy and mean energy
        totalEnergy = sum(quakesInNode.energyRelease)
        meanEnergy = totalEnergy / length(quakesInNode.energyRelease)

        # Push to dictionary
        cube_energy[value] = [totalEnergy, meanEnergy]
    end


    # STEP 2
    # Caclulate total energy and mean energy per each motif
    motif_energy = Dict()
    # Parse along the motifs
    for i=1:length(motifs[:,1])
        # initialize with zero
        total_energy_in_motif=0
        mean_energy_in_motif=0
        # parse each motif elements (3 for Triangles, 4 for Tetrahedrons)
        for j=1:length(motifs[1,:])
            # add to total energy and mean energy
            total_energy_in_motif += cube_energy[motifs[i,j]][1] 
            mean_energy_in_motif += cube_energy[motifs[i,j]][2]  
        end
        # Add to dictionary
        motif_energy[i] = [total_energy_in_motif,mean_energy_in_motif]
    end


    return motif_energy

end




#-------------------------TRIANGLES--------------------------------#

# Calculate the areas per motif
function area_triangles(motifs,df_cubes)
    areas = Dict()
    for i=1:length(motifs[:,1])
        x = Vector{Any}(undef,3)
        for j=1:3

            lat = df_cubes.cubeLatitude[motifs[i,j]]
            lon = df_cubes.cubeLongitude[motifs[i,j]]
            dep = df_cubes.cubeDepth[motifs[i,j]]

            x[j] = LLA(lat,lon,-1000*dep)

        end

        a = Geodesy.euclidean_distance(x[1],x[2]) / 1000
        b = Geodesy.euclidean_distance(x[2],x[3]) / 1000
        c = Geodesy.euclidean_distance(x[1],x[3]) / 1000;

        # calculate semiperimeter
        sp = (a+b+c)/2
                
        # Use Heron's formula Area = sqrt(semiperimeter(sp-a)(sp-b)(sp-c))
        # We have the areas of the triangle, append them to lists ! 

        A = sqrt(abs(sp*(sp-a)*(sp-b)*(sp-c)))

        areas[i] = A
    end

    return areas
end





#-------------------------TETRAHEDRONS--------------------------------#

# Calculate the volumes per motif
function volume_tetrahedrons(motifs,df_cubes)
    volumes = Dict()
    for i=1:length(motifs[:,1])
        # Establish the four points of the tetrahedron LLA(lat,lon,dep)
        x = Vector{Any}(undef,4)
        for j=1:4

            lat = df_cubes.cubeLatitude[motifs[i,j]]
            lon = df_cubes.cubeLongitude[motifs[i,j]]
            dep = df_cubes.cubeDepth[motifs[i,j]]

            x[j] = LLA(lat,lon,-dep)

        end

        # Calculate preliminary elements using tetrahedron sides (in km)
        W = Geodesy.euclidean_distance(x[1],x[2]) / 1000
        V = Geodesy.euclidean_distance(x[2],x[3]) / 1000
        U = Geodesy.euclidean_distance(x[1],x[3]) / 1000
        u = Geodesy.euclidean_distance(x[2],x[4]) / 1000
        v = Geodesy.euclidean_distance(x[1],x[4]) / 1000
        w = Geodesy.euclidean_distance(x[3],x[4]) / 1000

        # calculate elements that go into elements that go into Heron formula
        A = (w-U+v)*(U+v+w)
        B = (u-V+w)*(V+w+u)
        C = (v-W+u)*(W+u+v)
        
        a = (U-v+w)*(v-w+U)
        b = (V-w+u)*(w-u+V)
        c = (W-u+v)*(u-v+W)

        # elements that go into Heron formula
        p = sqrt(abs(a*B*C))
        q = sqrt(abs(b*C*A))
        r = sqrt(abs(c*A*B))
        s = sqrt(abs(a*b*c))
                
        # Use Heron's formula Area = sqrt(semiperimeter(sp-a)(sp-b)(sp-c))
        # We have the areas of the triangle, append them to lists ! 

        V = sqrt(abs((-p+q+r+s)*(p-q+r+s)*(p+q-r+s)*(p+q+r-s)))/(192*u*v*w)

        volumes[i] = V
    end

    return volumes
end