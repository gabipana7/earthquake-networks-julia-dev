using DataFrames, Geodesy

# This function splits the seismic region in equally sized cubes
function region_cube_split(df; side=5, energyRelease=false)

    # Minimum and maximum of every dimension
    minLat=minimum(df.Latitude)
    maxLat=maximum(df.Latitude)
    minLon=minimum(df.Longitude)
    maxLon=maximum(df.Longitude)
    minDepth=minimum(df.Depth)
    maxDepth=maximum(df.Depth)

    # Number of cubes needed for each dimension
    # Use ceil() to round up to integer (in order to encompass whole region)
    # (some dimensions will not be perfectly divisible with the side chosen for cubes)   

    # Latitude
    x0_lla = LLA(minLat,minLon,-minDepth)
    xf_lla = LLA(maxLat,minLon,-minDepth)
    lat_dist_in_km = Geodesy.euclidean_distance(xf_lla,x0_lla) / 1000 # eucl_dist returns in meters
    # xreal, without ceiling, used when assigning earthquakes to certain cubes
    xreal = lat_dist_in_km / side
    # Actual number of cubes (integer) on Latitude distance to cover whole volume
    x = ceil(Int, lat_dist_in_km / side)

    # Longitude
    y0_lla = LLA(minLat,minLon,-minDepth)
    yf_lla = LLA(minLat,maxLon,-minDepth)
    lon_dist_in_km = Geodesy.euclidean_distance(yf_lla,y0_lla) / 1000
    # yreal, without ceiling, used when assigning earthquakes to certain cubes
    yreal = lon_dist_in_km / side
    # Actual number of cubes (integer) on Longitude distance to cover whole volume
    y = ceil(Int, lon_dist_in_km / side)

    # Depth (already in km)
    # Number of cubes on depth distance to cover whole volume
    # zreal, without ceiling, used when assigning earthquakes to certain cubes
    zreal = (maxDepth-minDepth) / side
    # Actual number of cubes (integer) on Depth distance to cover whole volume
    z = ceil(Int, (maxDepth-minDepth) / side)


    # Calculating cube Index for every value in dataframe for each dimension

    # Use xreal / yreal / zreal
    # Value will span [1,x] or [1,y] or [1,z]
    xLatitude = [floor(Int,((i-minLat)*xreal/(maxLat-minLat))+1) for i in df.Latitude]
    yLongitude = [floor(Int,((i-minLon)*yreal/(maxLon-minLon))+1) for i in df.Longitude]
    zDepth = [floor(Int,((i-minDepth)*zreal/(maxDepth-minDepth))+1) for i in df.Depth]

    # Warning, use STRING for indexing because of naming problems when making graphs!
    cubeIndex = String[]
    for i in eachindex(xLatitude)
        push!(cubeIndex, string(xLatitude[i]-1+x*(yLongitude[i]-1)+x*y*(zDepth[i]-1)))
    end

    df.cubeIndex = cubeIndex


    # Getting the center of each cube Position

    # A half cube value in every direction
    x_half_cube = (maxLat-minLat)/(2*xreal)
    y_half_cube = (maxLon-minLon)/(2*yreal)
    # For z Should always be side/2 (side=5 => 2.5 km)
    z_half_cube = (maxDepth-minDepth)/(2*zreal)

    # Formula is: minimum + half value * (2n-1) where n is a cube index in that direction
    cubeLatitude=[round(minLat + x_half_cube*(2*i-1), digits = 4) for i in xLatitude]
    cubeLongitude=[round(minLon + y_half_cube*(2*n-1), digits = 4) for n in yLongitude]
    cubeDepth=[round(minDepth + z_half_cube*(2*n-1), digits = 4) for n in zDepth]


    if energyRelease
        energyRelease=[round(10^(4.4 + 1.5*i) * 0.000001, digits=4) for i in df.Magnitude];
        df.energyRelease = energyRelease
    end


    # Separate dataframe for cubes information (reduces data redundancy)
    df_cubes = DataFrame(
    # index of the cube
    cubeIndex = cubeIndex,
    # Index on each respective dimension
    xLatitude = xLatitude,
    yLongitude = yLongitude,
    zDepth = zDepth,
    # Position of cube center in every respective dimension
    cubeLatitude = cubeLatitude,
    cubeLongitude = cubeLongitude,
    cubeDepth = cubeDepth
    )
    
    # Do not forget to unique the df (data redundancy)
    unique!(df_cubes)

    # Return original df with cubeIndex added and the cubes info df
    # df_cubes has a 1:m relationship to df 
    return df, df_cubes

end