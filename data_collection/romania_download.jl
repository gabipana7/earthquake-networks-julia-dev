using CSV, DataFrames, Dates

url = "http://www.infp.ro/data/romplus.txt"
filename = "./romania.txt"
download(url,filename)


# Declare types of the needed columns
datetime = Vector{String}()
latitude, longitude, depth =  Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
magnitude =  Vector{Float64}()

# open the file
open(filename) do io
    # skip the first line
    line = readline(io)
    while true
        # start reading the lines
        line = readline(io)
        # detect the end of the line
        line == "" && break
        # push to the vectors, the characters based on the position in the text file
        # push!(date,  strip(line[1:10]))
        push!(datetime,  strip(line[1:23]))
        push!(latitude, parse(Float64,strip(line[38:45])))
        push!(longitude, parse(Float64, strip(line[48:56])))
        push!(depth, parse(Float64, strip(line[76:80])))
        push!(magnitude, parse(Float64, strip(line[108:110])))
    end
end

# Initialize dataframe
df = DataFrame(Datetime=datetime, 
                Latitude=latitude, Longitude=longitude, Depth=depth,
                Magnitude=magnitude)

# Turn to datetime 
dateformat = dateformat"yyyy/mm/dd  HH:MM:SS.s"
df.Datetime = DateTime.(df.Datetime, dateformat)

# Filter magnitudes
romania = df[df.Magnitude .> 0.0,:]

# Write to CSV
mkpath("./data/")
CSV.write("./data/romania.csv", romania)

# Delete original file
rm("romania.txt")