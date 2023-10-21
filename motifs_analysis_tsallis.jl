using CSV, DataFrames
using FileIO, Dates
using PyCall
using DelimitedFiles
using CairoMakie

include("./src/cubes.jl")
include("./src/network.jl")
include("./src/motifs_analysis.jl")

@pyimport powerlaw as powlaw
so = pyimport("scipy.optimize")


# Triangles Analysis
function analize_motifs_triangle_tsallis(region, weighted_by)
    # Read data
    path = "./data/"
    filepath = path * region * ".csv"
    df = CSV.read(filepath, DataFrame);

    motif="Triangle"
    if weighted_by == "totalenergy"
        weight_key = 1
    else 
        weight_key = 2
    end

    # Make path for results
    mkpath("./motifs/tsallis/$weighted_by/$region")

    # Based on parameter dependency, extract which cell_size lengths are the best:
    if region == "Romania"
        cell_sizes = [3.0, 3.5, 4.0, 4.5, 5.0];
        minimum_magnitudes = [0,1,2,3];
    elseif region == "Italy"
        cell_sizes = [4.0, 4.5, 5.0, 5.5, 6.0];
        minimum_magnitudes = [0,1, 2,3];
    elseif region == "California"
        cell_sizes = [1.0, 1.5, 2.0];
        minimum_magnitudes = [0,1,2,3];
    elseif region == "Japan"
        cell_sizes = [2.5, 3.0, 3.5, 4.0, 5.0];
        minimum_magnitudes = [2,3,4];
    end;

    for cell_size in cell_sizes
        for minimum_magnitude in minimum_magnitudes
            # Filter by magnitude
            df_filtered = df[df.Magnitude .> minimum_magnitude,:] 

            # Split into cubes
            df_filtered, df_filtered_cubes = region_cube_split(df_filtered,cell_size=cell_size,energyRelease=true);

            # Get the motif
            network_target_path = "./networks/$(region)/cell_size_$(string(cell_size))km/"
            motif_filename = "motif$(motif)_$(region)_cell_size_$(string(cell_size))km_minmag_$(string(minimum_magnitude)).csv"

            # motifs = CSV.read(network_target_path * motif_filename, DataFrame);
            motifs = readdlm(network_target_path * motif_filename, ',', Int64);
            
            # Energy and areas calculator
            motif_energy = total_mean_energy(motifs, df_filtered, df_filtered_cubes);
            areas = area_triangles(motifs, df_filtered_cubes);

            if weighted_by != "noweight"
                area_weight = []
                for key in keys(motif_energy)
                    # Used to filter out zeros and very small areas (triangles on the vertical for example)
                    if areas[key] > 1
                        push!(area_weight, areas[key]/motif_energy[key][weight_key])
                    end
                end
            else
                area_weight = []
                for key in keys(motif_energy)
                    # Used to filter out zeros and very small areas (triangles on the vertical for example)
                    if areas[key] > 1
                        push!(area_weight, areas[key]/1)
                    end
                end
            end

            # THE FIT
            x_ccdf_original_data, y_ccdf_original_data = powlaw.ccdf(area_weight)
            
            fit_tsallis = pyeval("""lambda fit: lambda a, b, c, d: fit(a, b, c, d)""")
            @. tsallis_ccdf(x, α, λ, c) = c*((1+x/(λ))^(-α))
            popt_tsallis, pcov_tsallis = so.curve_fit(fit_tsallis((x, α, λ, c)->tsallis_ccdf(x, α, λ, c)), x_ccdf_original_data, y_ccdf_original_data, bounds=(0, Inf), maxfev=3000)
            println("alpha= ", popt_tsallis[1],"\nlambda= ", popt_tsallis[2], "\nc= ", popt_tsallis[3])
        
            alpha = round(popt_tsallis[1], digits=4)
            lambda = round(popt_tsallis[2], digits=4)
            c = round(popt_tsallis[3], digits=4)

            set_theme!(Theme(fonts=(; regular="CMU Serif")))

            markers=[:utriangle, :circle]
            colors=[:lightblue, :lightgreen]
            line_colors=[:midnightblue, :green]
            i=1

            ########################################### ALL
            # CCDF of all data scattered 
            fig = Figure(resolution = (700, 600), font= "CMU Serif") 
            ax1 = Axis(fig[1,1], xlabel = L"A_{TE}", ylabel = L"P_A", xscale=log10, yscale=log10, ylabelsize = 28,
                xlabelsize = 28, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
                xticksize = 5, ytickalign = 1, yticksize = 5 , xlabelpadding = 10, ylabelpadding = 10, xticklabelsize=22, yticklabelsize=22)
    
            sc1 = scatter!(ax1, x_ccdf_original_data, y_ccdf_original_data,
                color=(colors[i], 0.75), strokewidth=0.05, strokecolor=(line_colors[i], 0.8), marker=markers[i], markersize=12)

            lines!(ax1, x_ccdf_original_data, tsallis_ccdf(x_ccdf_original_data, popt_tsallis[1], popt_tsallis[2], popt_tsallis[3]), label= L"\alpha=%$(alpha),\, \lambda=%$(lambda),\, c=%$(c)",
                color=(line_colors[i], 0.7), linewidth=2.5)

            # AXIS LEGEND
            axislegend(ax1, [sc1], [L"L=%$(cell_size), \, M_\mathrm{min}=%$(minimum_magnitude)"],
            position = :rt, bgcolor = (:grey90, 0.25), labelsize=18, titlesize=18);

            axislegend(ax1, position = :lb, bgcolor = (:grey90, 0.25), labelsize=18);

            save("./motifs/tsallis/$weighted_by/$region/motif$(motif)_$(region)_cell_size_$(string(cell_size))km_minmag_$(string(minimum_magnitude))_area_$weighted_by.png", fig, px_per_unit=7)
        end
    end
end

# Tetrahedron Analysis
function analize_motifs_tetrahedron_tsallis(region, weighted_by)
    # Read data
    path = "./data/"
    filepath = path * region * ".csv"
    df = CSV.read(filepath, DataFrame);

    motif="Tetrahedron"
    if weighted_by == "totalenergy"
        weight_key = 1
    else 
        weight_key = 2
    end

    # Make path for results
    mkpath("./motifs/$weighted_by/$region")

    # Based on parameter dependency, extract which cell_size lengths are the best:
    if region == "Romania"
        cell_sizes = [3.0, 3.5, 4.0, 4.5, 5.0];
        minimum_magnitudes = [0,1,2,3];
    elseif region == "Italy"
        cell_sizes = [4.0, 4.5, 5.0, 5.5, 6.0];
        minimum_magnitudes = [0,1, 2,3];
    elseif region == "California"
        cell_sizes = [1.0, 1.5, 2.0];
        minimum_magnitudes = [0,1,2,3];
    elseif region == "Japan"
        cell_sizes = [2.5, 3.0, 3.5, 4.0, 5.0];
        minimum_magnitudes = [2,3,4];
    end;

    for cell_size in cell_sizes
        for minimum_magnitude in minimum_magnitudes

            df_filtered = df[df.Magnitude .> minimum_magnitude,:] 

            # Split into cubes
            df_filtered, df_filtered_cubes = region_cube_split(df_filtered,cell_size=cell_size,energyRelease=true);

            # Get the motif
            network_target_path = "./networks/$(region)/cell_size_$(string(cell_size))km/"
            motif_filename = "motif$(motif)_$(region)_cell_size_$(string(cell_size))km_minmag_$(string(minimum_magnitude)).csv"

            # motifs = CSV.read(network_target_path * motif_filename, DataFrame);
            motifs = readdlm(network_target_path * motif_filename, ',', Int64);

            # Energy and volumes calculator
            motif_energy = total_mean_energy(motifs, df_filtered, df_filtered_cubes);
            volumes = volume_tetrahedrons(motifs, df_filtered_cubes);

            # Volumes weighted by total/mean energy
            if weighted_by != "noweight"
                # Volumes weighted by total/mean energy
                volume_weight = []
                for key in keys(motif_energy)
                    # Used to filter out zeros and very small volumes (triangles on the vertical for example)
                    if volumes[key] > 1
                        push!(volume_weight, volumes[key]/motif_energy[key][weight_key])
                    end
                end
            else
                # Volumes weighted by total/mean energy
                volume_weight = []
                for key in keys(motif_energy)
                    # Used to filter out zeros and very small volumes (triangles on the vertical for example)
                    if volumes[key] > 1
                        push!(volume_weight, volumes[key]/1)
                    end
                end
            end

            # THE FIT
            x_ccdf_original_data, y_ccdf_original_data = powlaw.ccdf(volume_weight)
            
            fit_tsallis = pyeval("""lambda fit: lambda a, b, c, d: fit(a, b, c, d)""")
            @. tsallis_ccdf(x, α, λ, c) = c*((1+x/(λ))^(-α))
            popt_tsallis, pcov_tsallis = so.curve_fit(fit_tsallis((x, α, λ, c)->tsallis_ccdf(x, α, λ, c)), x_ccdf_original_data, y_ccdf_original_data, bounds=(0, Inf), maxfev=3000)
            println("alpha= ", popt_tsallis[1],"\nlambda= ", popt_tsallis[2], "\nc= ", popt_tsallis[3])
        
            alpha = round(popt_tsallis[1], digits=4)
            lambda = round(popt_tsallis[2], digits=4)
            c = round(popt_tsallis[3], digits=4)

            set_theme!(Theme(fonts=(; regular="CMU Serif")))

            markers=[:utriangle, :circle]
            colors=[:lightblue, :lightgreen]
            line_colors=[:midnightblue, :green]
            i=1

            ########################################### ALL
            # CCDF of all data scattered 
            fig = Figure(resolution = (700, 600), font= "CMU Serif") 
            ax1 = Axis(fig[1,1], xlabel = L"A_{TE}", ylabel = L"P_A", xscale=log10, yscale=log10, ylabelsize = 28,
                xlabelsize = 28, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
                xticksize = 5, ytickalign = 1, yticksize = 5 , xlabelpadding = 10, ylabelpadding = 10, xticklabelsize=22, yticklabelsize=22)
    
            sc1 = scatter!(ax1, x_ccdf_original_data, y_ccdf_original_data,
                color=(colors[i], 0.75), strokewidth=0.05, strokecolor=(line_colors[i], 0.8), marker=markers[i], markersize=12)

            lines!(ax1, x_ccdf_original_data, tsallis_ccdf(x_ccdf_original_data, popt_tsallis[1], popt_tsallis[2], popt_tsallis[3]), label= L"\alpha=%$(alpha),\, \lambda=%$(lambda),\, c=%$(c)",
                color=(line_colors[i], 0.7), linewidth=2.5)

            # AXIS LEGEND
            axislegend(ax1, [sc1], [L"L=%$(cell_size), \, M_\mathrm{min}=%$(minimum_magnitude)"],
            position = :rt, bgcolor = (:grey90, 0.25), labelsize=18, titlesize=18);

            axislegend(ax1, position = :lb, bgcolor = (:grey90, 0.25), labelsize=18);

            save("./motifs/tsallis/$weighted_by/$region/motif$(motif)_$(region)_cell_size_$(string(cell_size))km_minmag_$(string(minimum_magnitude))_volume_$weighted_by.png", fig, px_per_unit=7)
    
        end
    end
end



analize_motifs_triangle_tsallis(ARGS[1],ARGS[2])
analize_motifs_tetrahedron_tsallis(ARGS[1],ARGS[2])


# ARGS[1]

# for region in ["Romania", "Italy", "California", "Japan"]
#     analize_motifs_triangle_tsallis(region, "totalenergy")
#     analize_motifs_tetrahedron_tsallis(region, "totalenergy")

#     analize_motifs_triangle_tsallis(region, "meanenergy")
#     analize_motifs_tetrahedron_tsallis(region, "meanenergy")
# end