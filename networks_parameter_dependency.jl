using CSV, DataFrames
using PyCall
using CairoMakie

include("./src/cubes.jl")
include("./src/network.jl")

@pyimport powerlaw as powlaw


###########################################################################################################################
###########################################################################################################################
# Parameter dependency connectivity on cell_size results function
function networks_parameter_dependency(region; magnitude_threshold=0.0)
    # Read data
    path = "./data/"
    filepath = path * region * ".csv"
    df = CSV.read(filepath, DataFrame);

    # Make path for results
    mkpath("./results/$region")

    # # Magnitude Threshold if you need it
    # magnitude_threshold = 0.0
    # df = df[df.Magnitude .> magnitude_threshold,:];

    # Cell sizes ranges from 0.5 to 20 with a 0.5 increment
    cube_cell_sizes= range(0.5, 20, step=0.5)

    degrees_alpha=[]
    degrees_sigma = []
    degrees_xmin=[]
    marker_size=[]
    

    for cell_size in cube_cell_sizes
        df, df_cubes = region_cube_split(df,cell_size=cell_size)
        MG = create_network(df, df_cubes)

        degrees=[]
        for i in 1:nv(MG)
            push!(degrees, get_prop(MG, i, :degree))
        end

        fit = powlaw.Fit(degrees);
        push!(degrees_alpha, fit.alpha)
        push!(degrees_sigma, fit.sigma)
        push!(degrees_xmin, fit.xmin)
        push!(marker_size, fit.power_law.KS(data=degrees))

    end

    results = DataFrame([cube_cell_sizes, degrees_alpha, degrees_sigma, degrees_xmin, marker_size], ["cell_size", "alpha", "sigma", "xmin", "KS"])
    CSV.write("./results/$region/$(region)_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size.csv", results, delim=",", header=true);

end


###########################################################################################################################


# Parameter dependency connectivity on cell_size plot function
function networks_parameter_dependency_plot_cairo(region; magnitude_threshold=0.0, goodness_of_fit=true)
    results = CSV.read("./results/$region/$(region)_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size.csv", DataFrame)

    if goodness_of_fit == true
        set_theme!(Theme(fonts=(; regular="CMU Serif")))
        fig = Figure(resolution = (1100, 400), font= "CMU Serif") ## probably you need to install this font in your system
        xlabels = [L"\text{cell\,size\,[km]}", L"\text{cell\,size\,[km]}"]
        ylabels = [L"\alpha", L"x_{min}"]
        ax = [Axis(fig[1, i], xlabel = xlabels[i], ylabel = ylabels[i], ylabelsize = 26,
            xlabelsize = 22, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
            xticksize = 5, ytickalign = 1, yticksize = 5 , xlabelpadding = 10, ylabelpadding = 10) for i in 1:2]
        scatter!(ax[1],results.cell_size, results.alpha; color = results.KS, colormap = reverse(cgrad(:redblue)),
            markersize = 13, marker = :circle, strokewidth = 0)
        xmin_scat = scatter!(ax[2],results.cell_size, results.xmin; color = results.KS, colormap = reverse(cgrad(:redblue)),
            markersize = 13, marker = :circle, strokewidth = 0)
        cbar = Colorbar(fig[1,3], xmin_scat, label="KS distance")
        save("./results/$region/$(region)_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size_goodness_fit.png", fig, px_per_unit=5)
    end

    set_theme!(Theme(fonts=(; regular="CMU Serif")))
    fig = Figure(resolution = (900, 400), font= "CMU Serif") ## probably you need to install this font in your system
    xlabels = [L"\text{cell\,size\,[km]}", L"\text{cell\,size\,[km]}"]
    ylabels = [L"\alpha", L"x_{min}"]
    ax = [Axis(fig[1, i], xlabel = xlabels[i], ylabel = ylabels[i], ylabelsize = 26,
        xlabelsize = 22, xgridstyle = :dash, ygridstyle = :dash, xtickalign = 1,
        xticksize = 5, ytickalign = 1, yticksize = 5 , xlabelpadding = 10, ylabelpadding = 10) for i in 1:2]
    scatter!(ax[1],results.cell_size, results.alpha;
        markersize = 13, marker = :circle, strokewidth = 0)
    scatter!(ax[2],results.cell_size, results.xmin;
        markersize = 13, marker = :circle, strokewidth = 0)
    save("./results/$region/$(region)_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size.png", fig, px_per_unit=5)
end

###########################################################################################################################
###########################################################################################################################

region_list = ["Romania", "Italy", "California", "Japan"]

for region in region_list
    networks_parameter_dependency(region)
    # networks_parameter_dependency_plot_cairo(region)
end




###########################################################################################################################
###########################################################################################################################


# Analysis of parameter dependency. Based on minimum KS, proper alpha and minimum xmin
function parameter_dependency_analysis(region; magnitude_threshold=0.0)
    results = CSV.read("./results/$region/$(region)_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size.csv", DataFrame);

    # Eliminate records if alpha is not in range [1.5,3.5]

    results_filter_alpha = results[(results.alpha .> 1.7) .& (results.alpha .< 3.2), :];

    # Sort by KS and keep first 10 only
    results_sorted_KS = sort!(results_filter_alpha, [:KS])
    results_filter_alpha_sorted_KS = first(results_sorted_KS, 10)

    # Sort by xmin and keep first 7 only
    results_filter_alpha_sorted_KS_sorted_xmin = sort!(results_filter_alpha_sorted_KS, [:xmin])
    best_fits_region = first(results_filter_alpha_sorted_KS_sorted_xmin, 7) 
    return(best_fits_region)
end


###########################################################################################################################

magnitude_threshold = 0.0
df = DataFrame([[0,0,0,0,0,0,0]],["par_dep"])
for region in ["romania","california","italy","japan"]
    par_dep_best_fits = DataFrame([[0,0,0,0,0,0,0]],[region])
    best_fits = parameter_dependency_analysis(region; magnitude_threshold)
    par_dep_best_fits = hcat(par_dep_best_fits, best_fits)
    df = hcat(df, par_dep_best_fits, makeunique=true)
end

CSV.write("./results/best_fits_all_regions_minmag_$(magnitude_threshold)_alpha_xmin_dependency_cell_size.csv", df, delim=",", header=true);

###########################################################################################################################

