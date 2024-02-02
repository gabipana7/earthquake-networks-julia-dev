using StatsBase, LinearAlgebra


function histogram_fit(fit_data)

    bstep = 2 * iqr(fit_data) * length(fit_data) ^ (-1/3)
    b = minimum(fit_data) : bstep : maximum(fit_data)
    h = StatsBase.fit(Histogram{Float64}, fit_data, b)


    # Collect bin edges and calculate middle between each 2 points
    xx = collect(h.edges[1])
    x=[]
    for i in eachindex(xx[1:end-1])
        push!(x,(xx[i]+xx[i+1])/2)
    end
    # Results x, the middle of each bin

    # Collect bin weights (number of counts in each bin)
    y = h.weights

    # Code for stoping at first zero (apply to both vectors)
    y_nozero=[]
    x_nozero=[]
    for i in eachindex(y)
        if y[i] == 0.0
            y_nozero = y[1:i-1]
            x_nozero = x[1:i-1]
            break
        end
    end

    # Normalize the bin weights (so they add up to 1)
    y_nozero_norm = LinearAlgebra.normalize(y_nozero);

    return h, x_nozero,y_nozero_norm
end


function power_law(x,a,b)
    return a*x .^(b)
end