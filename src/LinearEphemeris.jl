module LinearEphemeris

using ColorSchemes
using LaTeXStrings
using LinearAlgebra
using Plots, Plots.Measures
using Printf
using Random
using Statistics
using StatsBase

# pyplot()
# gr()

# exports
export BTJD
export DAY2SEC, SEC2DAY, DAY2MIN, MIN2DAY, DAY2HOUR, HOUR2DAY
export linear_error_prop, transit_epoch, linear_transit_time, linear_transit_times
export wls_nr
export fit_linear_ephemeris, statistics_fit, get_nboot, classical_bootstrap
export full_linear_ephemeris_analysis

BTJD = 2457000.0

DAY2SEC = 86400.0
SEC2DAY = 1.0 / DAY2SEC
DAY2MIN = 1440.0
MIN2DAY = 1.0 / DAY2MIN
DAY2HOUR = 24.0
HOUR2DAY = 1.0 / DAY2HOUR


function linfunc(q, m, x)
    return q + m * x
end

"""
    linear_error_prop(eq, em, x)

propagates error on x value as:
```math
\\sqrt{\\mathrm{eq}^2 + (x × \\mathrm{em})^2}
```
# Arguments:
- `eq::Float64`: error on intercept (reference time Tref)
- `em::Float64`: error on slope (reference period Pref)
- `x::Float64`: value at which compute the error
"""
function linear_error_prop(eq, em, x)
    return sqrt(eq^2 + (x * em)^2)
end


"""
    transit_epoch(T0, Tref, Pref)

computes the transit epoch (or transit number N) with respect to a linear ephemeris:
```math
    N =  \\mathrm{round}( (T_{0,\\mathrm{lin}} - T_\\mathrm{ref}) / P_\\mathrm{ref})  
```
"""
function transit_epoch(T0, Tref, Pref)
    return round((T0 - Tref) / Pref)
end


"""
    linear_transit_time(Tref, Pref, epoch)

computes the predicted transit time at specific epoch (transit number N) given a linear ephemeris:  
```math
    T_{0,\\mathrm{lin}, N} = T_\\mathrm{ref} + N \\times P_\\mathrm{ref}  
```
"""
function linear_transit_time(Tref, Pref, epoch)
    return linfunc(Tref, Pref, epoch)
end

"""
    linear_transit_times(Tref, Pref, epochs) 

as `linear_transit_time` but for vector/list of epochs

"""
function linear_transit_times(Tref, Pref, epochs)
    return linear_transit_time.(Tref, Pref, epochs)
end

# """
#     wls_nr(x, y, ey)

# WEIGHTED LEAST SQUARE AS IN [NUMERICAL RECIPES](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)

# # Arguments:
# - `x::Vector{Float64}`: independent variable  
# - `y::Vector{Float64}`: dependent variable  
# - `ey::Vector{Float64}`: measurement errors on dependent variable
# # Returns:
# - `intercept::tuple(fit::Float64, err::Float64)`: intercept with value (fit) and formal error (err)
# - `slope::tuple(fit::Float64, err::Float64)`: slope with value (fit) and formal error (err)
# """
# function wls_nr(x, y, ey)
#     # WEIGHTED LEAST SQUARE AS IN NUMERICAL RECIPES
#     # return intercept q and slope m of a straight linear_ephemeris
#     # of type y = q + m * x
#     #
#     # x == independent variable
#     # y == dependent variable
#     # ey == measurement errors on y
#     #
#     # q = tuple(fit, err)
#     # m = tuple(fit, err)

#     w = 1.0 ./ (ey .^ 2)
#     S = sum(w)
#     Sx = sum(w .* x)
#     Sy = sum(w .* y)
#     Sxx = sum(w .* (x .^ 2))
#     # Syy = sum( w .* (y .^2))
#     Sxy = sum(w .* x .* y)

#     Δ = S * Sxx - (Sx^2)
#     q = (Sxx * Sy - Sx * Sxy) / Δ
#     m = (S * Sxy - Sx * Sy) / Δ

#     err_q = sqrt(Sxx / Δ)
#     err_m = sqrt(S / Δ)

#     return (fit = q, err = err_q), (fit = m, err = err_m)
# end

"""
    wls_nr(x, y, ey)

WEIGHTED LEAST SQUARE AS IN [NUMERICAL RECIPES](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)

# Arguments:
- `x::Vector{Float64}`: independent variable  
- `y::Vector{Float64}`: dependent variable  
- `ey::Vector{Float64}`: measurement errors on dependent variable
# Returns:
- `intercept::tuple(fit::Float64, err::Float64)`: intercept with value (fit) and formal error (err)
- `slope::tuple(fit::Float64, err::Float64)`: slope with value (fit) and formal error (err)
"""
function wls_nr(x, y, ey)
    # WEIGHTED LEAST SQUARE AS IN NUMERICAL RECIPES
    # return intercept q and slope m of a straight linear_ephemeris
    # of type y = q + m * x
    #
    # x == independent variable
    # y == dependent variable
    # ey == measurement errors on y
    #
    # q = tuple(fit, err)
    # m = tuple(fit, err)

    # nx = len(x)
    w = 1.0 ./ (ey .^ 2)
    S = sum(w)
    Sx = dot(w, x)
    Sy = dot(w, y)
    
    t = (x .- (Sx/S)) ./ ey
    m = dot(t./ey, y)
    St2 = dot(t, t)

    m /= St2
    q = (Sy - (Sx*m))/S
    err_q = sqrt((1.0 + ((Sx*Sx)/(S*St2)))/S)
    err_m = sqrt(1.0/St2)

    return (fit = q, err = err_q), (fit = m, err = err_m)
end

"""
    wls_nr(x, y, ey)

WEIGHTED LEAST SQUARE AS IN [NUMERICAL RECIPES](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract)

# Arguments:
- `x::Vector{Float64}`: independent variable  
- `y::Vector{Float64}`: dependent variable  
- `ey::Vector{Float64}`: measurement errors on dependent variable
# Returns:
- `intercept::tuple(fit::Float64, err::Float64)`: intercept with value (fit) and formal error (err)
- `slope::tuple(fit::Float64, err::Float64)`: slope with value (fit) and formal error (err)
"""
function wls_nr(x, y)
    # WEIGHTED LEAST SQUARE AS IN NUMERICAL RECIPES
    # return intercept q and slope m of a straight linear_ephemeris
    # of type y = q + m * x
    #
    # x == independent variable
    # y == dependent variable
    #
    # q = tuple(fit, err)
    # m = tuple(fit, err)

    nx = len(x)
    S = nx
    Sx = sum(x)
    Sy = sum(y)

    t = x .- (Sx/S)
    m = dot(t, y)
    St2 = dot(t, t)
    
    m /= St2
    q = (Sy - (Sx*m))/S
    err_q = sqrt((1.0 + ((Sx*Sx)/(S*St2)))/S)
    err_m = sqrt(1.0/St2)

    t = y .- (q + m .* x)
    chi2 = dot(t,t)
    dof = nx-2
    sigdat = sqrt(chi2/dof)
    err_q = err_q*sigdat
    err_m = err_m*sigdat

    return (fit = q, err = err_q), (fit = m, err = err_m)
end


"""
    fit_linear_ephemeris(epochs, T0s, err_T0s)

as in `wls_nr` but with meaningful variable names
# Arguments:
- `epochs::Vector{Float64}`: transit epochs/numbers  
- `T0s::Vector{Float64}`: transit times  
- `err_T0s::Vector{Float64}`: measurement errors on T0s
# Returns:
- `Tref::Tuple(fit::Float64, err::Float64)`: reference time of the linear ephemeris with value (fit) and formal error (err)
- `Pref::Tuple(fit::Float64, err::Float64)`: period of the linear ephemeris with value (fit) and formal error (err)
"""
function fit_linear_ephemeris(epochs, T0s, err_T0s)

    Tref, Pref = wls_nr(epochs, T0s, err_T0s)

    return Tref, Pref
end

"""
    statistics_fit(y_obs, y_mod; ey_obs = nothing, n_fit = 2)

computes some statistics on the fitted linear ephemeris
# Arguments:
- `y_obs::Vector{Float64}` observed values (T0s)
- `y_mod::Vector{Float64}` model values (T0s,lin)
- `ey_obs::Vector{Float64} = nothing` errors on the observed values (if available)
- `n_fit::Int = 2` number of fitted parameters of the model (2 for the linear ephemeris)
# Returns:
- `stats::NamedTuple`: named tuple with various statistics based on the residuals
    - `stats.std::Float64` standard deviation of the residuals  
    - `stats.p68::Float64` 68.27th percentile of the absolute residuals
    - `stats.chi_square::Float64` `χ^2` like `= ∑((obs - mod) / err_obs)^2`
    - `stats.red_chi_square::Float64` `χ^2 / dof`
    - `stats.dof::Float64` degrees of freedom `dof = n_data-n_fit`
    - `stats.bic::Float64` [BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion)`= χ^2 + n_fit ln(n_data)`
    - `stats.aic::Float64` [AIC](https://en.wikipedia.org/wiki/Akaike_information_criterion)`= χ^2 + 2 n_fit`
    - `stats.lnc::Float64` log-Likelihood constant part given by the `2π` and `σ_i` terms
    - `stats.lnL::Float64` log-Likelihood
    - `stats.bic_lnL::Float64` [BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion)`= -2lnL + n_fit ln(n_data)`
    - `stats.aic_lnL::Float64` [AIC](https://en.wikipedia.org/wiki/Akaike_information_criterion)`= -2lnL + 2 n_fit`

"""
function statistics_fit(y_obs, y_mod; ey_obs = nothing, n_fit = 2)

    n_data = length(y_obs)
    dof = n_data - n_fit
    if isnothing(ey_obs)
        ey_obs = ones(n_data)
    end
    res = y_obs .- y_mod
    std = Statistics.std(res)
    p68 = StatsBase.percentile(abs.(res), 68.27)
    wres = res ./ ey_obs
    chi_square = sum(wres .^ 2)
    red_chi_square = chi_square / dof
    bic = chi_square + n_fit * log(n_data)
    aic = chi_square + 2.0 * n_fit
    lnc = n_data * (n_fit * log(2.0 * pi) + log(sum(ey_obs .^ 2)))
    lnL = -0.5 * (lnc + chi_square)
    bic_lnL = -2.0 * lnL + n_fit * log(n_data)
    aic_lnL = -2.0 * lnL + 2.0 * n_fit
    stats = (
        std = std,
        p68 = p68,
        chi_square = chi_square,
        red_chi_square = red_chi_square,
        dof = dof,
        bic = bic,
        aic = aic,
        lnc = lnc,
        lnL = lnL,
        bic_lnL = bic_lnL,
        aic_lnL = aic_lnL,
    )

    return stats
end


"""
    get_scaling_oc(Ad)

compute the scaling value and unit of the O-C diagram
# Arguments:
- `Ad::Float64` TTV amplitude or std or rms
# Returns:
- `scale::NamedTuple` access it with `.val` and `.lab` for value and unit, respectively
"""
function get_scaling_oc(Ad)

    Ah = Ad * DAY2HOUR
    Am = Ad * DAY2MIN
    As = Ad * DAY2SEC
    if As < 60.0
        u_v, u_s = DAY2SEC, "s"
    elseif 1.0 <= Am < 60.0
        u_v, u_s = DAY2MIN, "min"
    elseif 1.0 <= Ah < 24.0
        u_v, u_s = DAY2HOUR, "hours"
    else
        u_v, u_s = 1.0, "d"
    end

    scale = (val = u_v, lab = u_s)
    return scale
end


"""
    get_nboot(n)

Computes the minimum proposed number of bootstrap iteration as:
```math
    n_\\mathrm{boot} = n \\times \\mathrm{round}(Int, \\ln(n)^2)
```
# Arguments:
- `n::Int` number of observed transit times (T0s)
# Returns:
- `n_boot::Int` number of bootstrap iterations
"""
function get_nboot(n)
    return n * round(Int, log(n)^2)
end


"""
    classical_bootstrap(
        x,
        y, 
        q, 
        m; 
        ey_in = nothing, 
        nboot = nothing, 
        rng=Random.GLOBAL_RNG,
        return_distribution=true
    )

# Arguments:
- `x::Vector{Float64}`: independent variable  
- `y::Vector{Float64}`: dependent variable  
- `q::Float`: intercept (the true or best-fit value)  
- `m::Float`: slope (the true or best-fit value)  
- `ey_in::Vector{Float64} = nothing`: measurement errors on dependent variable (default is nothing)  
- `nboot::Int = nothing`: number of bootstrap iterations to run  
- `rng = Random.GLOBAL_RNG`: random number generator  
- `return_distribution::Bool = true`: return or not the bootstrap distribution of each parameter  
# Return:
- `tuple(median::Float64, rms=Float64, distribution=Vector{Float64})`: intercept bootstrap (median, rms, distribution); `return_distribution=false` will set `distribution=nothing`  
- `tuple(median::Float64, rms=Float64, distribution=Vector{Float64})`: slop bootstrap

"""
function classical_bootstrap(
    x,
    y,
    q,
    m;
    ey_in = nothing,
    nboot = nothing,
    rng = Random.GLOBAL_RNG,
    return_distribution = true,
)

    y_ok = linfunc.(q, m, x)
    r_ok = y .- y_ok
    scale = mean(r_ok)
    r_ok .-= scale
    n = length(x)
    if isnothing(nboot) || nboot < 1
        nboot = get_nboot(n)
    end

    if isnothing(ey_in)
        ey = ones(n)
    else
        ey = ey_in
    end

    q_boot = zeros(nboot)
    m_boot = zeros(nboot)

    for i_boot = 1:nboot
        idx = sample(rng, 1:n, n, replace = true)
        y_boot = y_ok + r_ok[idx]
        ey_boot = ey[idx]
        qx, mx = wls_nr(x, y_boot, ey_boot)
        q_boot[i_boot] = qx[1]
        m_boot[i_boot] = mx[1]
    end

    q_med = median(q_boot)
    q_rms = percentile(abs.(q_boot .- q), 68.27)

    m_med = median(m_boot)
    m_rms = percentile(abs.(m_boot .- m), 68.27)

    if return_distribution
        return (median = q_med, rms = q_rms, distribution = q_boot),
        (median = m_med, rms = m_rms, distribution = m_boot)
    else
        return (median = q_med, rms = q_rms, distribution = nothing),
        (median = m_med, rms = m_rms, distribution = nothing)
    end
end

"""
    custom_errorbar!(plt, x, y, ey; color = :auto, lw = 1.0, ls = :auto)

Create errorbars on plot withouth the cap.
"""
function custom_errorbar!(plt, x, y, ey; color = :auto, lw = 1.0, ls = :auto)

    for i in eachindex(x)
        xi = x[i]
        yl = y[i] - ey[i]
        yu = y[i] + ey[i]
        plot!(plt, [xi, xi], [yl, yu], color = color, lw = lw, ls = ls, label = "")
    end
    return
end


function plot_one_bootstrap_distribution(x_true, x_boot, x_label)

    println("plotting ", x_label)
    println("best   ", x_true)
    println("median ", x_boot.median)
    println("rms    ", x_boot.rms)
    println("nboot  ", length(x_boot.distribution))

    color = palette(:tab10)[1]
    x_plt = histogram(
        x_boot.distribution,
        bins=:auto,
        normalize=:pdf,
        color=:gray,
        linewidth=0.1,
        linecolor=:white,
        label="bootstrap"
    )

    vline!(
        x_plt,
        [x_true],
        linestyle = :solid,
        linecolor = :black,
        lw = 0.9,
        label = "best",
    )
    vline!(
        x_plt,
        [x_boot.median],
        linestyle = :solid,
        linecolor = color,
        lw = 0.9,
        label = "median",
    )
    vline!(
        x_plt,
        [x_true-x_boot.rms],
        linestyle = :dash,
        linecolor = color,
        lw = 0.8,
        label = L"\mathrm{best}\pm\mathrm{rms}",
    )
    vline!(
        x_plt,
        [x_true+x_boot.rms],
        linestyle = :dash,
        linecolor = color,
        lw = 0.8,
        label = "",
    )

    # xlabel!(x_plt, L"$T_\mathrm{ref}$")
    xlabel!(x_plt, x_label)
   
    return x_plt
end
    

function plot_full_bootstrap_distribution(q_true, q_boot, m_true, m_boot)

    q_plt = plot_one_bootstrap_distribution(q_true, q_boot, L"$T_\mathrm{ref}$")
    m_plt = plot_one_bootstrap_distribution(m_true, m_boot, L"$P_\mathrm{ref}$")

    plt = plot(
        q_plt, m_plt,
        layout=(2,1),
        grid = false,
        size = (1024, 1024),
        dpi = 300,
        thickness_scaling = 1,
        fontfamily = "Computer Modern",
        tickfontsize = 14,
        guidefontsize = 16,
        margin = 5mm,
    )

    return plt
end


"""
    full_linear_ephemeris_analysis(
        T0s,
        err_T0s;
        Tref_in = nothing,
        Pref_in = nothing,
        tscale = nothing,
        sources = nothing,
        bootstrap = true,
        nboot = nothing,
        return_distribution = true,
        do_plot = true,
        show_gui = true,
        plot_file = nothing,
        seed = 42,
    )

Compute the full linear ephemeris with all needed arguments.  
It prints the output on the standard output and do the plots if required.

# Arguments:
- `T0s::Vector{Float64}`: observed transit times vector  
- `err_T0s::Vector{Float64}`: measurement error vector on observed transit times  
- `Tref_in::Float64 = nothing`: input reference time for the linear ephemeris  
- `Pref_in::Float64 = nothing`: input reference period for the linear ephemeris  
- `tscale::Float64 = nothing`: scale time, just for plot  
- `sources::Vector{String} = nothing`: list/vector of sources of each transit time  
- `bootstrap::Bool = true`: run or not the bootstrap analysis  
- `nboot::Int = nothing`: number of bootstrap iterations if required  
- `return_distribution = true`: in case of bootstrap, return distribution (if also `do_plot`, it plots distribution)  
- `do_plot::Bool = true`: do or not the plot  
- `show_gui::Bool = true`: if `do_plot`, this flag open the GUI or not  
- `plot_file::String = nothing`: if `do_plot`, this is the full name of the output plot file  
- `seed::Int = 42`: seed of the random number  

# Returns:
nothing

"""
function full_linear_ephemeris_analysis(
    T0s,
    err_T0s;
    Tref_in = nothing,
    Pref_in = nothing,
    tscale = nothing,
    sources = nothing,
    bootstrap = true,
    nboot = nothing,
    return_distribution = true,
    do_plot = true,
    show_gui = true,
    plot_file = nothing,
    seed = 42,
)

    nT0s = length(T0s)
    if isnothing(Tref_in)
        Tref_in = T0s[nT0s÷2]
    end
    if isnothing(Pref_in)
        Pref_in = median(diff(T0s))
    end

    epo_in = transit_epoch.(T0s, Tref_in, Pref_in)

    Tref, Pref = fit_linear_ephemeris(epo_in, T0s, err_T0s)

    epochs = transit_epoch.(T0s, Tref.fit, Pref.fit)
    Tlin = linear_transit_time.(Tref.fit, Pref.fit, epochs)

    stats = statistics_fit(T0s, Tlin; ey_obs = err_T0s, n_fit = 2)
    oc_d = T0s .- Tlin
    Aoc_d = 0.5 * (maximum(oc_d) - minimum(oc_d))
    scl = get_scaling_oc(Aoc_d)

    if isnothing(sources)
        sources = ["" for i = 1:nT0s]
    end

    if isnothing(tscale) || tscale == 0.0
        tscale = 0.0
        xlscale = ""
    else
        xlscale = @sprintf(" - %.3f", tscale)
    end

    x_plt = Tlin .- tscale
    y_plt = oc_d .* scl.val
    ey_plt = err_T0s .* scl.val

    rng = MersenneTwister(seed)
    Random.seed!(rng, seed)

    println()
    println("--- Linear Ephemeris ---")
    println("WLS fit")
    println()
    @printf("Tref = %.7f +/- %.7f\n", Tref.fit, Tref.err)
    @printf("Pref = %.7f +/- %.7f\n", Pref.fit, Pref.err)
    println()
    println("stats:")
    println("n_data = $(nT0s) n_fit = 2")
    println("std(res) = $(stats.std) d = $(stats.std*scl.val) $(scl.lab)")
    println("percentile(|res|, 68.27th) = $(stats.p68) d = $(stats.p68*scl.val) $(scl.lab)")
    println(
        "χ^2 = $(stats.chi_square) χ^2/dof = $(stats.red_chi_square) dof = $(stats.dof)",
    )
    println("BIC = $(stats.bic) AIC = $(stats.aic)")
    println("lnL = $(stats.lnL) lnc = $(stats.lnc)")
    println("BIC_lnL = $(stats.bic_lnL) AIC_lnL = $(stats.aic_lnL)")
    println()
    @printf("TTV amp = %.7f d = %.7f %s\n", Aoc_d, Aoc_d * scl.val, scl.lab)
    println()

    if bootstrap
        if isnothing(nboot) || nboot < 0
            nboot = get_nboot(nT0s)
        end
        println("Running Bootstrap with nboot = $nboot")
        T_boot, P_boot = classical_bootstrap(
            epochs,
            T0s,
            Tref.fit,
            Pref.fit;
            ey_in = err_T0s,
            nboot = nboot,
            rng = rng,
            return_distribution = return_distribution,
        )
        @printf("Tref (bootstrap): median = %.7f rms = %.7f\n", T_boot.median, T_boot.rms)
        @printf("Pref (bootstrap): median = %.7f rms = %.7f\n", P_boot.median, P_boot.rms)
        println()
    end

    if do_plot

        # create the plot window
        plt = plot(
            grid = false,
            size = (1024, 1024),
            dpi = 300,
            thickness_scaling = 1,
            fontfamily = "Computer Modern",
            xlabel = L"\mathrm{BJD}_\mathrm{TDB}%$(xlscale)",
            ylabel = L"\mathrm{O-C\ (%$(scl.lab))}",
            # xtickfontsize=14, ytickfontsize=14,
            tickfontsize = 14,
            # xguidefontsize=16, yguidefontsize=16,
            guidefontsize = 16,
            margin = 5mm,
        )

        # plot the error propagation with bootstrap or formal error
        if bootstrap
            prop_ey = linear_error_prop.(T_boot.rms, P_boot.rms, epochs) .* scl.val
            println("Error propagation plot with bootstrap analysis")
            if return_distribution
                boot_plt = plot_full_bootstrap_distribution(Tref.fit, T_boot, Pref.fit, P_boot)
            end 
        else
            prop_ey = linear_error_prop.(Tref.err, Pref.err, epochs) .* scl.val
            println("Error propagation plot with formal error")
        end
        plot!(
            plt,
            x_plt,
            zeros(nT0s),
            ribbon = prop_ey,
            color = :gray,
            fillalpha = 0.4,
            lw = 0.0,
            label = "",
        )

        # define colors
        usrc = unique(sources)
        n_src = length(usrc)
        if n_src > 1
            c_src = LinRange(0.9, 0.1, n_src)
        else
            c_src = [0.5]
        end
        colors = ColorSchemes.get(ColorSchemes.nipy_spectral, c_src)
        # plot scatter with errorbars for each unique source
        for (i_s, ss) in enumerate(usrc)
            sel = findall(==(ss), sources)

            custom_errorbar!(
                plt,
                x_plt[sel],
                y_plt[sel],
                ey_plt[sel];
                color = colors[i_s],
                lw = 2.5,
                ls = :solid,
            )
            plot!(
                plt,
                x_plt[sel],
                y_plt[sel],
                seriestype = :scatter,
                label = "$(ss)",
                markersize = 5,
                markercolor = colors[i_s],
                markerstrokecolor = :white,
                markerstrokewidth = 0.5,
            )
        end
        # horizontal line at 0
        hline!(
            plt,
            [0],
            linestyle = :solid,
            linecolor = :black,
            lw = 0.9,
            label = "",
            z_order = 1,
        )
        # vertical line at Tref
        vline!(
            plt,
            [Tref.fit - tscale],
            linestyle = :dash,
            linecolor = :gray,
            lw = 0.5,
            la = 0.7,
            label = "Tref",
            z_order = 1,
        )

        
        # save the file if provided
        if !isnothing(plot_file)
            savefig(plt, plot_file)
            if bootstrap & return_distribution
                savefig(boot_plt, replace(plot_file, ".png" => "_TP_bootstrap.png"))
            end
        end
        if show_gui
            display(plt)
            if bootstrap & return_distribution
                display(boot_plt)
            end
        end

    end

    # TODO:
    # output of T0s and O-C
    println()
    println("# O-C SUMMARY")
    l = @sprintf(
        "# %5s %15s %10s %10s %10s %10s %15s %20s",
        "epoch",
        "T0s",
        "err_T0s",
        "OC_d",
        "OC_m",
        "OC_s",
        "T0s_lin",
        "source"
    )
    println(l)

    for (epo, T0, eT0, oc, lin, ss) in zip(epochs, T0s, err_T0s, oc_d, Tlin, sources)
        l = @sprintf(
            "  %5.0f %15.6f %10.6f %10.6f %10.3f %10.1f %15.6f %20s",
            epo,
            T0,
            eT0,
            oc,
            oc * DAY2MIN,
            oc * DAY2SEC,
            lin,
            ss
        )
        println(l)
    end

    println()

    # for (eT0, ePl) in zip(err_T0s, ey_plt)
    #     println(eT0, " ", ePl)
    # end

    return
end

# =============================================================================
# === EXAMPLE ===

function test()

    # "True" values
    Ptrue = 4.56 # d
    Ttrue = 2003.45 # fake reference time
    # "True" epochs and T0s
    epochs_true = [-3, -2, 0, 1, 3, 5, 6, 7, 8, 11, 12]
    nT0s = length(epochs_true)
    T0s_true = linear_transit_times(Ttrue, Ptrue, epochs_true)

    # define a TTV amplitude
    amp_true = 14.0 * MIN2DAY

    seed = 42
    rng = MersenneTwister(seed)
    Random.seed!(rng, seed)

    # create the observed T0s, err_T0s, and sources
    dT0s = Random.randn(rng, nT0s) .* amp_true
    T0s = T0s_true .+ dT0s
    err_T0s = 15.0 * SEC2DAY .+ rand(rng, nT0s) .* dT0s .* 3.0
    sources = vcat(["K2" for i = 1:4], ["TESS" for i = 5:nT0s])

    # define an initial value for Tref and Pref
    Tref_in = 2003.4563
    Pref_in = 4.5673

    BOOTSTRAP = true
    # nboot=3 * get_nboot(length(T0s))
    nboot = nothing
    do_plot = true
    plot_file = nothing

    full_linear_ephemeris_analysis(
        T0s,
        err_T0s;
        Tref_in = Tref_in,
        Pref_in = Pref_in,
        tscale = nothing,
        sources = sources,
        bootstrap = BOOTSTRAP,
        nboot = nboot,
        do_plot = do_plot,
        plot_file = plot_file,
        seed = seed,
    )

    return
end

end # module
