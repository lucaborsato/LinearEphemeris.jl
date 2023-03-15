# LinearEphemeris

See [Documentation](https://lucaborsato.github.io/LinearEphemeris.jl/).

## 1. Introduction  

A [Julia](https://julialang.org/) package for linear ephemerides fit of exoplanet transit times.  
It uses the Weighted Least Square (WLS) fit based on the [Numerical Recipes](https://ui.adsabs.harvard.edu/abs/1992nrfa.book.....P/abstract).  
It can evaluate the errors on the reference time (`Tref`, that is the intercept) and on the period (`Pref`, that is the slope) based on formal error of the WLS or with a bootstrap approach.  
It returns to the standard output a summary of the fit (with statistics) and each transit time (`T0`) with own corresponding `epoch` (or transit number), expected transit time (`Tlin`) from the linear ephem. and `O-C` (Observed - Calculated, where `O == T0` and `C == Tlin`), and the source (that is the telescope used to observe that transit).  

## 2. Install  

```
julia> using Pkg
julia> Pkg.clone("to/add/path/to/github/here")
```

## 3. Example  

Create/Read list of `T0s = [1234.9034, ..., 16789.3980]`, errors `err_T0s = [0.00034, ..., 0.00045]`.  
Define input `Tref`, `Pref`.  
Using default parameters, run as:  

```
using LinearEphemeris
full_linear_ephemeris_analysis(
    T0s,
    err_T0s;
    Tref_in = Tref,
    Pref_in = Pref
)
```
It will output to screen a summary of the fit with statistics, bootstrap error analysis, and it will open a plot gui.  

You can provide also a list of `sources = ["TESS", ..., "K2"] and it will plot each telescope/source with a different color.  

The following is a call of the function with the allowed parameters and default values:  

```
full_linear_ephemeris_analysis(
        T0s,
        err_T0s;
        Tref_in = nothing,
        Pref_in = nothing,
        tscale = nothing,
        sources = nothing,
        bootstrap = true,
        nboot = nothing,
        do_plot = true,
        show_gui = true,
        plot_file = nothing,
        seed = 42,
    )
```

where:  
`tscale` is a parameter to be subtracted to the x-axis, only for plot.  
`nboot` is the number of bootstrap iterations, if `nothing` it will computed automatically as `n (log(n))^2`, with `n` the number of `T0s`.  
`do_plot` a keyword to do the plot or not.  
`plot_file` is the file name (with complete path) of the output of the plot, if `do_plot = false` it will not be used.  
`seed` is the seed for the random number generator.  

There is a working example based on [Kepler- 9](https://ui.adsabs.harvard.edu/abs/2019MNRAS.484.3233B/abstract) in the [examples](../../examples/) folder, 
just run [Kepler-9.jl](../../examples/Kepler-9.jl) and it will print out the results and create the two O-C diagrams (one for each planet).  

---