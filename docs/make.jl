using Documenter
using LinearEphemeris

makedocs(
    sitename = "LinearEphemeris",
    pages = [
        "Home" => "index.md"
    ],
    format = Documenter.HTML(),
    modules = [LinearEphemeris]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/lucaborsato/LinearEphemeris.jl",
    devbranch = "main"
)
