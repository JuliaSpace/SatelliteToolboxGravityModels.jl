using Documenter
using SatelliteToolboxGravityModels

makedocs(
    modules = [SatelliteToolboxGravityModels],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolboxGravityModels.jl/stable/",
    ),
    sitename = "Satellite Toolbox Gravity Models",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Usage" => "man/usage.md",
        "API" => "man/api.md",
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolboxGravityModels.jl.git",
    target = "build",
)
