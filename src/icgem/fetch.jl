# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Function to fetch and store the ICGEM files.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export fetch_icgem_file

# Dictionary with the pre-configured gravity field models.
const _ONLINE_ICGEM_FILES = Dict(
    :EGM96   => "http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc",
    :EGM2008 => "http://icgem.gfz-potsdam.de/getmodel/gfc/c50128797a9cb62e936337c890e4425f03f0461d7329b09a8cc8561504465340/EGM2008.gfc",
    :JGM2    => "http://icgem.gfz-potsdam.de/getmodel/gfc/291f7d127f49fe3bcb4a633a06e8b50f461d4b13ff8e7f2046644a8148b57fc6/JGM2.gfc",
    :JGM3    => "http://icgem.gfz-potsdam.de/getmodel/gfc/a3375e01a717ac162962138a5e94f10466b71aa4a130d7f7d5b18ab3d5f90c3d/JGM3.gfc"
)

"""
    fetch_icgem_file(url::AbstractString; kwargs...) -> String
    fetch_icgem_file(model::Symbol; kwargs...) -> String

Fetch a ICGEM file from the `url` and return its file path to be parsed with the function
[`GravityModels.load`](@ref). If the file already exists, it will not be re-downloaded
unless the keyword `force = true` is passed.

A symbol can be passed instead the URL to fetch pre-configured gravity field models. The
supported values are:

- `:EGM96`: Earth Gravitational Model from 1996.
- `:EGM2008`: Earth Gravitational Model from 2008.
- `:JGM2`: Joint Gravity Model 2.
- `:JGM3`: Joint Gravity Model 3.

# Examples

```julia-repl
julia> fetch_icgem_file(:EGM96)
[ Info: Downloading the ICGEM file 'EGM96.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'...
"/Users/ronan.arraes/.julia/scratchspaces/bd9e9728-6f7b-4d28-9e50-c765cb1b7c8c/icgem/EGM96.gfc"

julia> fetch_icgem_file(:EGM96)
"/Users/ronan.arraes/.julia/scratchspaces/bd9e9728-6f7b-4d28-9e50-c765cb1b7c8c/icgem/EGM96.gfc"
```
"""
function fetch_icgem_file(model::Symbol; force::Bool = false)
    !haskey(_ONLINE_ICGEM_FILES, model) &&
        throw(ArgumentError("The model $model was not found in the pre-build dictionary."))

    return fetch_icgem_file(_ONLINE_ICGEM_FILES[model]; force = force)
end

function fetch_icgem_file(url::AbstractString; force::Bool = false)
    cache_dir = @get_scratch!("icgem")

    # We must be able to get the file name from the URL.
    filename = basename(url)
    isempty(filename) && throw(ArgumentError("We could not obtain the file name from the URL."))

    filepath = joinpath(cache_dir, filename)

    if isfile(filepath) && !force
        # If the file exists, just return the path.
        return filepath

    else
        # Otherwise, let's fetch the file from the URL.
        @info "Downloading the ICGEM file '$filename' from '$url'..."
        return Downloads.download(url, filepath)
    end
end
