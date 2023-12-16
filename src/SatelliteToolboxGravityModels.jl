module SatelliteToolboxGravityModels

using Crayons
using Dates
using Downloads
using Scratch

import Base: show

############################################################################################
#                                        Submodules
############################################################################################

include("./GravityModels/GravityModels.jl")
using .GravityModels
export GravityModels
export AbstractGravityModel

############################################################################################
#                                          Types
############################################################################################

include("./types.jl")

############################################################################################
#                                        Constants
############################################################################################

const _D = Crayon(reset = true)
const _B = crayon"bold"

############################################################################################
#                                         Includes
############################################################################################

include("./icgem/api.jl")
include("./icgem/compute.jl")
include("./icgem/fetch.jl")
include("./icgem/parse.jl")
include("./icgem/show.jl")

end # module SatelliteToolboxGravityModels
