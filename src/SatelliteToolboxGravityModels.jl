module SatelliteToolboxGravityModels

using Crayons
using Dates
using Downloads
using Scratch

import Base: show
import SatelliteToolboxBase: LowerTriangularStorage, RowMajor

############################################################################################
#                                        Submodules                                        #
############################################################################################

include("./GravityModels/GravityModels.jl")
using .GravityModels
import .GravityModels: _from_j2000_seconds, _to_j2000_seconds
export GravityModels
export AbstractGravityModel

############################################################################################
#                                          Types                                           #
############################################################################################

include("./types.jl")

############################################################################################
#                                        Constants                                         #
############################################################################################

const _D = string(Crayon(; reset = true))
const _B = string(crayon"bold")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./icgem/api.jl")
include("./icgem/compute.jl")
include("./icgem/fetch.jl")
include("./icgem/parse.jl")
include("./icgem/show.jl")

end # module SatelliteToolboxGravityModels
