module SatelliteToolboxGravityModels

using Crayons
using Dates
using Downloads
using Scratch

import Base: @boundscheck, @propagate_inbounds
import Base: show, throw_boundserror

############################################################################################
#                                        Submodules                                        #
############################################################################################

include("./GravityModels/GravityModels.jl")
using .GravityModels
export GravityModels
export AbstractGravityModel

############################################################################################
#                                          Types                                           #
############################################################################################

include("./types.jl")

############################################################################################
#                                        Constants                                         #
############################################################################################

const _D = string(Crayon(reset = true))
const _B = string(crayon"bold")

const _DT_J2000 = DateTime(2000, 1, 1, 12, 0, 0)

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./icgem/api.jl")
include("./icgem/compute.jl")
include("./icgem/fetch.jl")
include("./icgem/parse.jl")
include("./icgem/show.jl")

end # module SatelliteToolboxGravityModels
