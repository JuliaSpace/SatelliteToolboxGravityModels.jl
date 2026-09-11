module SatelliteToolboxGravityModels

using Dates
using Downloads
using Scratch

import Base: show
import SatelliteToolboxBase: EARTH_ANGULAR_SPEED, LowerTriangularStorage, RowMajor
import SatelliteToolboxBase:
    PrintedField, PrintedSection, format_value, print_tree, type_name

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
#                                         Includes                                         #
############################################################################################

include("./icgem/api.jl")
include("./icgem/compute.jl")
include("./icgem/fetch.jl")
include("./icgem/parse.jl")
include("./icgem/show.jl")

end # module SatelliteToolboxGravityModels
