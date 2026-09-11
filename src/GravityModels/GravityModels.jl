## Description #############################################################################
#
# Submodule to define the gravity model API.
#
############################################################################################

module GravityModels

using Dates
using SatelliteToolboxBase
using SatelliteToolboxLegendre
using StaticArrays

import SatelliteToolboxBase: LowerTriangularStorage, RowMajor

############################################################################################
#                                          Types                                           #
############################################################################################

include("./types.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./api.jl")
include("./accelerations.jl")
include("./gravitational_field_derivative.jl")
include("./potential.jl")
include("./time.jl")

end # module GravityModels
