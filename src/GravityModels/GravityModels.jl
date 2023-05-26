# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Submodule to define the gravity model API.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

module GravityModels

using Dates
using ReferenceFrameRotations
using SatelliteToolboxBase
using SatelliteToolboxLegendre
using StaticArrays

############################################################################################
#                                          Types
############################################################################################

include("./types.jl")

############################################################################################
#                                         Includes
############################################################################################

include("./api.jl")
include("./accelerations.jl")
include("./gravitational_field_derivative.jl")

end # module GravityModels
