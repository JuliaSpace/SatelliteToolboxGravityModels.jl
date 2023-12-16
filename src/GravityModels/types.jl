# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Types and structures for gravity models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export AbstractGravityModel

"""
    abstract type AbstractGravityModel{T<:Number}

Abstract data type of all gravity models.
"""
abstract type AbstractGravityModel{T<:Number} end
