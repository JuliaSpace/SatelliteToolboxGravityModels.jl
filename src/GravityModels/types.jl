## Description #############################################################################
#
# Types and structures for gravity models.
#
############################################################################################

export AbstractGravityModel

"""
    abstract type AbstractGravityModel{T<:Number, NT<:Val}

Abstract data type of all gravity models.
"""
abstract type AbstractGravityModel{T<:Number, NT<:Val} end
