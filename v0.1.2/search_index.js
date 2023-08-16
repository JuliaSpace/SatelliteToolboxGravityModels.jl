var documenterSearchIndex = {"docs":
[{"location":"lib/library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"lib/library/","page":"Library","title":"Library","text":"Documentation for SatelliteToolboxGravityModels.jl.","category":"page"},{"location":"lib/library/","page":"Library","title":"Library","text":"Modules = [SatelliteToolboxGravityModels, GravityModels]","category":"page"},{"location":"lib/library/#SatelliteToolboxGravityModels.fetch_icgem_file-Tuple{Symbol}","page":"Library","title":"SatelliteToolboxGravityModels.fetch_icgem_file","text":"fetch_icgem_file(url::AbstractString; kwargs...) -> String\nfetch_icgem_file(model::Symbol; kwargs...) -> String\n\nFetch a ICGEM file from the url and return its file path to be parsed with the function GravityModels.load. If the file already exists, it will not be re-downloaded unless the keyword force = true is passed.\n\nA symbol can be passed instead the URL to fetch pre-configured gravity field models. The supported values are:\n\n:EGM96: Earth Gravitational Model from 1996.\n:EGM2008: Earth Gravitational Model from 2008.\n:JGM2: Joint Gravity Model 2.\n:JGM3: Joint Gravity Model 3.\n\nExamples\n\njulia> fetch_icgem_file(:EGM96)\n[ Info: Downloading the ICGEM file 'EGM96.gfc' from 'http://icgem.gfz-potsdam.de/getmodel/gfc/971b0a3b49a497910aad23cd85e066d4cd9af0aeafe7ce6301a696bed8570be3/EGM96.gfc'...\n\"/Users/ronan.arraes/.julia/scratchspaces/bd9e9728-6f7b-4d28-9e50-c765cb1b7c8c/icgem/EGM96.gfc\"\n\njulia> fetch_icgem_file(:EGM96)\n\"/Users/ronan.arraes/.julia/scratchspaces/bd9e9728-6f7b-4d28-9e50-c765cb1b7c8c/icgem/EGM96.gfc\"\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#SatelliteToolboxGravityModels.icgem_coefficients-Union{Tuple{T}, Tuple{IcgemFile{T}, Int64, Int64, Dates.DateTime}} where T<:Number","page":"Library","title":"SatelliteToolboxGravityModels.icgem_coefficients","text":"icgem_coefficients(model::IcgemFile{T}, degree::Int, order::Int, t::DateTime) where T<:Number -> T, T\n\nCompute the ICGEM coefficients (Clm and Slm) of the model for the specified degree and order in the instant t.\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#SatelliteToolboxGravityModels.parse_icgem","page":"Library","title":"SatelliteToolboxGravityModels.parse_icgem","text":"parse_icgem(filename::AbstractString, T::DataType = Float64) -> IcgemFile\n\nParse the ICGEM file filename using the data type T.\n\nnote: Note\nT is converted to float to obtain the output type.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel","text":"abstract type AbstractGravityModel{T<:Number}\n\nAbstract data type of all gravity models.\n\n\n\n\n\n","category":"type"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.coefficient_norm","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.coefficient_norm","text":"coefficient_norm(model::AbstractGravityModel) where T<:Number -> Symbol\n\nReturn the normalization we must use in the spherical harmonics when computing the Legendre associated functions. The accepted values are:\n\n:full: Use full normalization.\n:schmidt: Use Schmidt quasi-normalization.\n:unnormalized: Do not perform normalization.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.coefficients","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.coefficients","text":"coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::DataTime) where T<:Number -> T, T\n\nReturn the Clm and Slm coefficients of the gravity model for the specified degree, order, and time.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.gravitational_acceleration-Union{Tuple{T}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector, Dates.DateTime}} where T<:Number","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.gravitational_acceleration","text":"gravitational_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number -> NTuple{3, T}\n\nCompute the gravitational acceleration [m / s²] represented in ITRF using the model in the position r [m], also represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.\n\nnote: Note\nGravitational acceleration is the acceleration caused by the central body mass only, i.e., without considering the centrifugal potential.\n\nKeywords\n\nmax_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.\n\nReturns\n\nT: The derivative of the gravitational field w.r.t. the radius (∂U/∂r).\nT: The derivative of the gravitational field w.r.t. the geocentric latitude (∂U/∂ϕ).\nT: The derivative of the gravitational field w.r.t. the longitude (∂U/∂λ).\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.gravitational_field_derivative-Union{Tuple{T}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector, Dates.DateTime}} where T<:Number","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.gravitational_field_derivative","text":"gravitational_field_derivative(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number -> NTuple{3, T}\n\nCompute the gravitational field derivative [SI] with respect to the spherical coordinates (∂U/∂r, ∂U/∂ϕ, ∂U/∂λ) using the model in the position r [m], represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.\n\ninfo: Info\nIn this case, ϕ is the geocentric latitude and λ is the longitude.\n\nKeywords\n\nmax_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.\n\nReturns\n\nT: The derivative of the gravitational field w.r.t. the radius (∂U/∂r).\nT: The derivative of the gravitational field w.r.t. the geocentric latitude (∂U/∂ϕ).\nT: The derivative of the gravitational field w.r.t. the longitude (∂U/∂λ).\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.gravity_acceleration-Union{Tuple{T}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector}, Tuple{SatelliteToolboxGravityModels.GravityModels.AbstractGravityModel{T}, AbstractVector, Dates.DateTime}} where T<:Number","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.gravity_acceleration","text":"gravity_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number -> NTuple{3, T}\n\nCompute the gravity acceleration [m / s²] represented in ITRF using the model in the position r [m], also represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.\n\nnote: Note\nGravity acceleration is the compound acceleration caused by the central body mass and the centrifugal force due to the planet's rotation.\n\nKeywords\n\nmax_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.\n\nReturns\n\nT: The derivative of the gravitational field w.r.t. the radius (∂U/∂r).\nT: The derivative of the gravitational field w.r.t. the geocentric latitude (∂U/∂ϕ).\nT: The derivative of the gravitational field w.r.t. the longitude (∂U/∂λ).\n\n\n\n\n\n","category":"method"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.gravity_constant","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.gravity_constant","text":"gravity_constant(model::AbstractGravityModel{T}) where T<:Number -> T\n\nReturn the gravity constant [m³ / s²] for the gravity model.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.load","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.load","text":"load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T\n\nLoad a gravity model of type T using the arguments args... and keywords kwargs....\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.maximum_degree","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.maximum_degree","text":"maximum_degree(model::AbstractGravityModel) -> Int\n\nReturn the maximum degree of the gravity model.\n\n\n\n\n\n","category":"function"},{"location":"lib/library/#SatelliteToolboxGravityModels.GravityModels.radius","page":"Library","title":"SatelliteToolboxGravityModels.GravityModels.radius","text":"radius(model::AbstractGravityModel{T}) where T<:Number -> T\n\nReturn the radius [m] for the gravity model.\n\n\n\n\n\n","category":"function"},{"location":"man/usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"CurrentModule = SatelliteToolboxGravityModels\nDocTestSetup = quote\n    using SatelliteToolboxGravityModels\nend","category":"page"},{"location":"man/usage/#Initialization","page":"Usage","title":"Initialization","text":"","category":"section"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"We can initialize a gravity model using the function:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"where the arguments and keywords depend on the gravity model type T. For ICGEM files, we must use T = IcgemFile and the following signature:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function GravityModels.load(::Type{IcgemFile}, filename::AbstractString, T::DataType = Float64)","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"where it loads the ICGEM file in the path filename converting the coefficients to the type T.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"We also provide a function to help downloading the ICGEM files:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function fetch_icgem_file(url::AbstractString; kwargs...)\nfunction fetch_icgem_file(model::Symbol; kwargs...)","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"It fetches a ICGEM file from the url and return its file path to be parsed with the function GravityModels.load. If the file already exists, it will not be re-downloaded unless the keyword force = true is passed.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"Notice that the functions downloads the files to a scratch space.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"A symbol can be passed instead the URL to fetch pre-configured gravity field models. The supported values are:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":":EGM96: Earth Gravitational Model from 1996.\n:EGM2008: Earth Gravitational Model from 2008.\n:JGM2: Joint Gravity Model 2.\n:JGM3: Joint Gravity Model 3.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"Finally, we can initialize, for example, the EGM96 model using:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"julia> egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96))\nIcgemFile{Float64}:\n      Product type : gravity_field\n       Model name  : EGM96\n  Gravity constant : 3.986004415e14\n            Radius : 6.3781363e6\n    Maximum degree : 360\n            Errors : formal\n       Tide system : tide_free\n              Norm : fully_normalized\n         Data type : Float64","category":"page"},{"location":"man/usage/#Gravitational-Field-Derivative","page":"Usage","title":"Gravitational Field Derivative","text":"","category":"section"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The following function:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function gravitational_field_derivative(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"computes the gravitational field derivative [SI] with respect to the spherical coordinates:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"fracpartial Upartial r fracpartial Upartial phi fracpartial Upartial lambda","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"using the model in the position r [m], represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"info: Info\nIn this case, phi is the geocentric latitude and lambda is the longitude.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The following keywords are available:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"max_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"julia> GravityModels.gravitational_field_derivative(egm96, [6378.137e3, 0, 0])\n(-9.814284376497435, 49.45906319417034, -115.71285105900459)","category":"page"},{"location":"man/usage/#Gravitational-Acceleration","page":"Usage","title":"Gravitational Acceleration","text":"","category":"section"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The gravitational acceleration is the acceleration caused by the central body mass only, i.e., without considering the centrifugal potential. We can compute it using the function:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function gravitational_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"where it returns the gravitational field acceleration [m / s²] represented in ITRF using the model in the position r [m], also represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The following keywords are available:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"max_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"julia> GravityModels.gravitational_acceleration(egm96, [6378.137e3, 0, 0])\n3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):\n -9.814284376497435\n -1.814210812013047e-5\n  7.754468615862334e-6","category":"page"},{"location":"man/usage/#Gravity-Acceleration","page":"Usage","title":"Gravity Acceleration","text":"","category":"section"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The gravity acceleration is the compound acceleration caused by the central body mass and the centrifugal force due to the planet's rotation. We can compute it using the function:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"function gravity_acceleration(model::AbstractGravityModel{T}, r::AbstractVector, time::DateTime = DateTime(\"2000-01-01\"); kwargs...) where T<:Number","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"where it computes the gravity acceleration [m / s²] represented in ITRF using the model in the position r [m], also represented in ITRF, at instant time. If the latter argument is omitted, the J2000.0 epoch is used.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"The following keywords are available:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"max_degree::Int: Maximum degree used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than the available number of   coefficients in the model, it will be clamped. If it is lower than 0, it will be set   to the maximum degree available. (Default = -1)\nmax_order::Int: Maximum order used in the spherical harmonics when computing the   gravitational field derivative. If it is higher than max_degree, it will be clamped.   If it is lower than 0, it will be set to the same value as max_degree.   (Default = -1)\nP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   coefficients, reducing the allocations. If it is nothing, the matrix will be created   when calling the function.\ndP::Union{Nothing, AbstractMatrix}: An optional matrix that must contain at least   max_degree + 1 × max_degree + 1 real numbers that will be used to store the Legendre   derivative coefficients, reducing the allocations. If it is nothing, the matrix will   be created when calling the function.","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"Thus, we can compute the gravity acceleration in Equator using the EGM96 model by:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"julia> egm96 = GravityModels.load(IcgemFile, fetch_icgem_file(:EGM96));\n\njulia> GravityModels.gravitational_acceleration(egm96, [6378.137e3, 0, 0])\n3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):\n -9.814284376497435\n -1.814210812013047e-5\n  7.754468615862334e-6","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"Whereas we can obtain the gravity acceleration at the poles by:","category":"page"},{"location":"man/usage/","page":"Usage","title":"Usage","text":"julia> GravityModels.gravitational_acceleration(egm96, [0, 0, 6356.7523e3])\n3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):\n -6.12152785935481e-5\n  0.0\n -9.83208158872835","category":"page"},{"location":"#SatelliteToolboxGravityModels.jl","page":"Home","title":"SatelliteToolboxGravityModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package implements the support of gravity models for the SatelliteToolbox.jl ecosystem. We can use it to obtain, for example, the gravitational acceleration to implement highly accurate numerical orbit propagators.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, we have the following functionalities:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Compute the gravity field derivative in spherical coordinates;\nCompute the gravitational acceleration; and\nCompute the gravity acceleration, taking into account the Earth's rotation rate.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package contains an API to allow the user to access any gravity model while computing the previous functions. We also provide in-built support for ICGEM files.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.install(\"SatelliteToolboxGravityModels\")","category":"page"},{"location":"man/api/#GravityModels-API","page":"API","title":"GravityModels API","text":"","category":"section"},{"location":"man/api/","page":"API","title":"API","text":"This document describes the API required for a gravity model. The user can add new models by overloading the functions listed here.","category":"page"},{"location":"man/api/#Structure","page":"API","title":"Structure","text":"","category":"section"},{"location":"man/api/","page":"API","title":"API","text":"All models require a structure with supertype AbstractGravityModel{T<:Number}, where T is the type of the coefficients in the model.","category":"page"},{"location":"man/api/#API-Functions","page":"API","title":"API Functions","text":"","category":"section"},{"location":"man/api/","page":"API","title":"API","text":"function coefficients(model::AbstractGravityModel{T}, degree::Int, order::Int, time::DataTime) where T<:Number -> T, T","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return the coefficients Clm and Slm of the gravity model for the specified degree, order, and time. Hence:","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"coefficients(model, 10, 8, DateTime(\"2023-06-19\"))","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"must return a Tuple{T, T} with the Clm and Slm, respectively, for the degree 10, order 8, and computed at day 2023-06-19.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"Note If the model has constant coefficients, the function must still accept the positional argument time, but it will be neglected.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"function coefficient_norm(model::AbstractGravityModel) where T<:Number -> Symbol","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return the normalization we must use in the spherical harmonics when computing the Legendre associated functions. The accepted values are:","category":"page"},{"location":"man/api/","page":"API","title":"API","text":":full: Use full normalization.\n:schmidt: Use Schmidt quasi-normalization.\n:unnormalized: Do not perform normalization.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"function gravity_constant(model::AbstractGravityModel{T}) where T<:Number -> T","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return the gravity constant [m³ / s²] for the gravity model.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"function load(::Type{T}, args...; kwargs...) where T<:AbstractGravityModel -> T","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return gravity model structure, which is loaded using the arguments args... and keywords kwargs....","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"function maximum_degree(model::AbstractGravityModel) -> Int","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return the maximum degree of the gravity model.","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"function radius(model::AbstractGravityModel{T}) where T<:Number -> T","category":"page"},{"location":"man/api/","page":"API","title":"API","text":"This function must return the radius [m] for the gravity model.","category":"page"}]
}
