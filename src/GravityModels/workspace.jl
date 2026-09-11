## Description #############################################################################
#
# Workspace with the buffers and precomputed coefficients used to evaluate the gravity
# models.
#
############################################################################################

"""
    struct Workspace{N, T <: AbstractFloat}

Store the buffers and the precomputed coefficients used to evaluate a gravity model,
avoiding allocations and the evaluation of square roots at every call. It must be created
with [`Workspace(model::AbstractGravityModel; kwargs...)`](@ref) and passed to the
evaluation functions using the keyword `workspace`.

The workspace supports computations up to the degree `max_degree` and order `max_order`
with element type `T`, which must be the type obtained by promoting the type of the model
coefficients, the element type of the position, and the type of the time. `N` is the
`Symbol` with the normalization of the model coefficients (`:full`, `:schmidt`, or
`:unnormalized`).

!!! warning

    The workspace holds mutable buffers. Hence, it must not be shared among threads that
    evaluate the model concurrently. Create one workspace per thread instead.

# Fields

- `coefficients::LegendreCoefficients{N, T}`: Precomputed recursion coefficients of the
    associated Legendre functions.
- `P::LowerTriangularStorage{RowMajor, T}`: Buffer to store the associated Legendre
    functions. It is overwritten at every evaluation and must not be modified directly.
- `dP::LowerTriangularStorage{RowMajor, T}`: Buffer to store the derivatives of the
    associated Legendre functions. It is overwritten at every evaluation and must not be
    modified directly.
- `max_degree::Int`: Maximum degree supported by the workspace.
- `max_order::Int`: Maximum order supported by the workspace.
"""
struct Workspace{N, T <: AbstractFloat}
    coefficients::LegendreCoefficients{N, T}
    P::LowerTriangularStorage{RowMajor, T}
    dP::LowerTriangularStorage{RowMajor, T}
    max_degree::Int
    max_order::Int
end

"""
    Workspace(model::AbstractGravityModel{Tm}; kwargs...) -> Workspace

Create a [`Workspace`](@ref) to evaluate the gravity `model` up to the degree `max_degree`
and order `max_order` using the element type `T`.

# Keywords

- `max_degree::Int`: Maximum degree supported by the workspace. If it is higher than the
    maximum degree of `model`, it will be clamped. If it is lower than 0, it will be set to
    the maximum degree of `model`.
    (**Default**: -1)
- `max_order::Int`: Maximum order supported by the workspace. If it is higher than
    `max_degree`, it will be clamped. If it is lower than 0, it will be set to the same
    value as `max_degree`.
    (**Default**: -1)
- `T::Type{<:AbstractFloat}`: Element type of the workspace, which must be the type
    obtained by promoting the type of the model coefficients, the element type of the
    position, and the type of the time used in the evaluations.
    (**Default**: `float(Tm)`)
"""
function Workspace(
    model::AbstractGravityModel{Tm};
    max_degree::Int = -1,
    max_order::Int = -1,
    T::Type{<:AbstractFloat} = float(Tm),
) where {Tm <: Number}
    n_max, m_max = _process_degree_and_order(model, max_degree, max_order)

    coefficients = LegendreCoefficients(coefficient_norm(model), n_max, m_max; T = T)
    P = LowerTriangularStorage{RowMajor, T}(n_max + 1)
    dP = LowerTriangularStorage{RowMajor, T}(n_max + 1)

    return Workspace(coefficients, P, dP, n_max, m_max)
end

function Base.show(io::IO, w::Workspace{N, T}) where {N, T}
    print(io, "Workspace{", repr(N), ", ", T, "}(", w.max_degree, ", ", w.max_order, ")")
    return nothing
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _process_degree_and_order(model::AbstractGravityModel, max_degree::Int, max_order::Int) -> Int, Int

Return the maximum degree and order used in the spherical harmonics computation of `model`
given the requested `max_degree` and `max_order`. If `max_degree` is negative or higher
than the maximum degree of `model`, it is clamped to the latter. If `max_order` is negative
or higher than the selected degree, it is set to the selected degree.

# Returns

- `Int`: Maximum degree used in the computation.
- `Int`: Maximum order used in the computation.
"""
function _process_degree_and_order(
    model::AbstractGravityModel, max_degree::Int, max_order::Int
)
    model_max_degree = maximum_degree(model)

    n_max =
        ((max_degree < 0) || (max_degree > model_max_degree)) ? model_max_degree :
        max_degree
    m_max = ((max_order < 0) || (max_order > n_max)) ? n_max : max_order

    return n_max, m_max
end

"""
    _check_workspace(workspace::Workspace, ::Type{RT}, n_max::Int, m_max::Int) -> Nothing

Check if `workspace` has the element type `RT` and supports the degree `n_max` and the
order `m_max`, throwing an `ArgumentError` otherwise.
"""
function _check_workspace(
    workspace::Workspace{N, T}, ::Type{RT}, n_max::Int, m_max::Int
) where {N, T, RT}
    (T !== RT) && throw(
        ArgumentError(
            "The workspace has element type $T but the computation requires $RT. Create the workspace with `T = $RT`.",
        ),
    )

    ((n_max > workspace.max_degree) || (m_max > workspace.max_order)) && throw(
        ArgumentError(
            "The workspace supports the maximum degree $(workspace.max_degree) and order $(workspace.max_order), but the computation requires degree $n_max and order $m_max.",
        ),
    )

    return nothing
end

"""
    _prepare_potential_inputs(model::AbstractGravityModel, ::Type{RT}, max_degree::Int, max_order::Int, workspace::Union{Nothing, Workspace}) -> Int, Int, Union{Val, LegendreCoefficients}, AbstractMatrix

Process the inputs of the gravitational potential computation of `model` with element type
`RT`, returning the degree and order used in the computation, the object that selects how
the associated Legendre functions are computed, and the matrix to store them.

The requested `max_degree` and `max_order` are clamped as described in
[`_process_degree_and_order`](@ref). If `workspace` is `nothing`, the matrix is allocated
and the Legendre functions are computed using the normalization of `model`. Otherwise, the
buffer and the precomputed coefficients of `workspace` are used, and the function throws an
`ArgumentError` if the workspace is not compatible with the computation.

# Returns

- `Int`: Maximum degree `n_max` used in the computation.
- `Int`: Maximum order `m_max` used in the computation.
- `Union{Val, LegendreCoefficients}`: Normalization or precomputed coefficients used to
    compute the associated Legendre functions.
- `AbstractMatrix`: Matrix `P` to store the associated Legendre functions.
"""
function _prepare_potential_inputs(
    model::AbstractGravityModel,
    ::Type{RT},
    max_degree::Int,
    max_order::Int,
    workspace::Union{Nothing, Workspace},
) where {RT}
    n_max, m_max = _process_degree_and_order(model, max_degree, max_order)

    if isnothing(workspace)
        legendre = coefficient_norm(model)
        P        = LowerTriangularStorage{RowMajor, RT}(n_max + 1)
    else
        _check_workspace(workspace, RT, n_max, m_max)
        legendre = workspace.coefficients
        P        = workspace.P
    end

    return n_max, m_max, legendre, P
end

"""
    _prepare_field_derivative_inputs(model::AbstractGravityModel, ::Type{RT}, max_degree::Int, max_order::Int, workspace::Union{Nothing, Workspace}) -> Int, Int, Int, Int, Union{Val, LegendreCoefficients}, AbstractMatrix, AbstractMatrix

Process the inputs of the gravitational field derivative computation of `model` with
element type `RT`, returning the degrees and orders used in the computation, the object
that selects how the associated Legendre functions are computed, and the matrices to store
them and their derivatives.

The requested `max_degree` and `max_order` are clamped as described in
[`_process_degree_and_order`](@ref). If `workspace` is `nothing`, the matrices are
allocated and the Legendre functions are computed using the normalization of `model`.
Otherwise, the buffers and the precomputed coefficients of `workspace` are used, and the
function throws an `ArgumentError` if the workspace is not compatible with the computation.

# Returns

- `Int`: Maximum degree `n_max` used in the computation.
- `Int`: Maximum order `m_max` used in the computation.
- `Int`: Maximum degree computed in `P`.
- `Int`: Maximum order computed in `P`, which is `m_max + 1` if `m_max < n_max` because
    the derivative computation requires one additional order.
- `Union{Val, LegendreCoefficients}`: Normalization or precomputed coefficients used to
    compute the associated Legendre functions.
- `AbstractMatrix`: Matrix `P` to store the associated Legendre functions.
- `AbstractMatrix`: Matrix `dP` to store the derivatives of the associated Legendre
    functions.
"""
function _prepare_field_derivative_inputs(
    model::AbstractGravityModel,
    ::Type{RT},
    max_degree::Int,
    max_order::Int,
    workspace::Union{Nothing, Workspace},
) where {RT}
    n_max, m_max = _process_degree_and_order(model, max_degree, max_order)

    # To compute the derivative when `m_max < n_max`, the matrix `P` must contain one order
    # more than `dP`. Otherwise, the algorithm would access regions with undefined numbers.
    n_max_P = n_max
    m_max_P = (n_max == m_max) ? m_max : m_max + 1

    if isnothing(workspace)
        legendre = coefficient_norm(model)
        P        = LowerTriangularStorage{RowMajor, RT}(n_max + 1)
        dP       = LowerTriangularStorage{RowMajor, RT}(n_max + 1)
    else
        _check_workspace(workspace, RT, n_max, m_max)
        legendre = workspace.coefficients
        P        = workspace.P
        dP       = workspace.dP
    end

    return n_max, m_max, n_max_P, m_max_P, legendre, P, dP
end

"""
    _legendre!(N::Val, P::AbstractMatrix, θ::Number, n_max::Int, m_max::Int) -> Nothing
    _legendre!(coefficients::LegendreCoefficients, P::AbstractMatrix, θ::Number, n_max::Int, m_max::Int) -> Nothing

Compute in `P` the associated Legendre functions `P_n,m[cos(θ)]` up to the degree `n_max`
and order `m_max` without the Condon-Shortley phase term, using either the normalization
`N` or the precomputed `coefficients`.
"""
function _legendre!(N::Val, P::AbstractMatrix, θ::Number, n_max::Int, m_max::Int)
    legendre!(N, P, θ, n_max, m_max; ph_term = false)
    return nothing
end

function _legendre!(
    coefficients::LegendreCoefficients, P::AbstractMatrix, θ::Number, n_max::Int, m_max::Int
)
    legendre!(P, θ, coefficients, n_max, m_max; ph_term = false)
    return nothing
end

"""
    _dlegendre!(N::Val, dP::AbstractMatrix, θ::Number, P::AbstractMatrix, n_max::Int, m_max::Int) -> Nothing
    _dlegendre!(coefficients::LegendreCoefficients, dP::AbstractMatrix, θ::Number, P::AbstractMatrix, n_max::Int, m_max::Int) -> Nothing

Compute in `dP` the first-order derivatives of the associated Legendre functions
`P_n,m[cos(θ)]` with respect to `θ` [rad] up to the degree `n_max` and order `m_max`
without the Condon-Shortley phase term, using either the normalization `N` or the
precomputed `coefficients`. The matrix `P` must contain the associated Legendre functions
computed with [`_legendre!`](@ref).
"""
function _dlegendre!(
    N::Val, dP::AbstractMatrix, θ::Number, P::AbstractMatrix, n_max::Int, m_max::Int
)
    dlegendre!(N, dP, θ, P, n_max, m_max; ph_term = false)
    return nothing
end

function _dlegendre!(
    coefficients::LegendreCoefficients,
    dP::AbstractMatrix,
    θ::Number,
    P::AbstractMatrix,
    n_max::Int,
    m_max::Int,
)
    dlegendre!(dP, θ, P, coefficients, n_max, m_max; ph_term = false)
    return nothing
end
