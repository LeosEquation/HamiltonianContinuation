import TaylorSeries: JetSpace, variables, constant_term!, order

"""
    TaylorWorkspace{U,N}

Stores the Taylor expansions and auxiliary multi-indices required during
a continuation process.

# Fields
- `x_jet_1::Vector{TaylorN{U}}`: First-order Taylor expansions of the
  position variables.
- `dx_jet_1::Vector{TaylorN{U}}`: First-order Taylor expansions of the
  corresponding vector field.
- `x_jet_2::Vector{TaylorN{U}}`: Second-order Taylor expansions of the
  position variables.
- `dx_jet_2::Vector{TaylorN{U}}`: Second-order Taylor expansions of the
  corresponding vector field.
- `ordtup::Matrix{NTuple{N,Int}}`: Matrix of multi-indices associated
  with all second-order partial derivatives. Each entry specifies the
  differentiation order with respect to each independent variable.

The concrete Taylor expansions depend on the continuation method and
the problem definition.
"""
struct TaylorWorkspace{U<:Real,N}
    x_jet_1::Vector{TaylorN{U}}
    dx_jet_1::Vector{TaylorN{U}}
    x_jet_2::Vector{TaylorN{U}}
    dx_jet_2::Vector{TaylorN{U}}
    ordtup::Matrix{NTuple{N,Int}}
end

"""
    build_taylor_workspace(::Type{C}, problem, space)

Builds the Taylor expansions and auxiliary multi-indices required by
the continuation method `C`.

Constructs first- and second-order Taylor variables using the provided
Taylor space, sets their constant terms to the initial conditions, and
evaluates the corresponding vector fields. It also builds the matrix
of multi-indices used to compute second-order partial derivatives.

The provided Taylor space must have the same number of independent
variables as the initial condition vector and a maximum expansion
order of at least two.

# Arguments
- `problem::AbstractProblem`: Problem defining the system to be continued.
- `space::JetSpace`: Taylor space used to construct the multivariate
  Taylor expansions.

# Returns
A `TaylorWorkspace` containing the first- and second-order Taylor
expansions of the position variables and vector field, together with
the corresponding second-order derivative multi-indices.

# Throws
- `DimensionMismatch`: If the number of independent variables in
  `space` does not match the length of the initial condition vector.
- `ArgumentError`: If the maximum expansion order of `space` is less
  than two.
"""
function build_taylor_workspace(
    problem::P,
    space::JetSpace
) where {P<:AbstractProblem}

    (; f!, x_ini, params, t0) = problem

    U = eltype(x_ini)
    n = length(x_ini)

    ordtup = [
        ntuple(k -> count(==(k), (i, j)), n)
        for i in 1:n, j in 1:n
    ]

    # Validate the Taylor space
    space.num_vars == n ||
        throw(DimensionMismatch(
            "The number of variables in the JetSpace must match the length of x_ini."
        ))

    order(space) >= 2 ||
        throw(ArgumentError(
            "The JetSpace must have a maximum expansion order of at least 2."
        ))

    # Build first- and second-order Taylor variables
    x_jet_1 = variables(U, space; order=1)
    x_jet_2 = variables(U, space; order=2)

    # Set the initial conditions
    for i in eachindex(x_ini)
        constant_term!(x_jet_1[i], x_ini[i])
        constant_term!(x_jet_2[i], x_ini[i])
    end

    # Allocate the corresponding vector fields
    dx_jet_1 = zero.(x_jet_1)
    dx_jet_2 = zero.(x_jet_2)

    # Evaluate the vector field
    f!(dx_jet_1, x_jet_1, params, t0)
    f!(dx_jet_2, x_jet_2, params, t0)

    return TaylorWorkspace(
        x_jet_1,
        dx_jet_1,
        x_jet_2,
        dx_jet_2,
        ordtup
    )
end