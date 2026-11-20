import TaylorSeries: get_numvars, ntuple

"""
    build_ordtup(problem::EquilibriumProblem)

Construct a matrix of multi-index tuples used to access mixed derivatives
directly from the Taylor polynomial coefficients.

For an equilibrium problem with `n` variables, the returned matrix has
dimensions `n × n`. Each entry `ordtup[i, j]` is a tuple whose components
represent the exponents associated with the mixed derivative with respect
to variables `i` and `j`.

These multi-indices allow mixed derivatives to be computed directly from
the Taylor expansion coefficients, without constructing additional
Taylor polynomials through successive differentiation.

Although this approach requires explicit bookkeeping of the multi-indices,
it significantly reduces memory allocations by avoiding the creation of
intermediate Taylor polynomial objects.

# Arguments
- `problem::EquilibriumProblem`: Equilibrium problem defining the number
  of independent variables.

# Returns
- `Matrix{NTuple}`: Matrix of multi-indices associated with the second-order
  mixed derivatives.

# Example
For a system with three variables, the matrix contains the multi-indices
corresponding to the derivatives with respect to each pair of variables.
For instance, the entry associated with variables `x` and `y` represents
the multi-index of the coefficient used to evaluate

    ∂²F / ∂x∂y

directly from the Taylor expansion.
"""
function build_ordtup(dimension::Int)

    return [ntuple(k -> count(==(k), (i, j)), get_numvars()) for i in 1:dimension, j in 1:dimension]

end