struct EquilibriumProblem{U<:Real,F,B,S,P} <: AbstractProblem
    f!::F
    bc!::B
    stop::S
    params::P
    t0::U
    dim_variables::Int
    dim_params::Int
    dim::Int
    x::Vector{U}
    dx::Vector{U}
    x_jet::Vector{TaylorN{U}}
    dx_jet::Vector{TaylorN{U}}
end

function EquilibriumProblem(;
    f!,
    bc!,
    stop,
    x_ini,
    params,
    t0,
    dim_params
)
    dim = length(x_ini)
    dim_variables = dim - dim_params

    dx = zero(x_ini)

    equilibrium_space = JetSpace()

    x_jet
    dx_jet

    return EquilibriumProblem(
        f!,
        bc!,
        stop,
        params,
        t0,
        dim_variables,
        dim_params,
        dim,
        x_ini,
        dx,
        x_jet,
        dx_jet
    )
end