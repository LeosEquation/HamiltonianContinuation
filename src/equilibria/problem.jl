struct EquilibriumProblem{
    U<:Real,H,DX,BC,STOP,
    F,J,LPF,LPJ,BPF,BPJ
}
    dimension::Int
    t::U
    init::Vector{U}
    params::Any
    space::Vector{JetSpace}
    H::H
    f!::DX
    bc!::BC
    stop::STOP
    System!::F
    Jacobian!::J
    LP_System!::LPF
    LP_Jacobian!::LPJ
    BP_System!::BPF
    BP_Jacobian!::BPJ
end

function equilibrium_problem(
    init::Vector{U},
    t::U,
    params,
    H,
    f!,
    bc!,
    stop
) where {U<:Real}

    dimension = length(init)
    ordtup = build_ordtup(dimension)

    space = [JetSpace(order=i, numvars=dimension) for i in 1:2]

    System! = build_eq_system(dimension)
    Jacobian! = build_eq_jacobian(dimension)

    LP_System! = build_lp_system(dimension)
    LP_Jacobian! = build_lp_jacobian(dimension, ordtup)

    BP_System! = build_bp_system(dimension)
    BP_Jacobian! = build_bp_jacobian(dimension, ordtup)

    return EquilibriumProblem(
        dimension,
        t,
        init,
        params,
        space,
        H,
        f!,
        bc!,
        stop,
        System!,
        Jacobian!,
        LP_System!,
        LP_Jacobian!,
        BP_System!,
        BP_Jacobian!
    )
end