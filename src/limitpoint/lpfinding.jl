function LPFinding!(
        f!,
        bc!,
        g,
        x::Array{U,1},
        params,
        v::Array{U,1},
        q::Array{U,1}, J::Array{U,2}, F::Array{U,1}, Δq::Array{U,1},
        dxeval::Array{U,1}, Jeval::Array{U,2},
        dy::Array{TaylorN{U},1}, yaux::Array{TaylorN{U},1},
        ordtup::Array{NTuple{N,T},2}, n::T,
        ite::T, atol::U, rtol::U) where {U<:Real,T<:Integer,N}

        q[1:n] .= x
        q[n+1:2*n-1] .= v

        for i in 1:n
                yaux[i][0][1] = q[i]
        end

        bc!(yaux, params, 0.0)

        f!(dy, yaux, params, zero(U))

        for i in 1:n-1
                for j in 1:n
                        Jeval[i, j] = dy[i][1][j]
                end
                dxeval[i] = dy[i][0][1]
        end

        LPSystem!(F, Jeval, dxeval, q, n)
        LPJacobian!(J, Jeval, dy, q, ordtup, n)

        if abs(det(J)) < (atol + rtol * norm(q, Inf))
                @warn("System Jacobian was been singular in LP bifurcation. Stopping continuation.")
                for qq in Δq
                        qq = NaN
                end
                return nothing
        end

        Δq .= J \ F

        k = 1

        while k <= ite && norm(Δq) > (atol + rtol * norm(norm(q, Inf), Inf))

                q .-= Δq

                if g(q, params, 0.0)
                        q .+= J \ F
                        break
                end

                for i in 1:n
                        yaux[i][0][1] = q[i]
                end

                bc!(yaux, params, 0.0)
                f!(dy, yaux, params, zero(U))

                for i in 1:n-1
                        for j in 1:n
                                Jeval[i, j] = dy[i][1][j]
                        end
                        dxeval[i] = dy[i][0][1]
                end

                LPSystem!(F, Jeval, dxeval, q, n)
                LPJacobian!(J, Jeval, dy, q, ordtup, n)

                if abs(det(J)) < (atol + rtol * norm(q, Inf))
                        @warn("System Jacobian was been singular in LP bifurcation. Stopping continuation.")
                        return nothing
                end


                Δq .= J \ F

                k += 1

        end

        for i in 1:n
                x[i] = q[i]
        end

        for i in 1:n-1
                v[i] = q[n+i]
        end

        return nothing
end



# function LPFinding!(f!, bc!, g, x::Array{U, 1}, params, v::Array{U, 1},
#                     q::Array{U, 1}, J::Array{U, 2}, F::Array{U, 1}, Δq::Array{U, 1},
#                     dxeval::Array{U, 1}, Jeval::Array{U, 2},
#                     dy::Array{TaylorN{U}, 1}, yaux::Array{TaylorN{U},1}, 
#                     ordtup::Array{NTuple{N, T}, 2}, n::T,
#                     ite::T, atol::U, rtol::U) where {U<:Real, T<:Integer, N}

#         q[1:n] .= x
#         q[n+1:2*n-1] .= v

#         for i in 1:n
#                 yaux[i][0][1] = q[i]
#         end

#         f!(dy, yaux, params, zero(U))

#         for i in 1:n-1
#                 for j in 1:n
#                         Jeval[i,j] = dy[i][1][j]
#                 end
#                 dxeval[i] = dy[i][0][1]
#         end

#         LPSystem!(F, Jeval, dxeval, q, n)
#         LPJacobian!(J, Jeval, dy, q, ordtup, n)

#         if abs(det(J)) < ( atol + rtol * norm(q1, Inf) )
#             @warn("System Jacobian was been singular in LP bifurcation. Stopping continuation.")
#             for qq in Δq
#                 qq = NaN
#                 end
#             return nothing
#         end

#         Δq .= J \ F 

#         k = 1

#         while k <= ite && norm(Δq) > ( atol + rtol * norm(norm(q1, Inf), Inf) )

#                 q .-= Δq

#                 bc!(q, params, 0.0)

#                 if g(q, params, 0.0)
#                         q .+= J \ F 
#                         break
#                 end

#                 for i in 1:n
#                         yaux[i][0][1] = q[i]
#                 end

#                 f!(dy, yaux, params, zero(U))

#                 for i in 1:n-1
#                         for j in 1:n
#                                 Jeval[i,j] = dy[i][1][j]
#                         end
#                         dxeval[i] = dy[i][0][1]
#                 end

#                 LPSystem!(F, Jeval, dxeval, q, n)
#                 LPJacobian!(J, Jeval, dy, q, ordtup, n)

#                 if abs(det(J)) < ( atol + rtol * norm(q1, Inf) )
#                         @warn("System Jacobian was been singular in LP bifurcation. Stopping continuation.")
#                         return nothing
#                 end


#                 Δq .= J \ F

#                 k += 1

#         end

#         for i in 1:n
#                 x[i] = q[i]
#         end

#         for i in 1:n-1
#                 v[i] = q[n+i]
#         end

#         return nothing

# end