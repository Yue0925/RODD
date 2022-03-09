
using Random
using PyPlot

Random.seed!(1234) 


T = 12
M = 4
E = [3 for _ in 1:T]
d = [rand(20:70) for _ in 1:T]
f = [10, 30, 60, 90]
e = [8, 6, 4, 2]
h = [1 for _ in 1:T]
p = zeros(Int, T, M)



"""
data : M, T, p, f, h, d, r, e, E
"""
function cplexSolve(r::Int64)
    model = Model(CPLEX.Optimizer) 

    @variable(model, x[1:T, 1:M] >= 0)
    @variable(model, y[1:T, 1:M], Bin)
    @variable(model, s[1:T] >= 0)

    @objective(model, Min,
        sum(p[t, m] * x[t, m] + f[m] * y[t, m] for t in 1:T, m in 1:M ) + sum(h .* s)
    )

    @constraint(model, [t in 2:T],
         sum(x[t, :])  - s[t] + s[t-1] == d[t]
    )

    @constraint(model, sum(x[1, :])  - s[1] == d[1])


    @constraint(model, [t in 1:T, m in 1:M],
        x[t,m] <= sum(d[t_] for t_ in t:T) * y[t,m]
    )

    @constraint(model, [t in 1:T-r],
        sum( (e[m] - E[t_]) * x[t_, m] for t_ in t:t+r, m in 1:M) <= 0
    )


    # solve the problem
    optimize!(model)

    exploredNodes = MOI.get(backend(model), MOI.NodeCount())
    solveTime = MOI.get(model, MOI.SolveTime())

    # status of model
    status = termination_status(model)
    isOptimal = status==MOI.OPTIMAL

    println("isOptimal ? ", isOptimal)
    println("solveTime : ", solveTime)

    if isOptimal
        total_cost = objective_value(model)
        println("total_cost = ", total_cost)
        @show JuMP.value.(x)
        x_star = JuMP.value.(x)

        mean_emission_carbon = sum(e[m] * x_star[t,m] for t in 1:T, m in 1:M) / sum(x_star)
        println(" mean_emission_carbon = ", mean_emission_carbon)

    end

end


function test()
    r = 3
    cplexSolve(r)
end


function analysis_r()
    # for r in 1:
        
    # end
end