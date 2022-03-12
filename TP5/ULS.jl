
using Random
import PyPlot; const plt = PyPlot

"""
tape in julia terminal : 

ENV["PYTHON"]=""
Pkg.build("PyCall")
"""


Random.seed!(1234) 

#function data()
global T = 12
global M = 4
global E = [3 for _ in 1:T]
global d = [rand(20:70) for _ in 1:T]
global f = [10, 30, 60, 90]
global e = [8, 6, 4, 2]
global h = [1 for _ in 1:T]
global p = zeros(Int, T, M)
#end


"""
parameter data : M, T, p, f, h, d, r, e, E
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

    set_silent(model) # turn off cplex output

    # solve the problem
    optimize!(model)

    exploredNodes = MOI.get(backend(model), MOI.NodeCount())
    solveTime = round(MOI.get(model, MOI.SolveTime()), digits = 2)

    # status of model
    status = termination_status(model)
    isOptimal = status==MOI.OPTIMAL

    println("isOptimal ? ", isOptimal)
    println("solveTime : ", solveTime)

    if isOptimal
        total_cost = round(objective_value(model), digits = 2)
        println("total_cost = ", total_cost)
        # @show JuMP.value.(x)
        x_star = JuMP.value.(x)

        mean_emission_carbon = round(sum(e[m] * x_star[t,m] for t in 1:T, m in 1:M) / sum(x_star), digits = 2)
        println(" mean_emission_carbon = ", mean_emission_carbon)
        return total_cost, mean_emission_carbon, solveTime
    else
        error("cplex doesn't find the optimal solution !")
    end

end

function test()
    r = 3
    cplexSolve(r)
end


function analysis_r()
    global d
    global T

    for test in 1:4
        seed = rand(1:10000)
        Random.seed!(seed) 
        d = [rand(20:70) for _ in 1:T]
        # @show d

        record_r = []
        record_cost = []
        record_emission = []
        record_times = []
        for r in 1:T
            total_cost, mean_emission_carbon, solveTime = cplexSolve(r)
            @show total_cost, mean_emission_carbon, solveTime
            append!(record_r, r)
            append!(record_cost, total_cost)
            append!(record_emission, mean_emission_carbon)
            append!(record_times, solveTime)
        end

        fig, ax = plt.subplots()
        ax.plot(record_r, record_cost, color="red", marker="o", linewidth=2.0, linestyle="--")
        ax.set_xlabel("L'intervalle r", fontsize=14)
        ax.set_ylabel("Coût total", color="red", fontsize=14)
        ax.set_title("L'évolution du coût total et l'émission carbone moyenne", fontsize=14)

        ax2=ax.twinx()
        ax2.plot(record_r, record_emission, color="blue", marker="o", linewidth=2.0, linestyle="--")
        ax2.set_ylabel("L'émission carbone moyenne", color="blue", fontsize=14)
        # plt.show()

        savefig("analysis_r_test$test" * ".png")
        # plt.close()
    end
end