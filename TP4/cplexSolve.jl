using Random
using PyPlot
import PyPlot; const plt = PyPlot

global solveTime, e1_size, e2_size, total_celles_not_cut, total_boundary_edges, obj_val
global m,n,t
global P=1
global w1=1
global w2=5
global L=3
global g=1.26157

"""
Linearization of binary quadratic problem
"""
function cplexSolveQuadratricRelaxed()
    global solveTime, e1_size, e2_size, total_celles_not_cut, total_boundary_edges, obj_val

    M = Model(CPLEX.Optimizer)

    @variable(M, x[1:m+2, 1:n+2], Bin)
    @variable(M, y[1:m+2, 1:n+2, 1:m+2, 1:n+2] >=0)

    A(i, j) = [(i+1, j), (i, j+1)]

    @objective(M, Max, w1 * sum(t[i-1][j-1] * (1 - x[i, j]) for i in 2:m+1, j in 2:n+1) + w2 * g * L *
        sum( sum( x[i, j] - y[i, j, i_, j_] + x[i_, j_] - y[i, j, i_, j_] for (i_, j_) in A(i, j)) for i in 1:m+1, j in 1:n+1) )


    for i in 1:m+1, j in 1:n+1
        for (i_, j_) in A(i, j)
            @constraint(M, -y[i, j, i_, j_] + x[i, j] + x[i_, j_] <= 1)
        end
    end


    @constraint(M, [i in 1:m+2], x[i, 1] == 0)
    @constraint(M, [i in 1:m+2], x[i, n+2] == 0)
    @constraint(M, [j in 1:n+2], x[1, j] == 0)
    @constraint(M, [j in 1:n+2], x[m+2, j] == 0)


    @constraint(M, [i in 1:m+1], -y[i, n+2, i+1, n+2] + x[i, n+2] + x[i+1, n+2] <= 1)

    @constraint(M, [j in 1:n+1], -y[m+2, j, m+2, j+1] + x[m+2, j] + x[m+2, j+1] <= 1)


    #TODO : imposed uncut number ≥ 60
    # @constraint(M, sum(x) ≥ 60)


    # solve the problem
    set_silent(M) # turn off cplex output
    optimize!(M)

    
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL


    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = round(MOI.get(M, MOI.SolveTime()), digits = 3)

    println("exploredNodes ", exploredNodes)
    println("solveTime ", solveTime)

    # display solution
    println("isOptimal ? ", isOptimal)
    println("GAP = ", GAP)


    if has_values(M) && isOptimal

        e1_size = round(Int, sum(t[i-1][j-1] * (1 - JuMP.value(x[i, j])) for i in 2:m+1, j in 2:n+1))
        println("Species size e1 = ", e1_size)


        e2_size = round( g * L * sum( sum( JuMP.value(x[i, j]) -
            JuMP.value(y[i, j, i_, j_]) + JuMP.value(x[i_, j_]) - JuMP.value(y[i, j, i_, j_]) for (i_, j_) in A(i, j)
            ) for i in 1:m+1, j in 1:n+1
        ), digits = 3)
        println("Species size e2 = ", e2_size)


        # println(JuMP.value.(y[1, :, :, :]))
        # println(value.(x)) 

        total_celles_not_cut = round(Int, sum(JuMP.value.(x)))
        println("the number of celles not cut : ", total_celles_not_cut)

        total_boundary_edges = round(Int, sum( sum( JuMP.value(x[i, j]) -
                JuMP.value(y[i, j, i_, j_]) + JuMP.value(x[i_, j_]) - JuMP.value(y[i, j, i_, j_]) for (i_, j_) in A(i, j)
            ) for i in 1:m+1, j in 1:n+1
        ))
        println("total boundary edges : ", total_boundary_edges)


        obj_val = round(objective_value(M), digits = 4)
        println("obj_val : ", obj_val)

    else
        error("cplex cannot find opt sol ! ")
    end

    
end



"""
MIP model
"""
function cplexSolveP1()
    global solveTime, e1_size, e2_size, total_celles_not_cut, total_boundary_edges, obj_val

    M = Model(CPLEX.Optimizer)

    @variable(M, x[1:m+2, 1:n+2], Bin)
    @variable(M, d[1:m+2, 1:n+2] >= 0, Int)

    A(i, j) = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]

    @objective(M, Max, w1 * sum(t[i-1][j-1] * (1 - x[i, j]) for i in 2:m+1, j in 2:n+1)
        + w2 * g * L * sum( 4 * x[i, j] - d[i, j] for i in 2:m+1, j in 2:n+1) )


    @constraint(M, [i in 2:m+1, j in 2:n+1], d[i, j] >= sum(x[k, l] for (k, l) in A(i, j)) - size(A(i, j), 1) * (1 - x[i, j]))


    @constraint(M, [i in 1:m], x[i, 1] == 0)
    @constraint(M, [i in 1:m], x[i, n+2] == 0)
    @constraint(M, [j in 1:n], x[1, j] == 0)
    @constraint(M, [j in 1:n], x[m+2, j] == 0)

    #TODO : imposed uncut number ≥ 60
    # @constraint(M, sum(x) ≥ 60)

    # solve the problem
    set_silent(M) # turn off cplex output
    optimize!(M)
    # println(solution_summary(M))
    
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL


    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = round(MOI.get(M, MOI.SolveTime()), digits = 3)

    println("exploredNodes ", exploredNodes)
    println("solveTime ", solveTime)

    # display solution
    println("isOptimal ? ", isOptimal)
    println("GAP = ", GAP)


    if has_values(M) && isOptimal

        e1_size = round(Int, sum(t[i-1][j-1] * (1 - JuMP.value(x[i, j])) for i in 2:m+1, j in 2:n+1))
        println("Species size e1 = ", e1_size)


        e2_size = round( g * L * sum( 4 * JuMP.value(x[i, j]) - JuMP.value(d[i, j]) for i in 2:m+1, j in 2:n+1), digits = 3)
        println("Species size e2 = ", e2_size)


        total_celles_not_cut = round(Int, sum(JuMP.value.(x)))
        println("the number of celles not cut : ", total_celles_not_cut)

        total_boundary_edges = round(Int, sum( 4 * JuMP.value(x[i, j]) - JuMP.value(d[i, j]) for i in 2:m+1, j in 2:n+1))
        println("total boundary edges : ", total_boundary_edges)


        obj_val = round(objective_value(M), digits = 4)
        println("obj_val : ", obj_val)

    else
        error("cplex cannot find opt sol ! ")
    end

    
end

function run()
    include("ExplForet_opl.dat")
    # include("ExplForet2_opl.dat")


    println("\n\n-----------------")
    println("------ P1 -------")
    println("-----------------")
    cplexSolveP1()


    println("\n\n-----------------")
    println("------ QLR ------")
    println("-----------------")
    cplexSolveQuadratricRelaxed()

end


function random_test()
    global solveTime, e1_size, e2_size, total_celles_not_cut, total_boundary_edges, obj_val
    global m,n,t

    record_grid = []
    models = ["MIP", "QLR"]
    record_times = Dict(m => [] for m in models)
    record_e1 = Dict(m => [] for m in models)
    record_e2 = Dict(m => [] for m in models)
    record_uncut = Dict(m => [] for m in models)
    record_boundaries = Dict(m => [] for m in models)
    record_objvalue = Dict(m => [] for m in models)

    for grid in 10:20
        m = grid
        n = grid
        t = [[] for _ in 1:m]

        for i in 1:m
            t[i] = [rand(60:100) for _ in 1:m]
        end

        append!(record_grid, grid)

        println("\n\n-----------------")
        println("------ P1 -------")
        println("-----------------")
        cplexSolveP1()

        append!(record_times["MIP"], solveTime)
        append!(record_e1["MIP"], e1_size)
        append!(record_e2["MIP"], e2_size)
        append!(record_uncut["MIP"], total_celles_not_cut)
        append!(record_boundaries["MIP"], total_boundary_edges)
        append!(record_objvalue["MIP"], obj_val)



        println("\n\n-----------------")
        println("------ QLR ------")
        println("-----------------")
        cplexSolveQuadratricRelaxed()

        append!(record_times["QLR"], solveTime)
        append!(record_e1["QLR"], e1_size)
        append!(record_e2["QLR"], e2_size)
        append!(record_uncut["QLR"], total_celles_not_cut)
        append!(record_boundaries["QLR"], total_boundary_edges)
        append!(record_objvalue["QLR"], obj_val)
    end

    for m in models
        plt.plot(record_grid, record_times[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the computing time between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("Computation time(s)", fontsize=12)
    savefig("res/times.png")
    plt.close()


    for m in models
        plt.plot(record_grid, record_e1[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the size of specie e1 between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("The number of specie e1", fontsize=12)
    savefig("res/e1.png")
    plt.close()


    for m in models
        plt.plot(record_grid, record_e2[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the size of specie e2 between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("The number of specie e2", fontsize=12)
    savefig("res/e2.png")
    plt.close()

    
    for m in models
        plt.plot(record_grid, record_uncut[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the total uncut parcels between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("The number of uncut parcels", fontsize=12)
    savefig("res/uncut.png")
    plt.close()


    for m in models
        plt.plot(record_grid, record_boundaries[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the boundary edges between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("The number of boundary edges", fontsize=12)
    savefig("res/boundary.png")
    plt.close()

    for m in models
        plt.plot(record_grid, record_objvalue[m], label = m, marker="o", linewidth=1.0, linestyle="--")
    end

    plt.legend(loc="upper left", fontsize=8)
    title("Comparaison the objective value between two models", fontsize=12)
    xlabel("Grid n x n", fontsize=12)
    ylabel("Objective value", fontsize=12)
    savefig("res/obj_value.png")
    plt.close()
    
end