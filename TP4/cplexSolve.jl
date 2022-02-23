

function cplexSolveQuadratricRelaxed()

    include("ExplForet_opl.dat")

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


    # solve the problem
    optimize!(M)

    
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL


    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = MOI.get(M, MOI.SolveTime())

    println(" exploredNodes ", exploredNodes)
    println(" solveTime ", solveTime)

    # display solution
    println("isOptimal ? ", isOptimal)
    println("GAP = ", GAP)


    if has_values(M)
        obj_val = objective_value(M)
        println("obj_val : ", obj_val)

        # println(JuMP.value.(y[1, :, :, :]))
        println(value.(x))
        println("the number of celles not cut : ", sum(JuMP.value.(x)))
    end

    
end




function cplexSolveP1()

    include("ExplForet_opl.dat")

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


    # solve the problem
    optimize!(M)
    println(solution_summary(M))
    
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL


    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = MOI.get(M, MOI.SolveTime())

    println(" exploredNodes ", exploredNodes)
    println(" solveTime ", solveTime)

    # display solution
    println("isOptimal ? ", isOptimal)
    println("GAP = ", GAP)


    if has_values(M)
        obj_val = objective_value(M)
        println("obj_val : ", obj_val)

        println("JuMP.value.(x) : ", JuMP.value.(x))
        println("the number of celles not cut : ", sum(JuMP.value.(x)))
    end

    
end

