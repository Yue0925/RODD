TOL = 0.00001


function cplexSOlve()

    M = Model(CPLEX.Optimizer) 

    @variable(M, x[1:N] >= 0, Int) # the numbre of chid
    @variable(M, 0 <= p[1:G, 1:A] <=1) # the proba of each allèle
    @variable(M, 0 <= t[1:G, 1:A] <= 1) # the proba

    @constraint(M, [g in 1:G, a in 1:A], p[g, a] >= t[g, a] - sum(x[i] for i in 1:N if individu[i][1][g] == [a, a]))

    θ_r(r) = init^((H-r)/(H-1))

    @constraint(M, [r in 1:H, g in 1:G, a in 1:A], log(θ_r(r)) + (1/θ_r(r)) * (t[g, a] - θ_r(r) ) >=
        sum(x[p] * log(0.5) for p in 1:N if sum(individu[p][1][g])%2==1))

    @constraint(M, sum(x[i] for i in 1:Nm) == N)
    @constraint(M, sum(x[i] for i in N-Nf+1:N) == N)

    @constraint(M, [i in 1:N], x[i] <= 2)

    @objective(M, Min, sum(p))

    # solve the problem
    optimize!(M)

    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    solveTime = MOI.get(M, MOI.SolveTime())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    println("isOptimal ? ", isOptimal)
    if isOptimal
        println("x = ", value.(M[:x]))

        println("Lower bound = ", objective_value(M))

        println("Log probabilites : ", value.(M[:p]))

        proba = zeros(G, A)
        for g in 1:G
            for a in 1:A
                if JuMP.value(p[g, a]) > TOL
                    proba[g, a] = (0.5)^(sum(JuMP.value(x[i]) for i in 1:N if sum(individu[i][1][g])%2==1))
                end
            end
        end

        println("Real probability = ",proba)

        println("True expected value = ", sum(proba))
    end

end