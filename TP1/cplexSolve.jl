TOL = 0.00001


function cplexSolve()
    include("ReserveNaturelles.txt")
    lines = size(proba, 1)
    probabilities = zeros(6, 10, 10)
    for l in 1:lines
        probabilities[round(Int, proba[l, 1]), round(Int, proba[l, 2]), round(Int, proba[l, 3])] = proba[l, 4]
    end

    M = Model(CPLEX.Optimizer) 

    @variable(M, x[1:6, 1:10, 1:10], Bin) # x[k, i, j] = 1, if espece k is in the cell(i, j)
    @variable(M, yc[1:10, 1:10], Bin) # yc[i, j] = 1, if cell(i, j) selected as zone central
    @variable(M, yt[1:10, 1:10], Bin) # yt[i, j] = 1, if cell(i, j) selected as zone tampon


    # each cell can be selected as at most one zone
    for i in 1:10
        for j in 1:10
            @constraint(M, yc[i, j] + yt[i, j] <= 1) 
        end
    end
    
    # apreantly, zone central cannot situates on the bord
    @constraint(M, sum(yc[1, 1:10]) == 0)
    @constraint(M, sum(yc[10, 1:10]) == 0)
    @constraint(M, sum(yc[1:10, 1]) == 0)
    @constraint(M, sum(yc[1:10, 10]) == 0)

    # each dangerours espece belongs to at least one zone
    for i in 1:6
        @constraint(M, sum(x[i, 1:10, 1:10]) >= 1)
    end

    for i in 1:10
        for j in 1:10
            # for each region
            @constraint(M, 3 * yc[i, j] >= sum(x[k, i, j] for k in 1:3)) # dangerours espece
            @constraint(M, 3 * (yc[i, j] + yt[i, j]) >= sum(x[k, i, j] for k in 4:6)) #commune espece
            # for k in 1:3
            #     @constraint(M, yc[i, j] >= x[k, i, j])
            # end
            # for k in 4:6
            #     @constraint(M, yc[i, j] + yt[i, j] >= x[k, i, j])
            # end
        end
    end
    
    # the survival probability constraint
    for k in 1:6
        @constraint(M, sum(log(1 - probabilities[k, i, j]) *  x[k, i, j] for i in 1:10 for j in 1:10) <= log(1 - alpha[k]))
    end

    # the geographical constraint
    for i in 2:9
        for j in 2:9
            @constraint(M, yc[i, j] >= sum(yc[p, q] + yt[p, q] for p in i-1:i+1 for q in j-1:j+1 if p!=i && q!=j) - 7)
        end
    end

    @objective(M, Min, sum(c[i][j] * (yc[i, j]+yt[i, j]) * 10 for i in 1:10 for j in 1:10))

    start = time()
    # solve the problem
    optimize!(M)

    computationTime = time() - start
    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    centrals = Array{Tuple{Int64, Int64}, 1}()
    tampons = Array{Tuple{Int64, Int64}, 1}()

    println("isOptimal ? ", isOptimal)
    if isOptimal
        println("objective value : ", objective_value(M))
        for i in 1:10
            for j in 1:10
                if JuMP.value(yc[i, j]) > TOL
                    append!(centrals, (i, j))
                end
                if JuMP.value(yt[i, j]) > TOL
                    append!(tampons, (i, j))
                end

            end
        end
    end
end