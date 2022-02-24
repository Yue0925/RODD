TOL = 0.00001

using Random
using PyPlot


function cplexSolve(output::String)

    fout = open(output, "w")

    lines = size(proba, 1)
    probabilities = zeros(k, n, m)
    for l in 1:lines
        probabilities[round(Int, proba[l, 1]), round(Int, proba[l, 2]), round(Int, proba[l, 3])] = proba[l, 4]
    end

    M = Model(CPLEX.Optimizer) 

    @variable(M, z[1:n, 1:m], Bin) # z[i, j] = 1, if cell(i, j) selected as zone central
    @variable(M, x[1:n, 1:m], Bin) # x[i, j] = 1, if cell(i, j) is protected
    
    # apreantly, zone central cannot situates on the bord
    @constraint(M, sum(z[1, 1:m]) == 0)
    @constraint(M, sum(z[n, 1:m]) == 0)
    @constraint(M, sum(z[1:n, 1]) == 0)
    @constraint(M, sum(z[1:n, m]) == 0)


    # a central plot is surrounded by protected plots
    @constraint(M, [i in 2:n-1, j in 2:m-1], z[i, j] >= sum(x[i_, j_] for i_ in i-1:i+1, j_ in j-1:j+1) - 8)

    # a plot cannot be selected as central zone if at least one of its neighbour is not protected
    for i in 2:n-1
        for j in 2:m-1
            @constraint(M, [i_ in i-1:i+1, j_ in j-1:j+1], z[i, j] <= x[i_, j_])
        end
        
    end

    
    # the survival probability constraint
    @constraint(M, [e in 1:p], sum(log(1 - probabilities[e, i, j]) * z[i, j] for i in 1:n, j in 1:m) <= log(1-alpha[e]))
    @constraint(M, [e in p+1:k], sum(log(1 - probabilities[e, i, j]) * x[i, j] for i in 1:n, j in 1:m) <= log(1-alpha[e]))


    @objective(M, Min, sum(c[i][j] * x[i, j] for i in 1:n, j in 1:m))

    # resolve
    optimize!(M)

    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    GAP = MOI.get(M, MOI.RelativeGap())
    solveTime = MOI.get(M, MOI.SolveTime())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL


    probability = zeros((k))

    println(fout, "isOptimal ? ", isOptimal)
    println(fout, "GAP : ", GAP)
    println(fout, "exploredNodes : ", exploredNodes)
    println(fout, "solveTime : ", round(solveTime, digits = 2))
    println(fout, "objective value : ", round(Int, objective_value(M)))

    if has_values(M)

        # proba of each espèce
        for e in 1:p
            π = 1
            for i in 1:n
                for j in 1:m
                    if JuMP.value(z[i, j]) > TOL
                        π *= 1 - probabilities[e, i, j]
                    end
                end
            end
            probability[e] = round(1 - π, digits = 2)
        end

        for e in p+1:k
            π = 1
            for i in 1:n
                for j in 1:m
                    if JuMP.value(x[i, j]) > TOL
                        π *= 1 - probabilities[e, i, j]
                    end
                end
            end
            probability[e] = round(1 - π, digits = 2)
        end

        println(fout, "probability : ", probability)
        println(fout)

        # display
        for i in 1:n
            for j in 1:m
                if JuMP.value(z[i, j]) > TOL
                    print(fout, '■')
                elseif JuMP.value(x[i, j]) > TOL
                    print(fout, '▪') 
                else
                    print(fout, '□')
                end
            end
            println(fout)
        end
    end

    close(fout)
    return solveTime, exploredNodes
end



function run()
    file = "ReserveNaturelles.txt"
    include(file)

    output = "res/cas4"
    cplexSolve(output)
end


function run_random_tests()
    Random.seed!(1234) 

    global record_size = [i*10 for i in 1:20]
    global record_time = []
    global record_nodes = []

    for i in 1:20
        global n = i * 10
        global m = n
        global p=3
        global q=3
        global k = 6
        global alpha = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        global proba = zeros(k*n, 4)
        l=0
        for e in 1:k
            for c in 1:n
                l +=1
                proba[l, 1] = e
                proba[l, 2] = rand(1:n)
                proba[l, 3] = rand(1:m)
                proba[l, 4] = round(rand(1:9) * 0.1, digits = 2)
            end
        end

        global c = [[] for _ in 1:n]
        for i in 1:n
            c[i] = [rand(1:10) for _ in 1:m]
        end
        instance = "instance$n" * "x$n"

        output = "res/random/" * instance
        solveTime, exploredNodes = cplexSolve(output)
        append!(record_time, solveTime)
        append!(record_nodes, exploredNodes)
    end

    plot(record_size, record_time)
    title("Computation time according to instances size")
    xlabel("size")
    ylabel("Time(s)")
    savefig("res/random/" * "time.png")
    close()


    plot(record_size, record_nodes)
    title("The number of explored nodes according to instances size")
    xlabel("size")
    ylabel("nodes")
    savefig("res/random/" * "nodes.png")
    close()
end