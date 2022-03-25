TOL = 0.00001

using Random
using PyPlot
import PyPlot; const plt = PyPlot


global exploredNodes, solveTime, LB, E

function cplexSOlve()
    global exploredNodes, solveTime, LB, E

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

    @constraint(M, [i in 1:N], x[i] <= 3)

    @objective(M, Min, sum(p))

    # solve the problem
    set_silent(M) # turn off cplex output
    optimize!(M)

    exploredNodes = MOI.get(backend(M), MOI.NodeCount())
    println("exploredNodes : ", exploredNodes)
    solveTime = round(MOI.get(M, MOI.SolveTime()), digits = 3)
    println("solveTime : ", solveTime)

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    println("isOptimal ? ", isOptimal)
    if isOptimal
        println("x = ", round.(Int, value.(M[:x])))

        println("Lower bound = ", round(objective_value(M), digits = 6))
        LB = round(objective_value(M), digits = 6)

        println("Log probabilites : ", round.(value.(M[:p]), digits = 6))

        proba = zeros(G, A)
        for g in 1:G
            for a in 1:A
                if JuMP.value(p[g, a]) > TOL
                    tmp = [JuMP.value(x[i]) for i in 1:N if sum(individu[i][1][g])%2==1]
                    if size(tmp, 1) == 0
                        proba[g, a] = 1
                    else
                        proba[g, a] = (0.5)^(sum(JuMP.value(x[i]) for i in 1:N if sum(individu[i][1][g])%2==1))
                    end
                    
                end
            end
        end

        println("Real probability = ", round.(proba, digits = 6))

        println("True expected value = ", round(sum(proba), digits = 6))
        E = round(sum(proba), digits = 6)
    else
        error("cplex cannot find opt sol !")
    end

end


function run()
    include("DivGenetique_opl.dat")
    cplexSOlve()

end


function analysis_picewises()
    include("DivGenetique_opl.dat")
    global H
    global exploredNodes, solveTime, LB, E

    fout = open("res/picewises.txt", "w")

    record_H = []
    record_E = []
    record_LB = []

    for pw in 1:20
        H = pw * 10
        println("H = ", H)
        cplexSOlve()

        append!(record_H, H)
        append!(record_LB, LB)
        append!(record_E, E)

        println(fout, H, " & ", exploredNodes, " & ", solveTime, " & ", LB, " & ", E, " \\\\")
    end

    close(fout)

    plt.plot(record_H, record_LB, label = "LB", marker="o", linewidth=1.0, linestyle="--")
    plt.plot(record_H, record_E, label = "E", marker="o", linewidth=1.0, linestyle="--")


    plt.legend(loc="lower right", fontsize=8)
    title("Impact of the number of picewises", fontsize=10)
    xlabel("number of picewises", fontsize=10)
    ylabel("Approximated / admissible expected value", fontsize=10)
    savefig("res/picewises.png")
    plt.close()

end



function analysis_population()
    include("DivGenetique_opl.dat")
    global exploredNodes, solveTime, LB, E
    global N, Nm, Nf, individu, G

    fout = open("res/population_size.txt", "w")

    for pop in 1:5
        N = 2 * pop
        Nm = pop
        Nf = pop

        println("N = ", N)
        seed = rand(1:10000)
        Random.seed!(seed) 
        individu = [ [[ [rand(1:2), rand(1:2)] for _ in 1:G]] for i in 1:N]

        cplexSOlve()

        println(fout, N, " & ", exploredNodes, " & ", solveTime, " & ", LB, " & ", E, " \\\\")

    end

    close(fout)

end


function analysis_genes()
    include("DivGenetique_opl.dat")
    global exploredNodes, solveTime, LB, E
    global N, individu, G

    fout = open("res/genes_number.txt", "w")

    for gene in 2:15
        G = gene

        println("G = ", G)

        seed = rand(1:10000)
        Random.seed!(seed) 
        individu = [ [[ [rand(1:2), rand(1:2)] for _ in 1:G]] for i in 1:N]

        cplexSOlve()

        println(fout, G, " & ", exploredNodes, " & ", solveTime, " & ", LB, " & ", E, " \\\\")

    end

    close(fout)


end