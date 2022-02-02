TOL = 0.00001

# distance euclidienne
dst(i, j, i_, j_) = sqrt( (i-i_)^2 + (j-j_)^2)

"""
Using cplex solver resolve P_λ problem
"""
function cplexSolve(λ)


    # modelization
    M = Model(CPLEX.Optimizer) 

    # variables 
    @variable(M, x[1:n, 1:m], Bin)
    @variable(M, d[1:n, 1:m, 1:n, 1:m], Bin)
    @variable(M, g)
    @variable(M, f)


    # constraint of air total
    @constraint(M, sum(x[i, j] for i in 1:n for j in 1:m) <= Amax)
    @constraint(M, sum(x[i, j] for i in 1:n for j in 1:m) >= Amin)

    # buget constraint
    @constraint(M, sum(c[i][j] *10 * x[i, j] for i in 1:n for j in 1:m) <= B)

    # the neighbourhood constraint
    for i in 1:n
        for j in 1:m
            @constraint(M, sum(d[i, j, i_, j_] for i_ in 1:n for j_ in 1:m if i!=i_ || j!=j_) == x[i, j])
            for i_ in 1:n
                for j_ in 1:m
                    if i==i_ && j==j_
                        continue
                    end
                    @constraint(M, d[i, j, i_, j_] <= x[i_, j_])
                end
            end
        end
    end

    # objective function is to maximize λ*g(x) - f(d)
    @constraint(M, g == sum(x[i, j] for i in 1:n for j in 1:m))
    @constraint(M, f == sum(dst(i, j, i_, j_) * d[i, j, i_, j_] for i in 1:n for j in 1:m for i_ in 1:n for j_ in 1:m if i!=i_ || j!=j_))
    @objective(M, Min, f - λ * g)

    # start a chronometer
    start = time()

    # solve the problem
    optimize!(M)

    computationTime = time() - start
    exploredNodes = MOI.get(backend(M), MOI.NodeCount())

    # status of model
    status = termination_status(M)
    isOptimal = status==MOI.OPTIMAL

    println("isOptimal ? ", isOptimal)
    println("time(s) = ", computationTime)
    println("exploredNodes = ", exploredNodes)

    v_λ = 0
    x_λ = zeros(n, m)
    d_λ = zeros(n, m, n, m)
    g_star = 0.0
    f_star = 0.0

    if isOptimal
        println("v(lambda) = ", objective_value(M))
        v_λ = objective_value(M)

        g_star = JuMP.value(g)
        f_star = JuMP.value(f)
        println("f(d) = ", f_star)
        println("g(x) = ", g_star)

        for i in 1:n
            for j in 1:m
                if JuMP.value(x[i, j]) > TOL
                    x_λ[i, j] = 1.0
                end
                for i_ in 1:n
                    for j_ in 1:m
                        if JuMP.value(d[i, j, i_, j_]) > TOL
                            d_λ[i, j, i_, j_] = 1.0
                        end 
                    end
                end
            end
        end
    end

    return isOptimal, v_λ, x_λ, d_λ, f_star, g_star
end


function Dinkelbach()
    include("MinFragmentation_opl.dat")
    
    # step 1
    λ = lambda
    ite = 1

    # step 2
    isOptimal, v_λ, x_λ, d_λ, f_star, g_star = cplexSolve(λ)

    while isOptimal && ( v_λ > TOL || v_λ < -TOL)
        ite += 1
        #λ = sum(dst(i, j, i_, j_) * d_λ[i, j, i_, j_] for i in 1:n for j in 1:m for i_ in 1:n for j_ in 1:m)/sum(x_λ[i, j] for i in 1:n for j in 1:m)
        λ = f_star/g_star
        isOptimal, v_λ, x_λ, d_λ, f_star, g_star = cplexSolve(λ)
    end 

    println("End of iteration with ", ite)
    println("v(lambda) = ", v_λ, " with lambda = ", λ)
    println("the total number celles selected = ", sum(x_λ))
    println("f(d) = ", f_star)
    println("g(x) = ", g_star)
end