using JuMP, Gurobi
function sample_problem(h)
    m = Model(solver = GurobiSolver(LogToConsole=0))
    ##
    # @variable(m, 0 <= x <= 2)
    # @variable(m, 0 <= y <= 30)
    @variable(m, 0 <= x <= h[1])
    @variable(m, 0 <= y <= h[2])

    @objective(m, Max, 5x + 3y)
    @constraint(m, 1x + 5y <= h[3])

    print(m)
    ##
    status = solve(m)

    println("Objective value: ", getobjectivevalue(m))
    println("x = ", getvalue(x))
    println("y = ", getvalue(y))
    return getobjectivevalue(m)
end
