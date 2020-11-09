using Pkg
Pkg.add("GAMS")
Pkg.add("JuMP")

using GAMS, JuMP
model = Model(GAMS.Optimizer)

@variable(model, 0 <= x <= 2)
@variable(model, 0 <= y <= 30)

@objective(model, Max, 5x + 3 * y)

@constraint(model, con, 1x + 5y <= 3)

optimize!(model)

set_optimizer_attribute(model, GAMS.Solver(), "CPLEX")
