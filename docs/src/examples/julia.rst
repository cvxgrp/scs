.. _julia_example:

Julia
======

The easiest way to use SCS in Julia is via JuMP:

.. code:: julia

  using JuMP, SCS
  items  = [:Gold, :Silver, :Bronze]
  values = Dict(:Gold => 5.0,  :Silver => 3.0,  :Bronze => 1.0)
  weight = Dict(:Gold => 2.0,  :Silver => 1.5,  :Bronze => 0.3)
  model = Model(SCS.Optimizer)
  @variable(model, 0 <= take[items] <= 1)  # Define a variable for each item
  @objective(model, Max, sum(values[item] * take[item] for item in items))
  @constraint(model, sum(weight[item] * take[item] for item in items) <= 3)
  optimize!(model)
  println(value.(take))
  # 1-dimensional DenseAxisArray{Float64,1,...} with index sets:
  #     Dimension 1, Symbol[:Gold, :Silver, :Bronze]
  # And data, a 3-element Vector{Float64}:
  #  1.0000002002226671
  #  0.4666659513182934
  #  1.0000007732744878

Read the `JuMP documentation <https://jump.dev/JuMP.jl/stable/>`_ and the 
`SCS.jl README <https://github.com/jump-dev/SCS.jl>`_ for more details.
