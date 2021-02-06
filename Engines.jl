### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ e145c39e-624a-11eb-392e-25c51ea78395
begin
	using PlutoUI, Dates, Gurobi, Cbc, Clp, GLPK, JuMP, Random, BenchmarkTools, DataStructures, LinearAlgebra, Plots
	project_name = "Optimization engines" 
	date = Dates.Date(Dates.now())
	company = "Rappi"
	objective = "to determine whether to acquire a commercial solver or develop in-house algorithms."
end;

# ╔═╡ 345a1756-624b-11eb-0e73-b99c01d7852d
md"""
# $project_name @ $company

The 
**objective** is $objective 

Last time modified: $date

"""

# ╔═╡ 3e591e1e-624b-11eb-3a49-e1c420a8a740
PlutoUI.TableOfContents()

# ╔═╡ e3877146-6495-11eb-3fde-bd2e6806a7ef
md"""
# Implementation

The VRP has been selected because represents one of the hardest for the operational process at Rappi.
"""

# ╔═╡ 750b95a0-6407-11eb-15b5-8b4a9805b7e8
md"""
## Capacitated VRP with Time Windows
"""

# ╔═╡ 505933d2-66f9-11eb-1cd1-d56e8c642c12
md"""
First, let's take a look at the most important structures of the model
"""

# ╔═╡ 384c8206-66f9-11eb-22bc-e3b7d3a8e05b
mutable struct Node
	ID   ::Int64
	lat  ::Float64
	long ::Float64
end

# ╔═╡ 3c12dd48-66f9-11eb-1d10-ad9884fa5827
mutable struct RT
	ID              ::Int64
	capacity        ::Int64
	max_travel_time ::Float64
end;

# ╔═╡ 438b08d6-66f9-11eb-2e44-256ad2bf5b11
mutable struct Order
    ID           ::Int64
	delivery     ::Node
    early_time   ::Float64
    late_time    ::Float64
    items_qty    ::Int64
    service_time ::Int64
end;

# ╔═╡ 47554dd2-66f9-11eb-1575-717f4997cec5
mutable struct Instance
    travel_times      ::Matrix{Float64}
    service_times     ::Vector{Int64}
    early_times       ::Vector{Int64}
    late_times        ::Vector{Int64}
    items_orders      ::Vector{Int64}
    RT_capacity       ::Int64
    RT_max_route_time ::Float64
end;

# ╔═╡ 4b27ed3e-66f9-11eb-011b-3d3e1427cd49
mutable struct Solution
    routes    ::Dict{Any,Any}
    objective ::Float64
end;

# ╔═╡ 94bbeeae-6407-11eb-2bf5-a510e938453c
md"""
We are given:

- Set of $M = \left \{ 1, ..., m \right \}$ clients
- 


**Assumptions**
- There are enought RTs to fulfill all orders: $Q_{orders} == Q_{RTs}$
- Homogenous RTs.
- Operational time: travel time + service time



**Decision variables**
- Binary $X_{i}$ indicates whether ...

**Objective**
The objective is to minimize the total . This cost breaks down into two components:
- Fixed cost of ...
- Cost of ...


$ MINIMIZE () $

**Contraints**
- Each ...

Ver comentarios de data generation!


"""

# ╔═╡ 0657a1be-66ad-11eb-233d-15f3f93307e4
md"""
## Data generation
"""

# ╔═╡ 5ca26cdc-66fa-11eb-0f4b-334f9ca51c80
distance(n1::Node, n2::Node) = sqrt((n1.lat - n2.lat)^2 + (n1.long - n2.long)^2)

# ╔═╡ 86141c1e-66fa-11eb-3adf-ffb61115f195
function calculate_time(nodes::Vector{Node}; minutes_constant = 3, digits = 2)
	
    n_nodes = length(nodes)
	store = n_nodes - 1 # start or source
	dummy = n_nodes # end or sink
    time = zeros(n_nodes, n_nodes)
	
    for i in 1:n_nodes-1
        for j in i+1:n_nodes
            d = distance(nodes[i], nodes[j])
            factor = 10^digits
            minutes = minutes_constant * (floor(factor * d) / factor)
            time[i,j] = minutes
            time[j,i] = minutes
        end
    end
	
    time[store, dummy] = 0
    time[dummy, store] = 0
	
    return time
	
end

# ╔═╡ c8db0efe-66fa-11eb-2c99-6f2089a20cf1
"""RBN => Random by now"""
function generate_data(size::Int64)#::Instance
	
	# Time windows
	# -----
	
	early_times = rand(1:50, size) # minutes until order gets ready + minutes from 										 store to customer's location. RBN.
	
	late_times = [x+rand(5:30) for x in early_times] # early + maximum time order 														   gets "cold". RBN.

	# Capacity
	# -----
	
	items_orders = rand(1:4, size) # items each order have or space volume of the 										 order. RBN.
	
	heaviest_order,  = maximum(items_orders)
	
	avg_weight_order = ceil(mean(items_orders))
	
	capacity = rand(heaviest_order:heaviest_order+avg_weight_order) # RT visits n																	      customers.

	# Locations
	# -----
	
	stores = repeat([Node(0, rand(0:.1:10), rand(0:.1:10))], 2) # store + dummy
	
	customers = [Node(i, rand(0:.1:10), rand(0:.1:10)) for i in 1:size] 
	
	nodes = vcat(customers, stores) # stores at the end

	# Travel and service times
	# -----

	travel_times = calculate_time(nodes, minutes_constant = 3)
	
	service_times = rand(1:5, size) # dist historical waiting time at store. RBN.

	worst_travel, avg_travel = maximum(travel_times), ceil(mean(travel_times))
	
	max_travel_time = rand(worst_travel:worst_travel+avg_travel)
	
	
	return Instance(travel_times, service_times, early_times, late_times, 
					items_orders, capacity, max_travel_time)
	
end

# ╔═╡ d4d31cba-66fa-11eb-2a5b-593e62edf09d
function instance_OR_tools()::Instance
	
	size = 16
	
	time_windows = [
		(7, 12),    # customer  1
		(10, 15),   # customer  2
		(16, 18),   # customer  3
		(10, 13),   # customer  4
		(0, 5),     # customer  5
		(5, 10),    # customer  6
		(0, 4),     # customer  7
		(5, 10),    # customer  8
		(0, 3),     # customer  9
		(10, 16),   # customer 10
		(10, 15),   # customer 11
		(0, 5),     # customer 12
		(5, 10),    # customer 13
		(7, 8),     # customer 14
		(10, 15),   # customer 15
		(11, 15),   # customer 16
	]
	
	early_times = [tw[1] for tw in time_windows]
	late_times = [tw[2] for tw in time_windows]
	
	coordinates = [
    (2.28, 0.00),   # location  1
    (9.12, 0.00),   # location  2
    (0.00, 0.80),   # location  3
    (1.14, 8.00),   # location  4
    (5.70, 1.60),   # location  5
    (7.98, 1.60),   # location  6
    (3.42, 2.40),   # location  7
    (6.84, 2.40),   # location  8
    (5.70, 4.00),   # location  9
    (9.12, 4.00),   # location 10
    (1.14, 4.80),   # location 11
    (2.28, 4.80),   # location 12
    (3.42, 5.60),   # location 13
    (6.84, 5.60),   # location 14
    (0.00, 6.40),   # location 15
    (7.98, 6.40),   # location 16
    (4.56, 3.20),   # location  0 - the depot
    (4.56, 3.20),   # location  0 - the depot dummy
	]
	
	nodes = [Node(1, coord[1], coord[2]) for coord in coordinates]
	travel_times = calculate_time(nodes, minutes_constant = 1)
	service_times = zeros(size)
	items_orders = [1, 1, 2, 4, 2, 4, 8, 8, 1, 2, 1, 2, 4, 4, 8, 8]
	capacity = 15
	max_travel_time = 35  
	
	return Instance(travel_times, service_times, early_times, late_times, 								items_orders, capacity, max_travel_time
					)
	
end

# ╔═╡ b860657c-67d6-11eb-0240-6b84b814c3a9
md"""
## Solution approaches

We will explore different approaches for solving the problem:
- OR-tools
- Mathematical model
- Heuristics

"""

# ╔═╡ 33674ca8-67d5-11eb-1d83-89ed81652979
md"""
### OR-tools
"""

# ╔═╡ 403d1c28-67d5-11eb-3a0f-e92303f07e3f
function or_tools_routing(instance::Instance)::Solution
	
	
end

# ╔═╡ f7c7f9c0-6632-11eb-27bb-e7a49dde68b8
md"""
### Mathematical model
"""

# ╔═╡ 2f455ffe-6634-11eb-36b5-d7d9c8e2decf
md"""
We will benchmark the following solvers:

- [Gurobi](https://www.gurobi.com): commercial solver for both LP and MILP.
- [Clp](https://github.com/coin-or/Clp): open-source for LP problems (COIN-OR).
- [Cbc](https://github.com/coin-or/Cbc): open-source for MILP problems (COIN-OR).
- [GLPK](https://www.gnu.org/software/glpk/): open source.
- [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio): IBM commercial solver for both LP and MILP.

Using [JuMP](https://github.com/jump-dev/JuMP.jl) as the mathematical modeling language.
"""

# ╔═╡ 804bcac8-6706-11eb-1d8a-756a5afde359
function format_routes(X, store::Int64, dummy_store::Int64, 												total_time::Matrix{Float64})
	
    x_val = JuMP.value.(X).data
    first_nodes = findall(x -> x > 0.9, x_val[store, :])
    routes_times = Dict()

    for i in first_nodes

        next = findfirst(x -> x > 0.9, x_val[i, :])
        route = [store, i, next]
        cc = total_time[store, i] + total_time[i, next]

        while next != dummy_store
			
            current = next
            next = findfirst(x -> x > 0.9, x_val[current, :])
            push!(route, next)
            cc += total_time[current, next]
			
        end
		
        routes_times[route] = cc
		
    end
	
	return routes_times
	
end

# ╔═╡ 48699cb4-6366-11eb-175d-8728ce809fea
function calculate_total_time(
		travel_times::Matrix{Float64}, service_times::Vector{Int64}, 						LC::UnitRange{Int64}, DL::UnitRange{Int64})

	for i in LC, j in LC
		if i in DL
			travel_times[i, j] += service_times[i]
		end
	end
	
	return travel_times
	
end

# ╔═╡ cda9bb12-6671-11eb-0494-476f9966741d
function optimal_routing(instance::Instance, solver::Module, solver_params::Pair...)
    
    # Parameters
    # -----

	n_locations = Base.size(instance.travel_times, 1) # delivery location for each order (customer's place) + store (source) + dummy (sink)
	LC = 1:n_locations
    
	n_deliveries = length(instance.service_times) # delivery location for each order (customer's place)
    DL = 1:n_deliveries

    n_RTs = n_deliveries # RTs available
    RT = 1:n_RTs

    store = n_locations - 1 # penultimate node
    dummy_store = deepcopy(n_locations) # last node of the route

    items_orders, capacity, max_travel_time = instance.items_orders, instance.RT_capacity, instance.RT_max_route_time

	early_times = [instance.early_times; 0; 0]
	late_times = [instance.late_times; 0; max_travel_time]
    
    total_time = calculate_total_time(instance.travel_times, instance.service_times, LC, DL)
	
	bigT = max_travel_time * 2
    bigQ = capacity + maximum(items_orders)
    

    # Formulation
    # -----

	model = Model(optimizer_with_attributes(solver.Optimizer, solver_params...))

	@variable(model, X[LC, LC], Bin) # whether location (n,n) is selected for the route or not.
	@variable(model, arrival[LC] >= 0) # arrival time at each node in the route.
	@variable(model, units[LC] >= 0) # units or volume RT carries for the route.

	@objective(model, Min, sum(total_time[i,j] * X[i,j] for i in LC, j in LC)) # minimize the total route time.

	@constraint(model, [i in DL], sum(X[i,j] for j in LC) == 1) # each delivery location must be visited exactly once.
	
	@constraint(model, [i in LC, j in DL], units[j] >= units[i] + items_orders[j] - bigQ * (1 - X[i,j])) # if selected, ...
		
	@constraint(model, [i in LC], 0 <= units[i] <= capacity) # units or volume carried by RT should not exceed capacity.
		
	@constraint(model, sum(X[store, j] for j in LC) <= n_RTs) # store is the beginning for all the routes (not all RTs have to be used).
		
	@constraint(model, [h in DL], sum(X[i, h] for i in LC) - sum(X[h,j] for j in LC) == 0) # for each delivery location, ...

	@constraint(model, sum(X[i, dummy_store] for i in LC) == sum(X[store, j] for j in LC))
		
	@constraint(model, [i in LC, j in LC], arrival[j] >= arrival[i] + total_time[i,j] - bigT * (1 - X[i,j])) # if node selected, then the arrival time should ...
		
	@constraint(model, [i in LC], early_times[i] <= arrival[i] <= late_times[i]) # arrival at node should be between time window.
		
	@constraint(model, [i in LC], X[i, store] == 0)
		
	@constraint(model, [i in LC], X[dummy_store, i] == 0) 
		
	@constraint(model, [i in LC], X[i, i] == 0) # should not select the diagonal of the matrix.

    # Optimize!
    # -----

    optimize!(model)
    
    # status = JuMP.termination_status(model)		
	dict_routes_times = format_routes(X, store, dummy_store, total_time)
    
	return Solution(dict_routes_times, objective_value(model))
	
end

# ╔═╡ 97222556-6366-11eb-03b2-e78c4247ccdf
md"""
### Heuristics
"""

# ╔═╡ 971ff504-6633-11eb-06d4-09b91dac0ac6
md"""
In broad terms, a two-phase process is implemented:

1. **Construction heuristic**: generate an initial solution to start from (greedy algorithm)

2. **Local-search metaheuristic**: iteratively improve the initial solution by applying several move operators.


**Check**
- Engine.jl (local)
- [Nodal.jl](https://github.com/phrb/NODAL.jl)
"""

# ╔═╡ 8ee45d66-67ce-11eb-0309-6b5c68f6f09c
function construction(data)
	
end

# ╔═╡ ab90d342-67ce-11eb-338c-a155fa468e83
function move_operator!(data)

end

# ╔═╡ 9d7b483c-67ce-11eb-359c-a38abb43f688
function local_search(initial_solution)
	
end

# ╔═╡ b9a63d6e-63f8-11eb-111f-afb7f1ae88c9
md"""
# Analysis
"""

# ╔═╡ 1972e386-6733-11eb-2e45-f50f5240c102
md"""
## Metrics
"""

# ╔═╡ 7743b4e8-6735-11eb-3ee7-117ee0d6b494
md"""
##### Technical
"""

# ╔═╡ 7bd83a4c-67d8-11eb-0864-574fbb0bcc35
md"""
Time

- The **minimum** is a robust estimator for the location parameter of the time distribution, and should not be considered an outlier.

- The **median** is as a robust measure of central tendency, should be relatively unaffected by outliers.

- The **mean** is a non-robust measure of central tendency, will usually be positively skewed by outliers

- The **maximum** should be considered a primarily noise-driven outlier, and can change drastically between benchmark trials.

Allocs
- k

Memory
- 

Solution quality
- j

Consistency
- h
"""

# ╔═╡ 7f46fae2-6735-11eb-3042-4360a215e0b1
md"""
##### Business
"""

# ╔═╡ 89b36da2-67d8-11eb-1bab-5f474524e5f5
md"""
Price
- average price found in website.

Long-term contract
- ...

Implementation cost

$C_{impl} = months \times people \times salary$
"""

# ╔═╡ 353c2780-6733-11eb-3bd0-25bf2e882dd8
engine_metrics = Dict("Gurobi" => Dict(), "Clp" => Dict(), "Cbc" => Dict(), 
					  "GLPK" => Dict(), "CPLEX" => Dict(), "Heuristic" => Dict());

# ╔═╡ 48f995ce-6633-11eb-0556-bf3da4f6392b
md"""
## Benchmark

It's time to technically analyze the performance of each solver. 

The idea is to **incrementally add new experiments and compare the individual technical metrics of each run and the cumulative results between engines**

Let's define:

- **Sample** as a complete measurement for a single solver. **Each sample is a set of ~200 individual runs**: `scenario x in for Gurobi`

- **Experiment** as a set of samples, one for each solver: `scenario x in Gurobi, Cbc and GLPK`.

Therefore, we have to define different scenarios with **different datasets and solver's parameters configurations**. 

By doing so, we will be able to measure **repeatability, reproducibility and consistency**.

The experiments were executed locally in a MacBook Pro (16-inch, 2019).
- Processor: 2.3 GHz 8-Core Intel Core i9
- Memory: 16 GB 2667 MHz DDR4
- Graphics: Intel UHD Graphics 630 1536 MB

"""

# ╔═╡ c44cb24c-6726-11eb-30fe-bf4ab0bf02b4
"""Given set of parameters, perform the benchmark for each solver"""
function benchmark_engines(data::Instance, engines::Vector{Module}, all_params::Any)

	suite = BenchmarkGroup(["engines"])

	for engine in engines
		engine_name = string(engine)
		e_params = all_params[engine_name]

		suite[engine_name] = @benchmark optimal_routing($data,$engine,$e_params...)
	end

	return suite
	
end

# ╔═╡ 215eb286-67c5-11eb-1a2e-fba123258c06
function add_experiment!(engine_metrics::Dict{String,Dict{Any,Any}}, 												experiment::BenchmarkGroup, engines::Vector{Module})

	for engine in string.(engines)
		
		allocs_estimated = experiment[engine].allocs
		memory_estimated = experiment[engine].memory
		
		if length(collect(engine_metrics[engine])) == 0
			engine_metrics[engine]["Allocs"] = [allocs_estimated]
			engine_metrics[engine]["Memory"] = [memory_estimated]
		else # not first time
			push!(engine_metrics[engine]["Allocs"], allocs_estimated)
			push!(engine_metrics[engine]["Memory"], memory_estimated)
			
		end
		
		for func in [minimum, median, mean, maximum]
			
			measure = func(experiment[engine]).time
			name = "Time-"*string(func)
						
			if name in collect(keys(engine_metrics[engine])) # not first time
				push!(engine_metrics[engine][name], measure)
			else
				engine_metrics[engine][name] = [measure]
			end
			
		end
		
	end
	
	return engine_metrics
	
end

# ╔═╡ 498502dc-6726-11eb-23f5-4151f147db6f
md"""

### Scenarios

Run an experiment for each scenario:

- Small instance until optimality (no solver params)
- Hard instance with 100 seconds execution limit (~1.6 min)
- Hard instance with .5 optimality gap
- ...
	
"""

# ╔═╡ 1c8a6d54-6725-11eb-0025-85cacfc3b2ab
md"""
### Instance selection

For randomly generated datasets, it's required to choose the **size**.

Delivery locations = $(@bind size Slider(1:300, default=10, show_value=true))
"""

# ╔═╡ b4ae58e8-6725-11eb-1414-05378b8960a4
begin
	random = generate_data(size);
	
	or_tools = instance_OR_tools();

	instances = [or_tools, random];
	
	slice_size = @bind index_instance Select(["1" => "OR Tools", 
											  "2" => "Random",
											 ])
	
end;

# ╔═╡ bc9cbe8c-6725-11eb-0f0b-e1a0ba7acb02
md"""
Choose between a randomly generated dataset or one provided by OR Tools.

$(slice_size)
"""

# ╔═╡ d81f84b4-6725-11eb-3b99-09d3f7fd51cc
data = instances[parse(Int64, index_instance)];

# ╔═╡ 7413ed26-6729-11eb-09ce-f558ae459b6f
md"""
### Parameters configuration

"""

# ╔═╡ 5e600f58-6735-11eb-3863-f7fd4a62cf35
md"""
Then, choose the parameter configuration.

##### 1. Time limit

- CPU threshold time.

##### 2. Optimality gap tolerance

- Measured by the difference between the best lower and upper bounds.
- Two types of difference: **absolute and relative**

$GAP_{abs} = B_{upper} - B_{lower}$

$GAP_{rel} = \frac{B_{upper} - B_{lower}}{B_{upper} + 10^{-10}}$

- The $10^{-10}$ constant is to avoid dividing by zero.

##### 3. Warmstart

- We might provide initial solutions to significantly reduce computational times
- This initial solutions can be produced by a heuristic algorithm, or generated in previous iterations.

"""

# ╔═╡ 7287f2d6-6729-11eb-0871-33c510dafd9e
begin
	time_limit = 100.0 # seconds
	absolute_gap = 0.1
	relative_gap = 0.1
	stdout = false
end;

# ╔═╡ 65f25ca8-67c6-11eb-093a-7f110d8168f2
md"""
Note that if a new parameter or engine is added, it should be considered in the mapped vectors.The parameter name depends on each solver.
"""

# ╔═╡ b76fdaac-67c5-11eb-12d1-81ef3b411ad6
begin
	log_level = stdout == true ? 1 : 0
	
	gurobi_params = "TimeLimit" => time_limit, "MIPGapAbs" => absolute_gap, 
					"MIPGap" => relative_gap, "OutputFlag" => log_level
	
	clp_params = "MaximumSeconds" => time_limit, "LogLevel" => log_level
	
	cbc_params = "seconds" => time_limit, "logLevel" => log_level, 
				  "allowableGap" => absolute_gap, "ratioGap" => relative_gap
	
	glpk_params = "tm_lim" => time_limit, "msg_lev" => log_level
	
	cplex_params = "CPX_PARAM_TILIM" => time_limit, 												   "CPX_PARAM_EPAGAP" => absolute_gap,
				   "CPX_PARAM_EPGAP" => relative_gap,
				   "CPX_PARAM_MIPDISPLAY" => log_level

	all_params = Dict("Gurobi" => gurobi_params,
					  "Clp" => clp_params,
					  "Cbc" => cbc_params,
					  "GLPK" => glpk_params,
					  "CPLEX" => cplex_params)
end;

# ╔═╡ 01680330-67c8-11eb-1a37-e34d1c3fb37c
function plot_vrp(sol_routes, sol_obj, durations)
    # nodes = 

    # n_nodes = length(nodes)
    # N = 1:n_nodes

    fig = figure() 
    for i in N 
        plot(nodes[i].cx, nodes[i].cy, marker=".", color="r", label=i)
        annotate(i, (nodes[i].cx, nodes[i].cy), textcoords="offset points", xytext=(2,2), ha="center")
    end
	
    for r in sol_routes
        for n in 1:length(r)-1
            i, j = r[n], r[n + 1]
            dx = nodes[j].cx - nodes[i].cx
            dy = nodes[j].cy - nodes[i].cy
            arrow(nodes[i].cx, nodes[i].cy, dx, dy, 
                color="b", head_width=0.5, length_includes_head=true)
        end
    end
	
    sol_obj = round(sol_obj * 100) / 100
    duration = round(duration * 100) /100
	
    title("Routes")
    savefig("Routes.svg", dpi=1000)
    close(fig)
    # n = length(cities)
    # inc(a) = a == n ? one(a) : a + 1
    # lineplot([cities[inc(j)][1] for j = 0:n], [cities[inc(j)][2] for j = 0:n])
end

# ╔═╡ abac3bc8-67ea-11eb-044e-c1d760fceb77
md"""
Execute a single run: $(@bind run_single PlutoUI.CheckBox(false))
"""

# ╔═╡ 753a39a6-6726-11eb-0444-afebdae52f2d
begin
	engine = Gurobi
	if run_single
		solution = optimal_routing(data, engine, all_params[string(engine)]...)
		# plot_vrp!()
	end
end

# ╔═╡ eba9c380-67c7-11eb-03a8-0b1b0228308c
md"""
Remember that:

- Node $(length(data.service_times) + 1) is the store or **source**.
- Node $(length(data.service_times) + 2) is final depot or **sink**.
"""

# ╔═╡ 2aab45fe-6758-11eb-2687-ad29d7f3a7a2
md"""
A **sample** consists in ~200 independent runs of `optimal_routing` with `benchmark` macro.

Execute a sample: $(@bind run_single_sample PlutoUI.CheckBox(false))
"""

# ╔═╡ c556da4e-67c7-11eb-2205-9585ae72e01c
begin
	if run_single_sample
		@benchmark optimal_routing(data, engine, all_params[string(engine)]...)
	end
end

# ╔═╡ 7e30b296-6724-11eb-21a0-b9ab1a61e0e5
md"""
### Run experiment $(@bind run_experiment PlutoUI.CheckBox(false))
"""

# ╔═╡ b97ede18-6709-11eb-2133-8b91accb17cd
begin
	engines = [Gurobi, GLPK] # Cbc, CPLEX, Clp (just LP)
	if run_experiment
		experiment = benchmark_engines(data, engines, all_params)		
		add_experiment!(engine_metrics, experiment, engines)
		experiment
	end
end

# ╔═╡ 9fef3a2a-6735-11eb-3854-7d09363e5865
md"""
### Determine weights

"""

# ╔═╡ 11ca652a-676d-11eb-1f19-a9360701370f
tech_weights = Dict("Allocs" => 4, "Memory" => 4, "Time-minimum" => 7, 
				    "Time-mean" => 6, "Time-median" => 5, "Time-maximum" => 6);

# ╔═╡ cb698a20-6771-11eb-2bbc-d5943bdb0319
business_weights = Dict("Business" => 0);

# ╔═╡ a974757e-6771-11eb-1db8-d5340f4224ba
"""
prices extracted from: https://ampl.com/products/standard-price-list/
- NOTA que el precio es dms grande en magnitud a lo otro...

"""
function business_metrics(engine, business_weights)
	
	people = 2

	if engine == "Heuristic"
		impl_cost = 3 * people * 2500 # months * people * mean salary (USD)
		return (impl_cost) * business_weights["Business"]

	elseif engine in ["Cbc", "Clp", "GLPK"] # open source solvers
		impl_cost = 2 * people * 2500 # months * people * mean salary (USD)
		return (impl_cost) * business_weights["Business"]

	elseif engine in ["Gurobi", "CPLEX"] # commercial solvers
		impl_cost = 2 * people * 2500 # months * people * mean salary (USD)
		licence_cost = 20000 * 2 # mean licence cost (USD) during impl. time
		return (impl_cost + licence_cost) * business_weights["Business"]

	end
	
end

# ╔═╡ b39c5e0c-6764-11eb-02d8-05d2f6b82ba1
function rank_engines(metrics::Dict{String,Dict{Any,Any}}, 												  tech_weights::Dict{String,Int64}, 												  business_weights::Dict{String, Int64})

	rank_values = Dict()
	
	for engine in collect(keys(metrics))

		tech_cum_sum = 0
		
		for measure in metrics[engine]

			metric_name = measure[1]
			
			if metric_name in collect(keys(tech_weights))
				tech_cum_sum += mean(measure[2]) * tech_weights[metric_name]
			else
				error("missing weight for key " * metric_name)
			end
			
		end
		
		business_cum_sum = business_metrics(engine, business_weights)
		
		rank_values[engine] = tech_cum_sum + business_cum_sum
		
	end
	
	return sort(collect(filter(x -> x[2] != 0, rank_values)), by=x->x[2]) 
	
end

# ╔═╡ 2d03b1fe-6769-11eb-3bff-451c6acbe6b5
function plot_benchmark!(engine_metrics, sorted_engines, tech_weights, engines)
	
	total_plots = []
	
	# Aggregated for each engine
	# -----
	solver_names = [i[1] for i in sorted_engines]
	solver_values = [i[2] for i in sorted_engines]
	aggregated = plot(solver_names, solver_values, leg=:topleft, label="Tech and business metrics", m=:o)
	xlabel!("Engine")
	ylabel!("Weighted sum")
	title!("Engine rank")
	
	push!(total_plots, aggregated)
	
	# Disaggregated for each tech metric and engine
	# -----
	
	for metric_name in collect(keys(tech_weights))
		
		p = plot()
		
		for engine in string.(engines)
			
			values = engine_metrics[engine][metric_name]	
			names = 1:length(values)
			
			plot!(names, values, leg=:topleft, label=engine, m=:o)
		end
		
		xlabel!("Experiment ID")
		ylabel!(metric_name)
		title!(metric_name * " for each engine")
		
		push!(total_plots, p)
	end
	
	return total_plots
	
end

# ╔═╡ 8ec22ec2-6778-11eb-1978-99a19f5156e8
begin
	sorted_engines = rank_engines(engine_metrics, tech_weights, business_weights)	
	boo_worth = missing
	update_metrics = false
	if length(sorted_engines) > 0
		
		boo_worth = first(sorted_engines)[1] == string(engine) ? true : false
		
		plots = plot_benchmark!(engine_metrics, sorted_engines, tech_weights, 										engines)

		export_plots = plot(plots..., layout = length(plots), linewidth = .5, 
						markersize = 2,legendfontsize = 3, xtickfontsize = 3, 								ytickfontsize = 3, xguidefontsize = 3, yguidefontsize = 3, 							titlefontsize = 6)

		savefig(export_plots, "Outputs/EnginesBenchmark.svg")	
		update_metrics = true
	end
end;

# ╔═╡ 835669a6-6633-11eb-3966-91411b7d6da1
md"""
# Takeaways
"""

# ╔═╡ fc4e0572-67ca-11eb-1ac4-6dbdfd4a4b26
md"""
## Is it worth?
"""

# ╔═╡ 465b58c8-67ca-11eb-3f0b-87ee6d3b1a12
begin
	n_experiments = 0
	if update_metrics
		if length(engine_metrics[string(engine)]) > 0
			n_experiments = length(engine_metrics[string(engine)]["Allocs"]);
		end
		engine_metrics
	end
end

# ╔═╡ fb82d2a4-67c9-11eb-3659-33c17b751fc3
md"""
Quantity of experiments: $(n_experiments)

Ranking:
"""

# ╔═╡ 1cc641f4-67c9-11eb-3d70-77a6792f7a6f
sorted_engines

# ╔═╡ 3163ca86-6707-11eb-0669-315ac2a30ee6
md"""
**Decision**: $(engine) is a good idea - $(boo_worth)

$(LocalResource("Outputs/EnginesBenchmark.svg", :width=>1000))
"""

# ╔═╡ 3b9a473c-6707-11eb-32dc-fd046fb57eb4
md"""
Hope you like it! Export the notebook to HTML as a report.

More info at: 
- sebastian.granda@rappi.com
- diego.strobl@rappi.com
"""

# ╔═╡ Cell order:
# ╟─e145c39e-624a-11eb-392e-25c51ea78395
# ╟─345a1756-624b-11eb-0e73-b99c01d7852d
# ╟─3e591e1e-624b-11eb-3a49-e1c420a8a740
# ╟─e3877146-6495-11eb-3fde-bd2e6806a7ef
# ╟─750b95a0-6407-11eb-15b5-8b4a9805b7e8
# ╟─505933d2-66f9-11eb-1cd1-d56e8c642c12
# ╠═384c8206-66f9-11eb-22bc-e3b7d3a8e05b
# ╠═3c12dd48-66f9-11eb-1d10-ad9884fa5827
# ╠═438b08d6-66f9-11eb-2e44-256ad2bf5b11
# ╠═47554dd2-66f9-11eb-1575-717f4997cec5
# ╠═4b27ed3e-66f9-11eb-011b-3d3e1427cd49
# ╟─94bbeeae-6407-11eb-2bf5-a510e938453c
# ╟─0657a1be-66ad-11eb-233d-15f3f93307e4
# ╟─5ca26cdc-66fa-11eb-0f4b-334f9ca51c80
# ╟─86141c1e-66fa-11eb-3adf-ffb61115f195
# ╟─c8db0efe-66fa-11eb-2c99-6f2089a20cf1
# ╟─d4d31cba-66fa-11eb-2a5b-593e62edf09d
# ╟─b860657c-67d6-11eb-0240-6b84b814c3a9
# ╟─33674ca8-67d5-11eb-1d83-89ed81652979
# ╟─403d1c28-67d5-11eb-3a0f-e92303f07e3f
# ╟─f7c7f9c0-6632-11eb-27bb-e7a49dde68b8
# ╟─2f455ffe-6634-11eb-36b5-d7d9c8e2decf
# ╟─cda9bb12-6671-11eb-0494-476f9966741d
# ╟─804bcac8-6706-11eb-1d8a-756a5afde359
# ╟─48699cb4-6366-11eb-175d-8728ce809fea
# ╟─97222556-6366-11eb-03b2-e78c4247ccdf
# ╟─971ff504-6633-11eb-06d4-09b91dac0ac6
# ╟─8ee45d66-67ce-11eb-0309-6b5c68f6f09c
# ╟─ab90d342-67ce-11eb-338c-a155fa468e83
# ╟─9d7b483c-67ce-11eb-359c-a38abb43f688
# ╟─b9a63d6e-63f8-11eb-111f-afb7f1ae88c9
# ╟─1972e386-6733-11eb-2e45-f50f5240c102
# ╟─7743b4e8-6735-11eb-3ee7-117ee0d6b494
# ╟─7bd83a4c-67d8-11eb-0864-574fbb0bcc35
# ╟─7f46fae2-6735-11eb-3042-4360a215e0b1
# ╟─89b36da2-67d8-11eb-1bab-5f474524e5f5
# ╠═353c2780-6733-11eb-3bd0-25bf2e882dd8
# ╟─48f995ce-6633-11eb-0556-bf3da4f6392b
# ╟─c44cb24c-6726-11eb-30fe-bf4ab0bf02b4
# ╟─215eb286-67c5-11eb-1a2e-fba123258c06
# ╟─498502dc-6726-11eb-23f5-4151f147db6f
# ╟─1c8a6d54-6725-11eb-0025-85cacfc3b2ab
# ╟─b4ae58e8-6725-11eb-1414-05378b8960a4
# ╟─bc9cbe8c-6725-11eb-0f0b-e1a0ba7acb02
# ╟─d81f84b4-6725-11eb-3b99-09d3f7fd51cc
# ╟─7413ed26-6729-11eb-09ce-f558ae459b6f
# ╟─5e600f58-6735-11eb-3863-f7fd4a62cf35
# ╠═7287f2d6-6729-11eb-0871-33c510dafd9e
# ╟─65f25ca8-67c6-11eb-093a-7f110d8168f2
# ╟─b76fdaac-67c5-11eb-12d1-81ef3b411ad6
# ╟─01680330-67c8-11eb-1a37-e34d1c3fb37c
# ╟─abac3bc8-67ea-11eb-044e-c1d760fceb77
# ╠═753a39a6-6726-11eb-0444-afebdae52f2d
# ╟─eba9c380-67c7-11eb-03a8-0b1b0228308c
# ╟─2aab45fe-6758-11eb-2687-ad29d7f3a7a2
# ╠═c556da4e-67c7-11eb-2205-9585ae72e01c
# ╟─7e30b296-6724-11eb-21a0-b9ab1a61e0e5
# ╠═b97ede18-6709-11eb-2133-8b91accb17cd
# ╟─9fef3a2a-6735-11eb-3854-7d09363e5865
# ╠═11ca652a-676d-11eb-1f19-a9360701370f
# ╠═cb698a20-6771-11eb-2bbc-d5943bdb0319
# ╟─b39c5e0c-6764-11eb-02d8-05d2f6b82ba1
# ╟─a974757e-6771-11eb-1db8-d5340f4224ba
# ╟─2d03b1fe-6769-11eb-3bff-451c6acbe6b5
# ╟─8ec22ec2-6778-11eb-1978-99a19f5156e8
# ╟─835669a6-6633-11eb-3966-91411b7d6da1
# ╟─fc4e0572-67ca-11eb-1ac4-6dbdfd4a4b26
# ╟─465b58c8-67ca-11eb-3f0b-87ee6d3b1a12
# ╟─fb82d2a4-67c9-11eb-3659-33c17b751fc3
# ╟─1cc641f4-67c9-11eb-3d70-77a6792f7a6f
# ╟─3163ca86-6707-11eb-0669-315ac2a30ee6
# ╟─3b9a473c-6707-11eb-32dc-fd046fb57eb4
