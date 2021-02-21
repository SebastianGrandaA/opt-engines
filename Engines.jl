### A Pluto.jl notebook ###
# v0.12.21

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
	using PlutoUI, Dates, Gurobi, Cbc, Clp, GLPK, JuMP, Random, BenchmarkTools, DataStructures, LinearAlgebra, Plots, LightXML, PyCall
	project_name = "Optimization engines" 
	date = Dates.Date(Dates.now())
	company = "Rappi"
	objective = "to determine whether is more convenient to acquire a commercial solver license or develop in-house optimization algorithms."
	git = md"Find in [Github](https://github.com/SebastianGrandaA/opt-engines)"
end;

# ╔═╡ 345a1756-624b-11eb-0e73-b99c01d7852d
md"""
# $project_name @ $company

The 
**objective** is $objective 

$(git)

Last time modified: $date

"""

# ╔═╡ 3e591e1e-624b-11eb-3a49-e1c420a8a740
PlutoUI.TableOfContents()

# ╔═╡ e3877146-6495-11eb-3fde-bd2e6806a7ef
md"""
# Implementation

The **VRP** has been selected for this analysis because it represents one of the hardest operational processes at Rappi.
"""

# ╔═╡ 750b95a0-6407-11eb-15b5-8b4a9805b7e8
md"""
## Capacitated VRP with Time Windows

For this first deliverable, an adaptation of the MTZ formulation has been developed. 

CAMBIAR
"""

# ╔═╡ 07b52084-6989-11eb-3019-5776e45a0a1b
begin
	mutable struct Node
		ID   ::Int64
		lat  ::Float64
		long ::Float64
	end
	
	mutable struct RT
		ID              ::Int64
		capacity        ::Int64
		max_travel_time ::Float64
	end

	mutable struct Order
		ID           ::Int64
		delivery     ::Node
		early_time   ::Float64
		late_time    ::Float64
		items_qty    ::Int64
		service_time ::Int64
	end

	mutable struct Instance
		travel_times      ::Matrix{Float64}
		service_times     ::Vector{Int64}
		early_times       ::Vector{Int64}
		late_times        ::Vector{Int64}
		items_orders      ::Vector{Int64}
		RT_capacity       ::Int64
		RT_max_route_time ::Float64
	end

	mutable struct Solution
		routes        ::Vector{Any}
		arrival_times ::Vector{Float64}
		units         ::Vector{Float64}
		objective     ::Any
		status        ::Any
	end
end;

# ╔═╡ 8603086e-7129-11eb-1bf1-37093e6afb28
md"""
There are multiple formulations. 

**We are given**

- Complete graph $G(V, A)$.

- Time $\forall (i,j) \in A: t_{ij}$.

- Set of customers $D = V \setminus \left \{ 0, n+1 \right \}$.

- Demand $\forall i \in D: d_{i}$.

- Service time $\forall i \in D: s_{i}$.

- Time windows $\forall i \in D: [e_{i}, l_{i}]$.

- Set of homogeneous couriers $RT$.

- Capacity $\forall i \in RT: C_{i}$.

- Maximum routing time $\forall i \in RT: T_{i}$.

**Define**

- Binary $X_{ijk}$: $1$ if arc $(i,j) \in A$ is selected for courier $k \in RT$; $0$ otherwise.

- Integer $T_{ik}$: arrival time of courier $k \in RT$ at node $i \in V$.


**Objectives**

We can implement different objectives:

1. For efficient routes, **minimize the total travel times**

$MIN \sum_{k \in RT} \sum_{(i,j) \in A} t_{ij} X_{ijk}$

2. For fast routes, **minimize the makespan (length of the longest route)**

$MIN$

3. For critical orders (s.a. black stores), **minimize the total arrival times**

$MIN$

4. For longer routes, **minimize the quantity of RTs** (or routes).

$MIN$

Then, apply hierarchically [1], trade them off [2] or according to the state of the operation:

1. First minimize the number of RTs to fulfill the orders and then minimize the total travel times using that ammount of RTs.

2. Minimize the total travel time with soft time windows and add a penalty for each minute late arriving at customer' location.

-----

**Constraints**

1. Each customer is assigned to exactly one route.

$C$

2. Source to sink path in $G$ for each RT $k$.

$C$

$C$

$C$

3. Time windows 

$C | BIGM_{ij} = max( b_{i} + s_{i} + t_{ij} - a_{j} ; 0)$ (la restriccion linealizada)

$C$

4. RT capacity

$C$



----

**Apprach**


**Facts**


- The linear relaxation of this three-index model provides very weak lower bounds.

- To produce better lower bounds, reformulate as a set partitioning model: Elementary Shortest Path Problem with Resource Constraints (ESPPRC). 

- The ESPPRC generates a set of all feasible routes for the VRPTW and decides if a route visits a certain customer considering two constraints: load $\leq$ capacity and arrival time $\in$ time window.

- However, ESPPRC model contains a huge number of variables (one per feasible route).

- To handle this huge number of variables, implement a column generation algorithm.


- Furthermore, the two-index model contains an exponential number of time windows and subtour elimination constraints. Solution: apply dinamically.



**Idea**

Solve the CVRPTW with a Branch-and-Price algorithm [1, 2] in which:

- The ESPPRC arises as the subproblem of finding feasible routes (columns) with negative reduced cost that are iteratively added to the restricted master problem.

- Solve this ESPPRC with pulse and labeling algorithms [3] and the monodirectional dynamic programming method [4].

- Linear relaxations of this are solved by column generation.



Limit the solver to one second in the ESPPRC because it is the bottle neck!
add all constraints in labelling?



**References**

[1] Toth & Vigo Chapter 5: https://epubs.siam.org/doi/10.1137/1.9781611973594.ch5
[2] https://doi.org/10.1287/trsc.2014.0582
[3] https://onlinelibrary.wiley.com/doi/abs/10.1002/net.20033
[4] Feillet, D., Dejax, P., Gendreau, M., Gueguen, C., 2004. An exact algorithm for the elementary shortest path problem with resource constraints: Application to some vehicle routing problems. 


Then, of course, we can dynamically solve the VRP with reoptimization approaches (see chapter 11.3.1.2 of [1])

To sum up:

Formulating the problem will give us the flexibility to directly switch between different objectives depending on the state of the operations and will reduce times of filters such as prospect validations by adding model constraints

"""

# ╔═╡ 94bbeeae-6407-11eb-2bf5-a510e938453c
md"""
[REEMPLAZADO POR LO DE ARRIBA]


**Assumptions**

- There are enought RTs to fulfill all orders: $Q_{orders} = Q_{RTs}$.

- Homogenous RTs.

- Operational time: travel time + service time.

- Symmetric case: time to $t_{i->j} = t_{j->i}$.


**We are given:**

- Set of $DL$ delivery locations or customer places.

- Set of $S$ source [1] and sink [2] locations for a given route.

- Set of $LC = \left \{ 1, ..., DL + S \right \}$ locations or total nodes in the graph.

- Set of $TT_{i, j}$ representing the travel time between node $i \in LC$ and $j \in LC$ calculated by Rappi maps.

- Set of $RT = \left \{ 1, ..., rt \right \}$ couriers available.

- Set of $d_{i}$ demands at location $i \in LC$ or the **quantity of items in the order.**

- Set of $(E, L)_{i}$ time window pairs at location $i \in LC$: early and late times respectively.

- Parameter $c_{i}$ as the capacity of courier $i \in RT$.

- Parameter $MTT_{i}$ as the maximum travelling time of courier $i \in RT$.

- Parameter $ST_{i}$ as the service time at location $i \in LC$.

- Big $BT = MTT_{i} \times 2$ for $i \in RT$.

- Big $BQ = c_{i} + max( d_{j} )$ for $i \in RT$ and $j \in LC$.


**Decision variables**

- Binary $X_{i, j}$ indicating whether arc $(i, j) \forall i \in LC, j \in LC$ is selected or not.

- Integer $ARR_{i}$ indicating the time of arrival at node $i \in LC$.

- Integer $UN_{i}$ indicating the units or volume carried by courier $i \in RT$.


**Objective**

The objective is to minimize the total route time. This time is calculated by the sum


$MIN \sum_{i = 1}^{LC}\sum_{j = 1}^{LC} X_{i, j} \times ( ST_{i} + TT_{i, j} )$


**Contraints**

$\forall i \in DL: \sum_{j}^{DL} X_{i, j} = 1$


$\forall i \in DL, j \in DL: UN_{j} \geq UN_{i} + d_{j} - BQ \times (1 - X_{i, j})$


$\forall i \in DL: 0 \leq UN_{i} \leq c_{i}$


$\sum_{j = 1}^{LC} X_{S[1], j} \leq RT$


$\forall j \in DL: \sum_{i=1}^{LC} X_{i, j} - \sum_{i=1}^{LC} X_{j, i} = 0$


$\sum_{i=1}^{LC} X_{i, S[2]} == \sum_{j=1}^{LC} X_{S[1], j}$


$\forall i \in LC, j \in LC: ARR_{j} \geq ARR_{i} + TT_{i, j} - BT \times (1 - X_{i, j})$

$\forall i \in LC: E_{i} \leq ARR_{i} \leq L_{i}$


$\forall i \in LC: X_{i, S[1]} = 0$


$\forall i \in LC: X_{S[2], i} = 0$


$\forall i \in LC: X_{i, i} = 0$

"""

# ╔═╡ 0657a1be-66ad-11eb-233d-15f3f93307e4
md"""
## Data generation

Three main data sources:

- Solomon, M., Algorithms for the vehicle routing and scheduling problems with time window constraints. Operations Research, 1987.

- OR-Tools.

- Randomly generated.


The computational complexity is determined by:
- Relative width of time windows: wide increase the feasible solution space.
- Distribution of locations: Clustered (C), Random (R) and mixed (RC)

**Note**: last two nodes should be the store (`n-1`) and the sink (`n`).
"""

# ╔═╡ 1a078ece-698a-11eb-0e77-edbb196ffcd6
distance(n1::Node, n2::Node) = sqrt((n1.lat - n2.lat)^2 + (n1.long - n2.long)^2);

# ╔═╡ 1d364b76-698a-11eb-2133-0f8f2c2259aa
function calculate_time(nodes::Vector{Node}; digits = 2)
	
    n_nodes = length(nodes)
	store = n_nodes - 1 # start or source
	dummy = n_nodes # end or sink
    time = zeros(n_nodes, n_nodes)
	
    for i in 1:n_nodes-1
        for j in i+1:n_nodes
            d = distance(nodes[i], nodes[j])
            time[i,j] = d
            time[j,i] = d
        end
    end
	
    time[store, dummy] = 0
    time[dummy, store] = 0
	
    return time
	
end;

# ╔═╡ f0ce7818-6989-11eb-0aba-4165a1700a27
"""RBN => Random by now"""
function generate_data(size::Int64)::Instance
	
	# Time windows
	# -----
	
	early_times = rand(1:50, size) # minutes until order gets ready + minutes from 										 store to customer's location. RBN.
	
	late_times = [x+rand(5:30) for x in early_times] # early + maximum time order 														   gets "cold". RBN.

	# Demand and capacity
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
	
end;

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
	travel_times = calculate_time(nodes)
	service_times = zeros(size)
	items_orders = [1, 1, 2, 4, 2, 4, 8, 8, 1, 2, 1, 2, 4, 4, 8, 8]
	capacity = 15
	max_travel_time = 35  
	
	return Instance(travel_times, service_times, early_times, late_times, 								items_orders, capacity, max_travel_time
					)
	
end;

# ╔═╡ 74925f1e-68c6-11eb-3169-cf8db4952332
to_int(x::String) = parse(Int, x);

# ╔═╡ 0f57e64a-68c7-11eb-1ea8-fbb1abfe9b6d
to_float(x::String) = parse(Float64, x);

# ╔═╡ d537230c-68be-11eb-2bd8-b9c93c0278f5
function solomon_dataset(path_name)::Instance
	
	xlm = parse_file(path_name)
	xroot = root(xlm)
	
	# Locations
	# -----
	xnodes = xroot["network"][1]["nodes"][1]["node"]
	nodes = Node[]
	
	for i in xnodes
		id = attribute(i, "id") |> to_int
		cx = find_element(i, "cx") |> content |> to_float
		cy = find_element(i, "cy") |> content |> to_float
		push!(nodes, Node(id, cx, cy))
	end
	
	source_sink = repeat([ nodes[1] ], 2)
	deleteat!(nodes, 1)
	nodes = vcat(nodes, source_sink) # add source and sink at the end
	
	travel_times = calculate_time(nodes)
	
	# RTs
	# -----
	xrt = xroot["fleet"][1]["vehicle_profile"][1]
	capacity = convert(Int64, find_element(xrt, "capacity") |> content |> to_float)
	max_travel_time = find_element(xrt, "max_travel_time") |> content |> to_float
	
	
	# Demand
	# -----
	service_times, early_times, late_times, items_orders = [], [], [], []
	xdem = xroot["requests"][1]["request"]
	
	for i in xdem
		push!(service_times, find_element(i, "service_time") |> content |> to_float)
		
		time_window = find_element(i, "tw")
		push!(early_times, find_element(time_window, "start") |> content |> to_int)
		push!(late_times, find_element(time_window, "end") |> content |> to_int)
		
		push!(items_orders, find_element(i, "quantity") |> content |> to_float)
    end
	
	return Instance(travel_times, service_times, early_times, late_times, 								items_orders, capacity, max_travel_time)
end;

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

As-is approach.
"""

# ╔═╡ 1316a098-7455-11eb-357a-992ca346cf65
ortools = pyimport("ortools.constraint_solver.pywrapcp");

# ╔═╡ 331f7c52-7455-11eb-168d-a3845d1a211f


# ╔═╡ 72c68b7a-7455-11eb-1e52-d7752324d087
manager = ortools.RoutingIndexManager(18, 20, 0)

# ╔═╡ 403d1c28-67d5-11eb-3a0f-e92303f07e3f
function solve_or_tools(instance::Instance)::Solution
	
	
end;

# ╔═╡ 0eb32b2e-7456-11eb-2e9d-e195fe724784


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

# ╔═╡ 1bc11276-72d6-11eb-16b6-2b2f66b3d3bd
"""
Returns the routes for a given solution. Routes that start at the source
"""
function get_routes(x_val, source::Int64, sink::Int64)
	
	first_nodes = @views findall(x -> x > 0.9, x_val[source, :])
	final_routes = []
	
	for i in first_nodes

		next = @views findfirst(x -> x > 0.9, x_val[i, :])
		
		if next isa Int64
		
			route = [source, i, next]

			while next != sink && next isa Int64

				current = next

				next = @views findfirst(x -> x > 0.9, x_val[current, :])

				if next isa Int64
					push!(route, next)
				end

			end

		else  
			route = [source, i]
		end
		

		push!(final_routes, route)

	end
	
	return final_routes
end;

# ╔═╡ 30268d08-72e4-11eb-2998-e5998b88afab
function find_cycle(perm_matrix, starting_ind = 1)
	
    cycle = [starting_ind]
    prev_ind = ind = starting_ind

    while true

		next_ind = findfirst(>(0.5), @views(perm_matrix[ind, 1:prev_ind-1]))

        if isnothing(next_ind)
        
			next_ind = findfirst(>(0.5), @views(perm_matrix[ind, prev_ind+1:end])) +
                       prev_ind

        end		
		
        next_ind == starting_ind && break
        push!(cycle, next_ind)
        prev_ind, ind = ind, next_ind

    end
	
    cycle
end;

# ╔═╡ b74c997a-72e4-11eb-3510-1583a16918ff
function get_cycles(X, source, sink)
	
	perm_matrix = JuMP.value.(X).data
    N = Base.size(perm_matrix, 1)
    remaining_inds = Set(1:N)
	
	# Botar depot
	cycles = get_routes(perm_matrix, source, sink)
	setdiff!(remaining_inds, unique(vcat(cycles...)))
	
    while length(remaining_inds) > 0
		
        cycle = find_cycle(perm_matrix, first(remaining_inds))
		
        push!(cycles, cycle)
		
        setdiff!(remaining_inds, cycle)
		
	end
	
   	return cycles
end;

# ╔═╡ a4aa3352-72d0-11eb-38e6-ad004408b054
has_subtour(route, source, sink) = route[1] != source && route[end] != sink;

# ╔═╡ fc91fb98-72cd-11eb-1246-7b72dc2044cc
function violate_subtour_capacity(route, source, sink, capacity, cum_demand)::Bool
	return has_subtour(route, source, sink) || cum_demand > capacity
end;

# ╔═╡ ac42497c-72d2-11eb-3e7b-794bd578e5ac
min_couriers(demand, capacity) = demand / capacity;

# ╔═╡ bcb00170-72ce-11eb-04eb-61c2f3a5d5a6
function violate_time_windows(route, late_times, times)::Bool
	
	cum_time = 0
	
	for i in route
		
		if i <= length(late_times) - 1 # only customers
			
			cum_time += times[i, i+1] # + service_times[i]
			
			if cum_time > late_times[i]
				
				return true
				
			end
		end
	end
	
	return false
end;

# ╔═╡ b4647d92-72d6-11eb-1ff9-77c3db9ae6f2
function remove_invalid(model, X, source, sink, capacity, demand, late_times, times)

    routes = get_routes(X, source, sink)
	
	if length(routes) > 1
		
		for route in routes

			zipped = zip(route[1:end-1], route[2:end])

			if violate_subtour_capacity(route, capacity, demand) # 5.17

				@info("Invalidating with soubtour or capacity constraints")

				thr = min_couriers(demand, capacity)

				@constraint(model, sum(X[i, j] for (i,j) in zipped) >= thr)

			elseif violate_time_windows(route, late_times, times) # 5.18

				@info("Invalidating with TW constraint")

				thr = length(route) - 1 # set of arcs in route 

				@constraint(model, sum(X[i, j] for (i,j) in zipped) <= thr)

			end

		end
		
		return true
		
	end
	
	return false
end;

# ╔═╡ 8b4bc38e-7196-11eb-315c-0922c360b5f6
"""
"""
function callback_funct(model, X, source, sink, capacity, demand, late_times, times)
	
	return function vrp_callback(cb_data)
		
        X_vals = callback_value.(Ref(cb_data), X)
		
        any(x -> !(x ≈ round(Int, x)), X_vals) && return # only integer solutions	

		routes = get_cycles(X_vals, source, sink)
		
		if length(routes) > 1
						
			for route in routes
					
				zipped = zip(route[1:end-1], route[2:end])
				
				thr = - 10
				
				cum_demand = sum([demand[i] for i in route if i <= length(demand)])
				
				if violate_subtour_capacity(route, source, sink, capacity, cum_demand) # 5.17
					
					@warn("Subtour or capacity violated by route: $route")
					
					thr = min_couriers(cum_demand, capacity)
					
					con = @build_constraint(sum(X[i, j] for (i,j) in zipped) >= thr)
					
				elseif violate_time_windows(route, late_times, times) # 5.18
					
					@warn("TW violated by route: $route")
					
					thr = length(route) - 1 # set of arcs in route 
					
					con = @build_constraint(sum(X[i, j] for (i,j) in zipped) <= thr)
				
				end
				
				if thr > 0
					
					@info("Submitting lazy constraint")
					
					MOI.submit(model, MOI.LazyConstraint(cb_data), con)
					
				end

			end
			
		else
			
			@info("Found full routes with no subtours.")
            
			return nothing
			
        end
    end
end;

# ╔═╡ 6b72d8f2-72d5-11eb-2b47-f9b3330f0188
md"""
el if i <= length() es redundante porque igual late_times tiene customers + 2

**Prueba solo con el callback de TW y con los dos!!**
"""

# ╔═╡ e0739df2-6f8e-11eb-0341-e5e67b838f21
"""

**Idea**

- Cutting plane method embedded into a Branch-and-Bound procedure. 
- Cutting planes are added in each node of the search three to thighten the LP relaxation as much as possible.
- Therefore, **time windows and capacity** constraints are added dynamically with callbacks (lazy constraints).


**Constraints**

"""
function lazy_model(instance::Instance, solver::Module, solver_params::Pair...)
	
    # Parameters
    # -----

	V_ = Base.size(instance.travel_times, 1) # nodes
	V = 1:V_
    
	D_ = length(instance.service_times) # customers
    D = 1:D_

    source = V_ - 1
    sink = V_
	
    RT_ = D_ # RTs
    RT = 1:RT_

    d,C,T = instance.items_orders, instance.RT_capacity, instance.RT_max_route_time

	early_times = [instance.early_times; 0; 0]
	late_times = [instance.late_times; 0; T]
    time = instance.travel_times

    # Formulation
    # -----

	model = Model(optimizer_with_attributes(solver.Optimizer, solver_params...))

	@variable(model, X[V, V] >= 0, Int)

	@objective(model, Min, sum(time[i,j] * X[i,j] for i in V, j in V)) 
	
	@constraint(model, [i in D], sum(X[i,j] for j in V) == 1)
	@constraint(model, [j in D], sum(X[i,j] for i in V) == 1)
	
	# separar
	# @constraint(model, [k in D], sum(X[i,k] for i in V) - sum(X[k,j] for j in V)==0)
	
	@constraint(model, X[source, sink] <= RT_ - 1)
	@constraint(model, [i in D, j in D], X[i,j] + X[j,i] <= 1)
	
	@constraint(model, sum(X[source, j] for j in V) == RT_)
	
	@constraint(model, sum(X[i, sink] for i in V) == sum(X[source, j] for j in V))
	
	@constraint(model, [i in V], X[i, source] == 0)
	@constraint(model, [i in V], X[sink, i] == 0) 
	@constraint(model, [i in V], X[i,i] == 0)
	
	
    # Optimize!
    # -----
	
	callback = callback_funct(model, X, source, sink, C, d, late_times, time)	
	
	MOI.set(model, MOI.LazyConstraintCallback(), callback)
	
	boo_remove = true
	
    # while boo_remove

	optimize!(model)
	
	# boo_remove = remove_invalid(model, X, source, sink, C, d, late_times, time)
		
    # end
	
	status = JuMP.termination_status(model)
	
	final_route = get_routes(JuMP.value.(X).data, source, sink)
	
	return Solution(final_route, [], [], objective_value(model), status)
	
end;

# ╔═╡ 48c406c2-71c6-11eb-0d89-c90f3642a39d
"""
Validations
-----

[1] Arrival time should happen in time window
[2] Total travel times should be less than max travel time
[3] Demand carried should be less than capacity
[4] Each customer is visited once.
[5] Objective value of solution is sum of route times
[6] All customers fulfilled.


"""
function validate_solution(solution, data)
	routes = solution.routes
	arrivals = solution.arrival_times
	units = solution.units
	objective = solution.objective
	
	times = data.travel_times
	early = data.early_times
	late = data.late_times
	demand = data.items_orders
	max_capacity = data.RT_capacity
	max_time = data.RT_max_route_time
	
	global_time = 0

	for route in routes
		
		total_time, total_demand = 0, 0
		
		for c in 1:length(route) - 1

			total_time += times[route[c], route[c+1]]
			
			if route[c] <= length(data.service_times) # only customers
				
				@assert arrivals[route[c]] <= late[route[c]] # 1
			
			end
			
		end
		
		total_demand = units[route[end-2]]
		
		@assert total_time <= max_time # 2

		@assert total_demand <= max_capacity # 3
		
		global_time += total_time
		
	end
	
	all_nodes = vcat(routes...)

	@assert length(all_nodes) - 2*(length(routes)-1)==length(unique(all_nodes)) # 4

	@assert round(global_time) == round(objective) # 5

	@assert unique(sort(all_nodes)) == collect(1:length(arrivals)) # 6
	
end;

# ╔═╡ a0cd2ef2-71c6-11eb-1efa-15c088d4b8d0
# validate_solution(lm, deepcopy(data))

# ╔═╡ 39d05d1e-7234-11eb-0338-abc06b474347
# @benchmark lazy_model(deepcopy(data), engine, all_params[string(engine)]...)

# ╔═╡ 3c8221d6-718b-11eb-1465-171acd331599
# mtzm = mtz_model(deepcopy(data), engine, all_params[string(engine)]...)

# ╔═╡ d3ee2f86-71c9-11eb-33bf-9391c8ec534c
# validate_solution(mtzm, deepcopy(data))

# ╔═╡ 64376310-71ab-11eb-2719-eff1cb5b7c7f
# @benchmark mtz_model(deepcopy(data), engine, all_params[string(engine)]...)

# ╔═╡ cda9bb12-6671-11eb-0494-476f9966741d
"""
MTZ
-----

**Constraints**

(1) Each customer is visited only once.
(2) Flow in == flow out
(3) Availability on RT quantity.

(4) Subtour elimination and connectivity.
(5) Capacity should not be exceeded.
(6) Routes should be from source -> sink.

(7) Schedule feasibility: arrival at node j should be after visiting node i
(8) Time windows.

(9-10) Source and sink nodes should be visited first and last respectively
(11) Avoid cycling

"""
function mtz_model(instance::Instance, solver::Module, solver_params::Pair...)

    # Parameters
    # -----

	V_ = Base.size(instance.travel_times, 1) # nodes
	V = 1:V_
    
	D_ = length(instance.service_times) # customers
    D = 1:D_

    source = V_ - 1
    sink = V_
	
    RT_ = D_ # RTs
    RT = 1:RT_

    d,C,T = instance.items_orders, instance.RT_capacity, instance.RT_max_route_time

	early_times = [instance.early_times; 0; 0]
	late_times = [instance.late_times; 0; T]
    time = instance.travel_times
	
	bigT = T * 2
    bigQ = C + maximum(d)
    

    # Formulation
    # -----

	model = Model(optimizer_with_attributes(solver.Optimizer, solver_params...))
	
	@variable(model, X[V, V], Bin)
	@variable(model, A[V] >= 0) # arrival time
	@variable(model, U[V] >= 0) # units carried in route

	@objective(model, Min, sum(time[i,j] * X[i,j] for i in V, j in V)) 
	
	@constraint(model, [i in D], sum(X[i,j] for j in V) == 1)
	@constraint(model, [k in D], sum(X[i,k] for i in V) - sum(X[k,j] for j in V)==0)
	@constraint(model, sum(X[source, j] for j in V) <= RT_)
	
	@constraint(model, [i in V, j in D], U[j] >= U[i] + d[j] - bigQ * (1 - X[i,j]))
	@constraint(model, [i in V], 0 <= U[i] <= C)
	@constraint(model, sum(X[i, sink] for i in V) == sum(X[source, j] for j in V))
	
	@constraint(model, [i in V, j in V], A[j] >= A[i] + time[i,j] - bigT*(1-X[i,j]))
	@constraint(model, [i in V], early_times[i] <= A[i] <= late_times[i])
	
	@constraint(model, [i in V], X[i, source] == 0)
	@constraint(model, [i in V], X[sink, i] == 0) 
	@constraint(model, [i in V], X[i,i] == 0)
	
	
    # Optimize!
    # -----

    optimize!(model)
    
    status = JuMP.termination_status(model)
	
	final_rout = get_routes(JuMP.value.(X).data, source, sink)
	final_arr = @views collect(JuMP.value.(A))
	final_un = @views collect(JuMP.value.(U))
	
	return Solution(final_rout, final_arr, final_un, objective_value(model), status)

end;

# ╔═╡ bbb41a3a-6f8e-11eb-02c9-fbd9f8088eca
function optimal_routing(instance::Instance, solver::Module, solver_params::Pair...)
	
	callback_feature = ["Gurobi", "CPLEX", "GLPK"]
	
	if string(solver) in callback_feature
		
		@warn("Using model with lazy constraints in " * string(solver))

		lazy_model(instance, solver, solver_params...)
		
	else
		
		@warn("Using model with normal constraints in " * string(solver))
		
		mtz_model(instance, solver, solver_params...)
		
	end	
end;

# ╔═╡ 97222556-6366-11eb-03b2-e78c4247ccdf
md"""
### Heuristics
"""

# ╔═╡ 971ff504-6633-11eb-06d4-09b91dac0ac6
md"""
In broad terms, a two-phase process is implemented:

1. **Construction heuristic**: generate an initial solution to start from (greedy algorithm)

2. **Local-search metaheuristic**: iteratively improve the initial solution by applying several move operators.

"""

# ╔═╡ 8ee45d66-67ce-11eb-0309-6b5c68f6f09c
function construction(data)
	
end;

# ╔═╡ ab90d342-67ce-11eb-338c-a155fa468e83
function move_operator!(data)

end;

# ╔═╡ 9d7b483c-67ce-11eb-359c-a38abb43f688
function local_search(initial_solution)
	
end;

# ╔═╡ a69f3478-6f92-11eb-38cd-edb17d49d2eb
function solve_heuristic()
	
end;

# ╔═╡ b9a63d6e-63f8-11eb-111f-afb7f1ae88c9
md"""
# Analysis
"""

# ╔═╡ 1972e386-6733-11eb-2e45-f50f5240c102
md"""
## Metrics

The best solver would be the one with the **minimum weighted sum** of the following technical and business metrics.
"""

# ╔═╡ 7743b4e8-6735-11eb-3ee7-117ee0d6b494
md"""
##### Technical

These metrics will be calculated multiple times and compared between engines at the end of the **benchmark experiments** (see next section). 
"""

# ╔═╡ 7bd83a4c-67d8-11eb-0864-574fbb0bcc35
md"""
Time

- The **minimum** is a robust estimator for the location parameter of the time distribution, and should not be considered an outlier.

- The **median** is as a robust measure of central tendency, should be relatively unaffected by outliers.

- The **mean** is a non-robust measure of central tendency, will usually be positively skewed by outliers

- The **maximum** should be considered a primarily noise-driven outlier, and can change drastically between benchmark trials.


Allocs

- Units of estimated memory allocations.

Memory

- Bytes of estimated memory usage.

Solution quality

- Objective value of output model solution.

Consistency

- Behaviour between benchmark experiments.

"""

# ╔═╡ 7f46fae2-6735-11eb-3042-4360a215e0b1
md"""
##### Business

These metrics are estimated depending on the engine' license terms and implementation.
"""

# ╔═╡ 89b36da2-67d8-11eb-1bab-5f474524e5f5
md"""

Price

- Estimated price found in website.

Long-term contract

- Years of agreement' validation.

Implementation cost

- Estimated as: $C_{impl} = months \times people \times salary$

"""

# ╔═╡ 53bdc282-68d1-11eb-3188-5d845aba50fc
function initialize_metrics()
	return Dict("Gurobi" => Dict(), "Clp" => Dict(), "Cbc" => Dict(), "GLPK" => 				Dict(), "CPLEX" => Dict(), "Heuristic" => Dict(), "OR-Tools"=>Dict())
end;

# ╔═╡ 353c2780-6733-11eb-3bd0-25bf2e882dd8
engine_metrics = initialize_metrics();

# ╔═╡ 48f995ce-6633-11eb-0556-bf3da4f6392b
md"""
## Benchmark

Now, we'll technically analyze the performance of each solver. 

The idea is to **incrementally add new experiments** and, at the end, **compare the individual and cumulative metrics between engines**.


Let's define:

- **Scenario** as the combination of a data instance and solver parameters.

- **Run** as an optimization procedure to solve the scenario with a given engine. The output is a solution for the problem.

- **Sample** as an execution of ~200 independent runs. The output is a complete set of measurements for that scenario optimized by the engine.

- **Experiment** as a set of samples for different solvers.


This way we ensure **repeatability and reproducibility** in the experiments.

The experiments were executed locally in a MacBook Pro (16-inch, 2019).
- Processor: 2.3 GHz 8-Core Intel Core i9
- Memory: 16 GB 2667 MHz DDR4
- Graphics: Intel UHD Graphics 630 1536 MB

"""

# ╔═╡ c44cb24c-6726-11eb-30fe-bf4ab0bf02b4
"""Given set of parameters, perform the benchmark for each solver"""
function benchmark_engines(data::Instance, engines::Vector{Module}, all_params::Any)

	suite, objective_values = BenchmarkGroup(["engines"]), Dict()
	
	for engine in engines
		
		engine_name = string(engine)
		
		e_params = all_params[engine_name]

		objective_values[engine_name] = optimal_routing(deepcopy(data), engine, 														e_params...).objective

		suite[engine_name] = @benchmark optimal_routing(deepcopy(data), engine, 																e_params...)

	end

	return suite, objective_values
	
end;

# ╔═╡ 215eb286-67c5-11eb-1a2e-fba123258c06
"""
Updates `engine_metrics` with the new information of the experiment.
"""
function add_experiment!(engine_metrics::Dict{String,Dict{Any,Any}}, 										 experiment::BenchmarkGroup, engines::Vector{Module}, 								 objective_values::Dict{Any, Any})

	for engine in string.(engines)
		
		allocs_estimated = experiment[engine].allocs
		memory_estimated = experiment[engine].memory
		objective_value = objective_values[engine]
		
		if length(collect(engine_metrics[engine])) == 0
			engine_metrics[engine]["Allocs"] = [allocs_estimated]
			engine_metrics[engine]["Memory"] = [memory_estimated]
			engine_metrics[engine]["OV"] = [objective_value]
			
		else # not first time
			push!(engine_metrics[engine]["Allocs"], allocs_estimated)
			push!(engine_metrics[engine]["Memory"], memory_estimated)
			push!(engine_metrics[engine]["OV"], objective_value)
		end
		
		for func in [minimum, median, mean, maximum]
			
			raw_measure = func(experiment[engine]) # miliseconds
			
			measure_t = raw_measure.time
			name_t = "Time-"*string(func)
			
			measure_gc = raw_measure.gctime
			name_gc = "GCTime-"*string(func)
						
			if name_t in collect(keys(engine_metrics[engine])) # not first time
				push!(engine_metrics[engine][name_t], measure_t) # t
				push!(engine_metrics[engine][name_gc], measure_gc) # gc
			else
				engine_metrics[engine][name_t] = [measure_t]
				engine_metrics[engine][name_gc] = [measure_gc]
			end
			
		end
		
	end
	
	return engine_metrics
	
end;

# ╔═╡ 498502dc-6726-11eb-23f5-4151f147db6f
md"""

### Scenarios

We'll run an experiment for each scenario:

|   Source  	| Qty 	| Size 	|    Params   	|
|:---------:	|:---:	|:----:	|:-----------:	|
|  OR Tools 	|  2  	|  16  	|     None    	|
|  Solomon  	|  10 	|  25  	|     None    	|
|  Solomon  	|  10 	|  50  	|  10s or .1% 	|
|  Solomon  	|  10 	|  100 	|  30s or .2% 	|
| Homberger 	|  5  	|  200 	|  50s or .5% 	|
| Homberger 	|  5  	|  400 	| 100s or .6% 	|
| Homberger 	|  5  	|  600 	| 100s or .6% 	|
| Homberger 	|  5  	|  800 	| 100s or .6% 	|
| Homberger 	|  5  	| 1000 	| 100s or .6% 	|

`Qty` is the number of different instances in the `Source`, `Size` is the quantity of nodes in the instance' network and `Params` the time and gap thresholds for the solver.
"""

# ╔═╡ 1c8a6d54-6725-11eb-0025-85cacfc3b2ab
md"""
### Instance selection

For randomly generated instances: $(@bind size Slider(1:300, default=10, show_value=true))

"""

# ╔═╡ 8b2d6f44-6992-11eb-2c77-25e14955515c
md"""
For Solomon datasets, there are many `.xml` files at the `Input` folder.
"""

# ╔═╡ b0fbb7f2-68d0-11eb-3d20-f5f015d3b4d1
path_name = joinpath("Inputs", "solomon-1987", "C101_100.xml");

# ╔═╡ a85fed86-68bf-11eb-2bca-9985e71090cb
md"""
Select instance: $(@bind index_instance Select(["2" => "OR Tools", 
								     	  		"1" => "Solomon",
										  		"3" => "Random",
											 ]))
"""

# ╔═╡ b4ae58e8-6725-11eb-1414-05378b8960a4
begin
	if index_instance == "1"
		data = solomon_dataset(path_name)
		
	elseif index_instance == "2"
		data = instance_OR_tools()
		
	elseif index_instance == "3"
		data = generate_data(size)

	else
		error("Select a dataset")
	end
	
end;

# ╔═╡ 23f824d6-7452-11eb-188a-49f801f9d60e
data

# ╔═╡ c8c0906a-68cb-11eb-0ed0-7341082f74a0
md"""
Quantity of nodes: $(Base.size(data.travel_times, 1))
"""

# ╔═╡ 82073998-6bbf-11eb-158c-455db68ee805
data

# ╔═╡ 7413ed26-6729-11eb-09ce-f558ae459b6f
md"""
### Parameters configuration

"""

# ╔═╡ 5e600f58-6735-11eb-3863-f7fd4a62cf35
md"""
Finally, choose the parameters configuration.

##### 1. Time limit

- CPU threshold time.

##### 2. Optimality gap tolerance

- Measured by the difference between the best lower and upper bounds.


- Absolute: $GAP_{abs} = B_{upper} - B_{lower}$


- Relative: $GAP_{rel} = \frac{B_{upper} - B_{lower}}{B_{upper} + 10^{-10}}$


- The $10^{-10}$ constant is to avoid dividing by zero.

##### 3. Warmstart

- We might provide initial solutions to significantly reduce computational times.


- This initial solutions can be produced by a heuristic algorithm, or generated in previous iterations.

"""

# ╔═╡ 7287f2d6-6729-11eb-0871-33c510dafd9e
begin
	time_limit = 10.0 # seconds
	absolute_gap = .01
	relative_gap = .01
	stdout = false
end;

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
end;

# ╔═╡ abac3bc8-67ea-11eb-044e-c1d760fceb77
md"""
Execute a single run: $(@bind run_single PlutoUI.CheckBox(false))
"""

# ╔═╡ af883f14-6a77-11eb-0bac-0be99ccd9413
engine = Gurobi;

# ╔═╡ 3c273b04-718b-11eb-1263-2def968b9848
lm = lazy_model(deepcopy(data), engine, all_params[string(engine)]...)

# ╔═╡ 753a39a6-6726-11eb-0444-afebdae52f2d
begin
	if run_single
		try
			solution = optimal_routing(deepcopy(data), engine, 																all_params[string(engine)]...)
			# plot_vrp!()
			
		catch err
			solution = Solution(Dict(missing=>missing), missing, err) 
		end

	end
end

# ╔═╡ eba9c380-67c7-11eb-03a8-0b1b0228308c
md"""
Recall:
- Node $(length(data.service_times) + 1) is the store or **source**.
- Node $(length(data.service_times) + 2) is final depot or **sink**.
"""

# ╔═╡ 2aab45fe-6758-11eb-2687-ad29d7f3a7a2
md"""
Execute a sample: $(@bind run_single_sample PlutoUI.CheckBox(false))
"""

# ╔═╡ c556da4e-67c7-11eb-2205-9585ae72e01c
begin
	if run_single_sample

		@benchmark optimal_routing(deepcopy(data), engine, 																		all_params[string(engine)]...)
		
	end
end

# ╔═╡ 7e30b296-6724-11eb-21a0-b9ab1a61e0e5
md"""
### Run experiment $(@bind run_experiment PlutoUI.CheckBox(false))
"""

# ╔═╡ cf83b88a-68aa-11eb-07e1-1ddc5da72910
engines = [Gurobi]; # , CPLEX, Clp (just LP)

# ╔═╡ b97ede18-6709-11eb-2133-8b91accb17cd
begin
	if run_experiment
		experiments, objective_values = benchmark_engines(deepcopy(data), engines, 																all_params)
		add_experiment!(engine_metrics, experiments, engines, objective_values)
	end
end;

# ╔═╡ 9fef3a2a-6735-11eb-3854-7d09363e5865
md"""
### Determine weights

"""

# ╔═╡ 11ca652a-676d-11eb-1f19-a9360701370f
tech_weights = Dict("OV"=>7, "Time-minimum" => 7, "Time-mean" => 6, 
					"Time-median" => 5, "Time-maximum" => 6, "Allocs" => 4, 
					"Memory" => 4, "GCTime-minimum" => 2, "GCTime-mean" => 1, 							"GCTime-median" => 1, "GCTime-maximum" => 2);

# ╔═╡ cb698a20-6771-11eb-2bbc-d5943bdb0319
business_weights = Dict("Business" => 0);

# ╔═╡ a974757e-6771-11eb-1db8-d5340f4224ba
"""
prices extracted from: https://ampl.com/products/standard-price-list/
- NOTA que el precio es dms grande en magnitud a lo otro...

"""
function business_metrics(engine, business_weights; people = 2)
	
	if engine == "Heuristic"
		impl_cost = 3 * people * 2500 # months * people * mean salary (USD)
		return (impl_cost) * business_weights["Business"]
		
	elseif engine == "OR-Tools"
		impl_cost = 0 * people * 2500 # months * people * mean salary (USD)
		return (impl_cost) * business_weights["Business"]

	elseif engine in ["Cbc", "Clp", "GLPK"] # open source solvers
		impl_cost = 2 * people * 2500 # months * people * mean salary (USD)
		return (impl_cost) * business_weights["Business"]

	elseif engine in ["Gurobi", "CPLEX"] # commercial solvers
		impl_cost = 2 * people * 2500 # months * people * mean salary (USD)
		licence_cost = 20000 * 2 # mean licence cost (USD) during impl. time
		return (impl_cost + licence_cost) * business_weights["Business"]
		
	else
		error("An engine is not beeing considered in business_metrics")

	end
	
end;

# ╔═╡ b39c5e0c-6764-11eb-02d8-05d2f6b82ba1
"""
Aggregates tech and business metrics and ranks each engine.
"""
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
	
end;

# ╔═╡ 2d03b1fe-6769-11eb-3bff-451c6acbe6b5
"""
Only plots the tech metrics and the aggreated ones because they depend on the experiment run. Business metrics are constant along experiments.
"""
function plot_benchmark!(engine_metrics, sorted_engines, tech_weights, engines)
	
	total_plots = []
	
	# Aggregated for each engine
	# -----
	solver_names = [i[1] for i in sorted_engines]
	solver_values = [i[2] for i in sorted_engines]
	aggregated = bar(solver_names, solver_values, leg=:topright, label="Cumulative")
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
		title!(metric_name)
		
		push!(total_plots, p)
	end
	
	return total_plots
	
end;

# ╔═╡ 8ec22ec2-6778-11eb-1978-99a19f5156e8
begin
	boo_worth = missing
	
	if run_experiment	
		
		sorted_engines = rank_engines(engine_metrics, tech_weights, 													  business_weights)	

		if length(sorted_engines) > 0

			boo_worth = first(sorted_engines)[1] == string(engine) ? true : false

			plots = plot_benchmark!(engine_metrics, sorted_engines, tech_weights, 										engines)

			export_plots = plot(plots..., layout = length(plots), linewidth = .5, 
							markersize = 2,legendfontsize = 3, xtickfontsize = 3, 								ytickfontsize = 3, xguidefontsize = 3, 
							yguidefontsize = 3, titlefontsize = 6)

			savefig(export_plots, "Outputs/EnginesBenchmark.svg")	

		end
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
	n_experiments = missing
	if run_experiment
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

# ╔═╡ 87e40234-68ad-11eb-30b8-df265cd112f7
run_experiment == true ? sorted_engines : missing

# ╔═╡ 3163ca86-6707-11eb-0669-315ac2a30ee6
md"""
**Decision**: $(engine) is a good idea - $(boo_worth)

$(LocalResource("./Outputs/EnginesBenchmark.svg", :width=>1000))
"""

# ╔═╡ f5e81bda-6993-11eb-11eb-5bf3ee416fc9
md"""

## Next steps

`lazy_model()`

- Replace MTZ sub-tours constraints with **lazy constraints as [callbacks](https://jump.dev/JuMP.jl/v0.21.1/callbacks/index.html#Available-solvers-1).**

- Implement Branch and Price framework.

- Ultima corroboracion con Routing.jl VSC

- Resources: Docs Diego, Toth, J-B VRP, [Gurobi TSP](https://www.gurobi.com/documentation/9.0/examples/tsp_py.html), [Gurobi VRP](https://support.gurobi.com/hc/en-us/community/posts/360057640171-VRP-model-is-infeasible), 

----

`solve_or_tools()`

- Implement OR-tools engine and benchmark

-  Resources: .py Diego and [docs](https://developers.google.com/optimization/routing/vrptw)

----

`Data generation`

- Add bigger datasets such as [Homberger](https://www.sintef.no/projectweb/top/vrptw/homberger-benchmark/) and [Russell](https://neo.lcc.uma.es/vrp/vrp-instances/capacitated-vrp-with-time-windows-instances/).

----

`solve_heuristic()`

- Develop heuristic algorithms. Checkout [Nodal.jl](https://github.com/phrb/NODAL.jl) and Engine.jl (local).

----

`Analysis`

- Ver la propuesta de meter todos los post validadores al modelo y benchmark

- Resources: [ppt](https://docs.google.com/presentation/d/13Nn6rNmHrv0Vgm8Os_clfuyPAFaRJwNMd9Boi0XuIl4/edit#slide=id.g8fd89a81d9_0_233)

----

`Others`

- Plot solutions.

- [Nextmv](https://docs.nextmv.io/latest/)

- Heterogeneous RTs: different capacity and travel times.

- New features: open and pickup/delivery

- Check python implementations: [1](https://github.com/chkwon/vrpy), [2](https://github.com/chkwon/foodora-routing-problem).

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
# ╟─07b52084-6989-11eb-3019-5776e45a0a1b
# ╟─8603086e-7129-11eb-1bf1-37093e6afb28
# ╟─94bbeeae-6407-11eb-2bf5-a510e938453c
# ╟─0657a1be-66ad-11eb-233d-15f3f93307e4
# ╟─1a078ece-698a-11eb-0e77-edbb196ffcd6
# ╟─1d364b76-698a-11eb-2133-0f8f2c2259aa
# ╟─f0ce7818-6989-11eb-0aba-4165a1700a27
# ╟─d4d31cba-66fa-11eb-2a5b-593e62edf09d
# ╟─74925f1e-68c6-11eb-3169-cf8db4952332
# ╟─0f57e64a-68c7-11eb-1ea8-fbb1abfe9b6d
# ╟─d537230c-68be-11eb-2bd8-b9c93c0278f5
# ╟─b860657c-67d6-11eb-0240-6b84b814c3a9
# ╟─33674ca8-67d5-11eb-1d83-89ed81652979
# ╠═23f824d6-7452-11eb-188a-49f801f9d60e
# ╠═1316a098-7455-11eb-357a-992ca346cf65
# ╠═331f7c52-7455-11eb-168d-a3845d1a211f
# ╠═72c68b7a-7455-11eb-1e52-d7752324d087
# ╠═403d1c28-67d5-11eb-3a0f-e92303f07e3f
# ╠═0eb32b2e-7456-11eb-2e9d-e195fe724784
# ╟─f7c7f9c0-6632-11eb-27bb-e7a49dde68b8
# ╟─2f455ffe-6634-11eb-36b5-d7d9c8e2decf
# ╟─bbb41a3a-6f8e-11eb-02c9-fbd9f8088eca
# ╠═b4647d92-72d6-11eb-1ff9-77c3db9ae6f2
# ╠═1bc11276-72d6-11eb-16b6-2b2f66b3d3bd
# ╠═30268d08-72e4-11eb-2998-e5998b88afab
# ╠═b74c997a-72e4-11eb-3510-1583a16918ff
# ╠═a4aa3352-72d0-11eb-38e6-ad004408b054
# ╠═fc91fb98-72cd-11eb-1246-7b72dc2044cc
# ╠═ac42497c-72d2-11eb-3e7b-794bd578e5ac
# ╠═bcb00170-72ce-11eb-04eb-61c2f3a5d5a6
# ╠═8b4bc38e-7196-11eb-315c-0922c360b5f6
# ╠═6b72d8f2-72d5-11eb-2b47-f9b3330f0188
# ╠═e0739df2-6f8e-11eb-0341-e5e67b838f21
# ╟─48c406c2-71c6-11eb-0d89-c90f3642a39d
# ╠═3c273b04-718b-11eb-1263-2def968b9848
# ╠═a0cd2ef2-71c6-11eb-1efa-15c088d4b8d0
# ╠═39d05d1e-7234-11eb-0338-abc06b474347
# ╠═3c8221d6-718b-11eb-1465-171acd331599
# ╠═d3ee2f86-71c9-11eb-33bf-9391c8ec534c
# ╠═64376310-71ab-11eb-2719-eff1cb5b7c7f
# ╟─cda9bb12-6671-11eb-0494-476f9966741d
# ╟─97222556-6366-11eb-03b2-e78c4247ccdf
# ╟─971ff504-6633-11eb-06d4-09b91dac0ac6
# ╟─8ee45d66-67ce-11eb-0309-6b5c68f6f09c
# ╟─ab90d342-67ce-11eb-338c-a155fa468e83
# ╟─9d7b483c-67ce-11eb-359c-a38abb43f688
# ╟─a69f3478-6f92-11eb-38cd-edb17d49d2eb
# ╟─b9a63d6e-63f8-11eb-111f-afb7f1ae88c9
# ╟─1972e386-6733-11eb-2e45-f50f5240c102
# ╟─7743b4e8-6735-11eb-3ee7-117ee0d6b494
# ╟─7bd83a4c-67d8-11eb-0864-574fbb0bcc35
# ╟─7f46fae2-6735-11eb-3042-4360a215e0b1
# ╟─89b36da2-67d8-11eb-1bab-5f474524e5f5
# ╟─53bdc282-68d1-11eb-3188-5d845aba50fc
# ╠═353c2780-6733-11eb-3bd0-25bf2e882dd8
# ╟─48f995ce-6633-11eb-0556-bf3da4f6392b
# ╟─c44cb24c-6726-11eb-30fe-bf4ab0bf02b4
# ╟─215eb286-67c5-11eb-1a2e-fba123258c06
# ╟─498502dc-6726-11eb-23f5-4151f147db6f
# ╟─1c8a6d54-6725-11eb-0025-85cacfc3b2ab
# ╟─8b2d6f44-6992-11eb-2c77-25e14955515c
# ╠═b0fbb7f2-68d0-11eb-3d20-f5f015d3b4d1
# ╟─a85fed86-68bf-11eb-2bca-9985e71090cb
# ╟─b4ae58e8-6725-11eb-1414-05378b8960a4
# ╟─c8c0906a-68cb-11eb-0ed0-7341082f74a0
# ╟─82073998-6bbf-11eb-158c-455db68ee805
# ╟─7413ed26-6729-11eb-09ce-f558ae459b6f
# ╟─5e600f58-6735-11eb-3863-f7fd4a62cf35
# ╠═7287f2d6-6729-11eb-0871-33c510dafd9e
# ╟─b76fdaac-67c5-11eb-12d1-81ef3b411ad6
# ╟─01680330-67c8-11eb-1a37-e34d1c3fb37c
# ╟─abac3bc8-67ea-11eb-044e-c1d760fceb77
# ╠═af883f14-6a77-11eb-0bac-0be99ccd9413
# ╟─753a39a6-6726-11eb-0444-afebdae52f2d
# ╟─eba9c380-67c7-11eb-03a8-0b1b0228308c
# ╟─2aab45fe-6758-11eb-2687-ad29d7f3a7a2
# ╟─c556da4e-67c7-11eb-2205-9585ae72e01c
# ╟─7e30b296-6724-11eb-21a0-b9ab1a61e0e5
# ╠═cf83b88a-68aa-11eb-07e1-1ddc5da72910
# ╟─b97ede18-6709-11eb-2133-8b91accb17cd
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
# ╟─87e40234-68ad-11eb-30b8-df265cd112f7
# ╟─3163ca86-6707-11eb-0669-315ac2a30ee6
# ╟─f5e81bda-6993-11eb-11eb-5bf3ee416fc9
# ╟─3b9a473c-6707-11eb-32dc-fd046fb57eb4
