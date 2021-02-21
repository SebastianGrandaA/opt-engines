"""Vehicles Routing Problem (VRP) with Time Windows."""
import math

from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp


def dist(coor1, coor2):
    return math.sqrt(sum((coor1[k] - coor2[k]) ** 2 for k in range(2)))*100


def create_data_model(nodes, V, depot='0'):
    """Stores the data for the problem."""
    data = {}
    data['time_matrix'] = [[dist(nodes[i]['coor'], nodes[j]['coor']) for j in nodes] for i in nodes]
    data['time_windows'] = [tuple(x*100 for x in nodes[i]['tw']) for i in nodes]
    data['num_vehicles'] = len(V)
    data['depot'] = depot
    return data


def get_solution(manager, routing, solution, V):
    """Prints solution on console."""
    time_dimension = routing.GetDimensionOrDie('Time')
    total_time = 0

    for vehicle_id in V:
        index = routing.Start(vehicle_id)

        while not routing.IsEnd(index):
            V[vehicle_id].append(manager.IndexToNode(index))
            index = solution.Value(routing.NextVar(index))

        time_var = time_dimension.CumulVar(index)
        V[vehicle_id].append(manager.IndexToNode(index))
        total_time += solution.Min(time_var)

    return V, total_time/100


def or_tools_main(v, nodes, depot='0'):
    """Solve the VRP with time windows."""
    # Instantiate the data problem.
    data = create_data_model(nodes, v, depot=depot)

    # Create the routing index manager.
    manager = pywrapcp.RoutingIndexManager(len(data['time_matrix']),
                                           data['num_vehicles'], int(data['depot']))

    # Create Routing Model.
    routing = pywrapcp.RoutingModel(manager)


    # Create and register a transit callback.
    def time_callback(from_index, to_index):
        """Returns the travel time between the two nodes."""
        # Convert from routing variable Index to time matrix NodeIndex.
        from_node = manager.IndexToNode(from_index)
        to_node = manager.IndexToNode(to_index)
        return data['time_matrix'][from_node][to_node]

    transit_callback_index = routing.RegisterTransitCallback(time_callback)

    # Define cost of each arc.
    routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

    # def demand_callback(from_index):
    #     """Returns the demand of the node."""
    #     # Convert from routing variable Index to demands NodeIndex.
    #     from_node = manager.IndexToNode(from_index)
    #     return data['demands'][from_node]
    #
    # demand_callback_index = routing.RegisterUnaryTransitCallback(
    #     demand_callback)
    # routing.AddDimensionWithVehicleCapacity(
    #     demand_callback_index,
    #     0,  # null capacity slack
    #     data['vehicle_capacities'],  # vehicle maximum capacities
    #     True,  # start cumul to zero
    #     'Capacity')

    # Add Time Windows constraint.
    time = 'Time'
    routing.AddDimension(
        transit_callback_index,
        10000,  # allow waiting time
        10000,  # maximum time per vehicle
        False,  # Don't force start cumul to zero.
        time)
    time_dimension = routing.GetDimensionOrDie(time)

    # Add time window constraints for each location except depot.
    for location_idx, time_window in enumerate(data['time_windows']):
        if location_idx == 0:
            continue
        index = manager.NodeToIndex(location_idx)
        time_dimension.CumulVar(index).SetRange(time_window[0], time_window[1])
    # Add time window constraints for each vehicle start node.
    for vehicle_id in range(data['num_vehicles']):
        index = routing.Start(vehicle_id)
        time_dimension.CumulVar(index).SetRange(data['time_windows'][0][0],
                                                data['time_windows'][0][1])

    # Instantiate route start and end times to produce feasible times.
    for i in range(data['num_vehicles']):
        routing.AddVariableMinimizedByFinalizer(
            time_dimension.CumulVar(routing.Start(i)))
        routing.AddVariableMinimizedByFinalizer(
            time_dimension.CumulVar(routing.End(i)))

    # Setting first solution heuristic.
    search_parameters = pywrapcp.DefaultRoutingSearchParameters()
    search_parameters.first_solution_strategy = (
        routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC)

    # Solve the problem.
    solution = routing.SolveWithParameters(search_parameters)
    return [manager, routing, solution]



if __name__ == '__main__':
    main()