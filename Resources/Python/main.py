
import random
import time

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from OR_tools import or_tools_main, get_solution
from VRPTW import VRP


def draw_graph(nodes, solution):
    G = nx.Graph()

    nodes_coor = {node: nodes[node]['coor'] for node in nodes}
    G.add_nodes_from(nodes_coor.keys())

    for n in nodes_coor.keys():
        G.nodes[n]['nodes'] = nodes_coor[n]

    for vehicle in solution:
        for edge in zip(solution[vehicle], solution[vehicle][1:]):
            if vehicle == 0:
                color = 'blue'
            elif vehicle == 1:
                color = 'green'
            elif vehicle == 2:
                color = 'black'
            else:
                color = 'red'

            G.add_edge(edge[0], edge[1], color=color)

    colors = [G[i][j]['color'] for v in solution for i, j in zip(solution[v], solution[v][1:])]
    nx.draw(G, nodes_coor, with_labels=True, node_size=500, edge_color=colors)
    plt.show()


def solve_gurobi(V, nodes, depot='0'):
    vrptw = VRP(V, nodes, depot=depot, Q_var=True)
    vrptw.init_model()
    now = time.time()
    solution, objval = vrptw.optim()
    print(f"Gurobi's solution: {solution}, {objval}, {time.time()-now} seg")
    return solution


def solve_ortools(V, nodes, depot='0'):
    vrptw = or_tools_main(V, nodes, depot)
    now = time.time()
    solution = get_solution(vrptw[0], vrptw[1], vrptw[2], V)
    print(f"OR tool's solution: {solution}")
    print(f"OR tool time's solution: {time.time() - now} seg")

    return solution

if __name__ == '__main__':

    N_veh = 20
    V = {v: [] for v in range(N_veh)}
    N = 20
    np.random.seed(1)

    # Nodes = {
    #     str(n): {
    #         'coor': tuple(np.random.uniform(-10, 10, 2)) if n > 0 else (0, 0),
    #         'tw': (0, random.randint(15, 30) if n > 0 else 10000)
    #     } for n in range(N)
    # }

    time_windows = [
        (0, 1000),  # store
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

    coordinates = [
        (4.56, 3.20),   # location  0 - the depot
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
	]

    Nodes = {
        str(n): {
            'coor': coordinates[n],
            'tw': time_windows[n]
        } for n in range(len(coordinates))
    }

    
    # sol_gurobi = solve_gurobi(V, Nodes)

    sol_or_tools = solve_ortools(V, Nodes)


    #draw_graph(Nodes, sol_gurobi)
    # draw_graph(Nodes, {v: [str(sol) for sol in sol_or_tools[0][v]] for v in V})


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
