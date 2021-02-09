#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC

# Solve a traveling salesman problem on a randomly generated set of
# points using lazy constraints.   The base MIP model only includes
# 'degree-2' constraints, requiring each node to have exactly
# two incident edges.  Solutions to this model may contain subtours -
# tours that don't visit every city.  The lazy constraint callback
# adds new constraints to cut them off.
import copy
import sys
import math
import random
from itertools import combinations
import gurobipy as gp
from gurobipy import GRB


class VRPModel(gp.Model):
    def __init__(self, v, nodes, depot):
        gp.Model.__init__(self)
        self._nodes = copy.deepcopy(nodes)
        self._V = v
        self._depot = depot

    @property
    def nodes(self):
        return self._nodes

    @property
    def V(self):
        return self._V

    @property
    def depot(self):
        return self._depot


class VRP:
    def __init__(self, v, nodes, depot='0', Q_var=False):
        self.m = VRPModel(v, nodes, depot)
        #self.m.setParam('MIPGap', 10 ^ -1)
        if Q_var:
            self.m.nodes.update({f'aux{veh}': {'coor': (0, 0), 'tw': (0,1000)} for veh in v})

        self.m._vars = self.m.addVars(self.m.nodes, self.m.nodes, self.m.V, vtype=GRB.BINARY, name='x')
        self.T_vars = self.m.addVars(self.m.nodes, self.m.V, vtype=GRB.CONTINUOUS, name='T')
        self.dist = {
            (i, j): math.sqrt(
                sum((self.m.nodes[i]['coor'][k] - self.m.nodes[j]['coor'][k]) ** 2 for k in range(2))
            ) for i in self.m.nodes for j in self.m.nodes
        }

    def init_model(self):
        self.set_flow_constraints()
        #self.set_tw_constraints()
        self.set_depot_constraints()
        self.m.Params.lazyConstraints = 1
        self.m.setObjective(gp.quicksum(
            self.m._vars[i, j, v] * self.dist[(i, j)] for i in self.m.nodes for j in self.m.nodes for v in self.m.V
        ), GRB.MINIMIZE)
        self.m.update()

    def optim(self):
        self.m.optimize(self.subtourelim)
        self.m.write('out.mps')
        self.m.write('out.sol')
        vals = self.m.getAttr('x', self.m._vars)
        solution = {
            v: get_cycle(
                gp.tuplelist(
                    (i, j) for i, j, v1 in vals.keys() if vals[i, j, v1] > 0.5 and v1 == v
                ),
                begin=self.m.depot
            )[0]
            for v in self.m.V
        }
        return {v: [x for x in solution[v] if 'aux' not in x] + [self.m.depot] for v in self.m.V}, self.m.ObjVal

    def set_flow_constraints(self):
        self.m.addConstrs((self.m._vars.sum(i, '*', '*') == 1 for i in self.m.nodes if i != self.m.depot), name='out')
        self.m.addConstrs((
            (self.m._vars.sum(i, '*', v) - self.m._vars.sum('*', i, v) == 0)
            for i in self.m.nodes if i != self.m.depot
            for v in self.m.V
        ), name='in')
        self.m.addConstrs((self.m._vars.sum(i, i, v) == 0 for i in self.m.nodes for v in self.m.V), name='out')

    def set_tw_constraints(self):
        self.m.addConstrs((
            self.T_vars[i, v] + self.dist[i, j] - self.T_vars[j, v] <= (1-self.m._vars[i, j, v]) * 10000
            for i in self.m.nodes for j in self.m.nodes for v in self.m.V
        ), name='tw_flow')
        self.m.addConstrs((
            self.T_vars[i, v] >= self.m.nodes[i]['tw'][0] for i in self.m.nodes for v in self.m.V
        ), name='tw')
        self.m.addConstrs((
            self.T_vars[i, v] <= self.m.nodes[i]['tw'][1] for i in self.m.nodes for v in self.m.V
        ), name='tw')

    def set_depot_constraints(self):
        self.m.addConstrs((self.m._vars.sum(self.m.depot, '*', v) == 1 for v in self.m.V), name='depot_out')
        self.m.addConstrs((self.m._vars.sum('*', self.m.depot, v) == 1 for v in self.m.V), name='depot_in')

    @staticmethod
    def subtourelim(model, where):
        if where == GRB.Callback.MIPSOL:
            vals = model.cbGetSolution(model._vars)
            for v in model.V:
                selected = gp.tuplelist(
                    (i, j) for i, j, v1 in model._vars.keys() if vals[i, j, v1] > 0.5 and v1 == v
                )
                tour, unvisited = get_cycle(selected)
                tour_edges = list(zip(tour, tour[1:])) + [(tour[-1], tour[0])]

                if len(unvisited) > 0:
                    for v1 in model.V:
                        if model.depot in unvisited:
                            model.cbLazy(gp.quicksum(model._vars[i, j, v1] for i, j in tour_edges) <= len(tour)-1)
                        else:
                            model.cbLazy(
                                gp.quicksum(
                                    model._vars[i, edge[1], v1] for i in unvisited for edge in selected if edge[0] == i
                                )
                                <= len(unvisited) - 1
                            )
                else:
                    pass


def get_cycle(edges, begin=False):
    unvisited = []

    for edge in edges:
        if edge[0] not in unvisited:
            unvisited.append(edge[0])
        if edge[1] not in unvisited:
            unvisited.append(edge[1])

    if begin:
        cycle = [begin]
        unvisited.remove(begin)
    else:
        cycle = [unvisited.pop(0)]

    while edges:  # true if list is non-empty
        next_edge = next(edge for edge in edges if edge[0] == cycle[-1])

        if next_edge[1] == cycle[0]:
            return cycle, unvisited

        cycle.append(next_edge[1])
        unvisited.remove(next_edge[1])
        edges.remove(next_edge)

    raise Exception("Solution it's not a cycle")
