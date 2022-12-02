# Test of L1 minimization to detect errors in graph

# Author: Etienne Sandier

# We import cvxpy
import cvxpy as cp
# We import the networkX package for graph analysis
import networkx as nx
# import numpy
import numpy as np

# create a graph with weights
G = nx.DiGraph()
G.add_edges_from(
    [
        (0,1,{'weight':1}),
        (1,2, {'weight':0}), 
        (2,0,{'weight':-1}),
        (0,3,{'weight':1}),
        (3,4,{'weight':0}), 
        (4,0,{'weight':-1}),
        (1,4,{'weight':2}),
    ]
)
edge_index = {}
for i, edge in enumerate(G.edges):
    edge_index[edge] = i
    # print("index %i, edge (%i, %i), weight %i"%(i, *edge, G.edges[edge]['weight']))

H = G.to_undirected()
cycle_basis = nx.cycle_basis(H)

weights = cp.Variable(len(G.edges))

constraints = []
for cycle in cycle_basis:
    LHS = 0
    for i in range(len(cycle)):
        edge = (cycle[i], cycle[(i+1) % len(cycle)])
        # print(edge)
        if edge in G.edges:
            LHS += weights[edge_index[edge]]
        else:
            LHS -= weights[edge_index[edge[::-1]]]
    constraints.append(LHS == 0)
weight_list = np.empty(len(G.edges))
for edge in G.edges:
    weight_list[edge_index[edge]] = G.edges[edge]["weight"]

objective = cp.Minimize(cp.norm(weights - weight_list, 1))
        
prob = cp.Problem(objective, constraints)

result = prob.solve(verbose = True)#, solver = cp.ECOS)

for i, edge in enumerate(G.edges):
    print("index %i, edge (%i, %i), initial weight %i, solution %f"%(i, *edge, weight_list[i], weights.value[i]))

# print(cycle_basis)    
# print(constraints)


# for i, edge in enumerate(G.edges):
#     print("solution: %i, initial weight: %i"%(weights.value[i], weight_list[i]))


