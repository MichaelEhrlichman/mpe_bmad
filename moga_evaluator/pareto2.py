#!/usr/bin/env python

import numpy as np

popfile = 'objective_report.last'
indexes = (1,2,3,4)
#popfile = 'moga_results.last'
#indexes = (0,19,20,21,22,23)

def pareto_dumb(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
        is_efficient[i] = np.all(np.any(costs[:i]>c, axis=1)) and np.all(np.any(costs[i+1:]>c, axis=1))
    return is_efficient

def pareto(costs, return_mask = True):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :param return_mask: True to return a mask
    :return: An array of indices of pareto-efficient points.
        If return_mask is True, this will be an (n_points, ) boolean array
        Otherwise it will be a (n_efficient_points, ) integer array of indices.
    """
    is_efficient = np.arange(costs.shape[0])
    n_points = costs.shape[0]
    next_point_index = 0  # Next index in the is_efficient array to search for
    while next_point_index<len(costs):
        nondominated_point_mask = np.any(costs<costs[next_point_index], axis=1)
        nondominated_point_mask[next_point_index] = True
        is_efficient = is_efficient[nondominated_point_mask]  # Remove dominated points
        costs = costs[nondominated_point_mask]
        next_point_index = np.sum(nondominated_point_mask[:next_point_index])+1
    if return_mask:
        is_efficient_mask = np.zeros(n_points, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient

npop = len(open(popfile).readlines())

population = np.zeros((npop,4))
i=0
with open(popfile,'r') as f:
	for l in f:
		sl = l.split()
		population[i] = [sl[n] for n in indexes]
		i+=1

#nondominated = pareto_dumb(population[:,:])
nondominated = pareto(population[:,:],False)

print(len(nondominated))

#for i in range(npop):
#	print('{} {}'.format(population[i,0],pareto_mask[i]))
	
