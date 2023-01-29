from constraint import *
from collections import defaultdict

problem = Problem()

problem.addVariable('9-2', [0, 1])
problem.addVariable('10-2', [0, 1])
problem.addVariable('11-0', [0, 1])
problem.addVariable('11-1', [0, 1])
problem.addVariable('11-2', [0, 1])

ls = ['9-2', '10-2', '11-0', '11-1', '11-2']
problem.addConstraint(ExactSumConstraint(1), ['9-2', '10-2'])
problem.addConstraint(ExactSumConstraint(2), ['9-2', '10-2', '11-2'])
problem.addConstraint(ExactSumConstraint(1), ['11-0', '11-1'])
problem.addConstraint(ExactSumConstraint(3), ls)

solutions = problem.getSolutions()
sol_dict = {'9-2':0, '10-2':0, '11-0':0, '11-1':0, '11-2':0}
for s in solutions:
    for var, value  in s.items():
        sol_dict[var] += value

for var, count in sol_dict.items():
    sol_dict[var] = count / len(solutions)


print(sol_dict)
