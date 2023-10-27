#Dwave imports
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod
from __future__ import print_function

import dimod
import math
import sys
import copy

from dimod.generators.constraints import combinations
from hybrid.reference import KerberosSampler

#Impots for QUBO probelm
import numpy as np
from numpy import linalg as LA
import pandas as pd
import matplotlib as plt
import itertools as it
from itertools import product
from matplotlib import pyplot as plt
import math 



#Dwave imports
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod
from __future__ import print_function

import dimod
import math
import sys
import copy

from dimod.generators.constraints import combinations
from hybrid.reference import KerberosSampler

#Imports for QUBO probelm
import numpy as np
from numpy import linalg as LA
import pandas as pd
import matplotlib as plt
import itertools as it
from itertools import product
from matplotlib import pyplot as plt
import math

#check rows and cols:
def check_OLS(A,k_k):
    set_k_k = set(k_k)
    pass_flag = True
    f = 0
    for l in range(A.shape[0]):
        for i in range(A.shape[1]):
            if len(set(A[l,i,:])) != A.shape[1]:
            
                print(f'layer {l}, row {i} fails')
                pass_flag = False
                f += 1
                #print(len(set(A[0,i,:])), ' : ' ,A.shape[1])
            else:
                print(f'layer {l}, row {i} passes') 


        for j in range(A.shape[2]):
            if len(set(A[l,:,j])) != A.shape[2]:
                #print(len(set(A[0,i,:])), ' : ' ,A.shape[2])
                print(f'layer {l}, column {j} fails')
                pass_flag = False
                f += 1

            else:
                print(f'layer {l}, column {j} passes') 
    print(f, ' rows or columns have failed')
    #Check tuples
    if len(set_k_k) == len(k_k):
        print('All tuples are unique')
    else:
        print('Non-unique tuples exist.', 'Total tuples: ', len(k_k), ' Unique tuples: ' , len(set_k_k))
    return pass_flag


### s_i_j_tuple formulation: 

def linear_bqm(digits,cols,rows):

    # Set up
    L = range(digits)
    K = range(digits)
    X = range(cols)
    Y = range(rows)
    
    bqm = dimod.BinaryQuadraticModel({}, {}, 0.0, dimod.BINARY)
    bqm_dict = {}

    print('BQM started')

    # (1): Every l,k has a i,j
    print('(1):')
    for l in L:
        for k in K:
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (x,y) in it.product(X,Y)], 1))

    # (2): 
    print('(2):')
    for y in Y:
        for l in L:
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (x,k) in it.product(X,K)], 1))
            #print([f's_{x}_{y}_{k}_{l}' for (x,k) in it.product(X,K)])

    # (3): 
    print('(3):')
    for y in Y:
        for k in K:
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (x,l) in it.product(X,L)], 1))
            #print([f's_{x}_{y}_{k}_{l}' for (x,l) in it.product(X,L)])

    # (4):
    print('(4):')
    for x in X:
        for l in L:
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (y,k) in it.product(Y,K)], 1))
            #print([f's_{x}_{y}_{k}_{l}' for (y,k) in it.product(Y,K)])

    # (5): 
    print('(5):')
    for x in X:
        for k in K:
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (y,l) in it.product(Y,L)], 1))
            #print([f's_{x}_{y}_{k}_{l}' for (y,l) in it.product(Y,L)])

    # (6): every box has an l,k
    for x in X:
        for y in Y: 
            bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (l,k) in it.product(L,K)], 1))
  
    # # x*y binary vars

    bqm.update(combinations([f's_{x}_{y}_{k}_{l}' for (x,y,l,k) in it.product(X,Y,L,K)], cols * rows))
    
    #print(bqm)


    solution = KerberosSampler().sample(bqm, num_reads = 10)
    df = solution.to_pandas_dataframe()
    best_solution = solution.first.sample
    solution_list = [k for k, v in best_solution.items() if v == 1]
            #solutions.append(label.split('_')[1:])
    solutions = [label.split('_')[1:] for label in solution_list]

    simplified_coords = []
    k_k = []
    for s in solutions:
        k_k.append((int(s[2]),int(s[3])))

    A = np.zeros((2, x_card, y_card))
    for sol in solutions:
        for l in range(2):
            A[l,int(sol[0]),int(sol[1])] = int(sol[2+l])
    return df, A, k_k, solutions

digits = 7 #possible digits
x_card = 7 #x length
y_card = x_card #y length

true_solutons = []
a = False
count = 0
while a == False:
    df, A, k_k, solutions = linear_bqm(digits, x_card, y_card)
    print('BQM complete')
    a = check_OLS(A,k_k)
    print('itter: ', count)
    count += 1
    if a == True:
        true_solutons.append([df, A, k_k])
        print('True solution : ' , A)



