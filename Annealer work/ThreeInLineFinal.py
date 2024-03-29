#Dwave imports
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite
import dimod

#Impots for QUBO probelm
import numpy as np
from numpy import linalg as LA
import pandas as pd
import matplotlib as plt
import itertools
from itertools import product
from matplotlib import pyplot as plt
import math
 
def sample_to_img(df):
    num = df[['energy']].idxmin()[0]
    bit_string = np.zeros(n**2)
    for i in range(df.shape[1]-2):
        if df.iloc[num][f's_{i}'] == 1:
            bit_string[i] = 1
    make_img_bit_string(bit_string)


def make_img_bit_string(b):
    plt.imshow(b.reshape((n,n)))
    plt.show()
    print(b)


def line(point1,point2):
    # ax+by+c=0 from (p1-p)(p1-p2)'=0
    a=point2[1]-point1[1]; b=point1[0]-point2[0]; c=point1[1]*point2[0]-point1[0]*point2[1]
    f=math.gcd(a,b,c)
    f=1 if f==0 else -f if a<0 else f if a>0 else -f if b<0 else f
    return (a//f,b//f,c//f)


def get_linear_sets(n):
    #Get list of cart coords
    xs = range(n)
    ys = range(n)
    X = product(xs,ys)
    X = list(X)
    lines = [] #store lines
    binary_sets = []  #store sets of binary variables corresponding to each line
    line_set_dict = {} #store binary variables and lines
    for k in range(len(X)): 
        for d in range(k+1,len(X)): #iterate over points after the kth, prevents repeating the same operation
            line_ = (line(X[k],X[d])) #define a line
            if line_ in lines:
                continue #check if this line has been done
            else:
                lines.append(line_) #add new line
                my_set = []
                my_set.append(X[k]) #add initial points
                my_set.append(X[d])
                x = X[d][0] #set new x,y values
                y = X[d][1]
                dx = line_[1] #define dx,dy
                dy = -line_[0]
                #print('(x,y)_1 = ', (x,y))
                if dx < 0:
                    dx = -dx
                    dy = -dy
                if dx == 0:
                    dy = abs(dy)
                #print('(dx,dy) = ', (dx,dy))
                x += dx #find next point of the line
                y += dy
                #my_set.append((x,y)) 
                l = 1
                while 0 <= x < n and 0 <= y < n:
                    l += 1
                    my_set.append((x,y)) #add the next point on the line until the boundary is exceeded
                    #print(f'(x,y)_{l} = ', (x,y))
                    
                    x += dx
                    y += dy
                    
                    
                #my_set = set(my_set)
                binary_sets.append(my_set) #add the points to a list
                #print(line_, ':' ,my_set)
                line_set_dict.update({line_ : my_set}) #The will be a tuple represnting a line as the key and a set of points as the value
    return line_set_dict

def get_interation_strengths(n):    
    line_set_dictionary = get_linear_sets(n)
    sets = list(line_set_dictionary.values())
    X_tri_int = np.zeros((n**2,n**2,n**2),dtype = 'int32')
    for i in range(len(sets)):
        s = sets[i]
        if len(s) > 2:
            ind_set = [a for a in itertools.combinations([a for a in range(len(s))], 3)]
            #print(ind_set,s)
            for A in ind_set:
                x = s[A[0]][0]*n + s[A[0]][1] #x = in + j 
                y = s[A[1]][0]*n + s[A[1]][1]
                z = s[A[2]][0]*n + s[A[2]][1]
                #print('(x,y,z) = ' , (x,y,z))
                X_tri_int[x,y,z] += 1
    return X_tri_int
               
n = 5

X_tri_int = get_interation_strengths(5)
X_tri_int[0,1,2]

#define binary decision vector.
x = [f's_{i}' for i in range(n**2)]
int = [(f's_{i}',f's_{j}',f's_{k}') for (i,j,k) in itertools.combinations(range(n**2),3)]

int_dict = {}
for i in range(n**2):
    int_dict.update({(f's_{i}',f's_{i}') : -1})
for (i,j,k) in itertools.combinations(range(n**2),3):
    if X_tri_int[i,j,k] != 0:
        int_dict.update({(f's_{i}',f's_{j}',f's_{k}') : X_tri_int[i,j,k]})

df = dimod.ExactPolySolver().sample_hubo(int_dict).to_pandas_dataframe()

sample_to_img(df)





check_OLS(A,k_k)