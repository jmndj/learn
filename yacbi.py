#!/usr/bin/python
import math
import numpy as np
def error(x,y):
    sum =0 
    for i in range(len(x)):
        sum = sum + (x[i]-y[i])**2
    return math.sqrt(sum)
t = 1
x=[1,1,1]
y=np.zeros(3)
z=np.array([[8,-3,2,20],[4,11,-1,33],[2,1,4,12]])
while (t != 0 ):
    for i in range(len(x)):
        p = 0
        for j in range(len(x)):
            if j !=i:
                p = p - z[i][j]*x[j]+z[i][3]
        y[i]= p/z[i][i]
    t = error(x,y)
    for i in range(len(x)):
        x[i]=y[i]
    print(y)   
