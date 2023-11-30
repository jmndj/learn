#!/usr/bin/python
def interpolation_lagrange(x,y,u):
    n=len(x)
    l=[]
    for i in range(0,n):
        v=1
        for j in range(0,n):
            if j == i:
                continue
            v=v*(u-x[j])/(x[i]-x[j])
        l.append(v)
    v=0
    for i in range(0,n):
        v=v+y[i]*l[i]
    return v
x=[100,121,144]
y=[10,11,12]
print(interpolation_lagrange(x,y,115))        




    
