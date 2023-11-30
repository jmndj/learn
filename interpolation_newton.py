#!/usr/bin/python
def interpolation_newton(x,y,u):
    n=len(x)
    z=[]
    for i in range(0,n):
        z.append(y[i])
    for i in range(1,n):
        for j in range(i,n):
            y[j]=(z[j]-z[j-1])/(x[j]-x[j-i])
        for j in range(i,n):
            z[j]=y[j]
    v=0
    for i in range(n-1,-1,-1):
        v=v*(u-x[i])+y[i]
    return v,y
#x=[1,3,4,7]
#y=[0,2,15,12]
x=[100,121,144]
y=[10,11,12]
print(interpolation_newton(x,y,120)) 


