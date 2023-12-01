import sympy as sp
import math
def fun(x,y):
    result = 0
    for i in range(len(x)):
        result = (x[i]+y[i])**2+result
    return math.sqrt(result)
x1,x2,a,y1,y2 = sp.symbols('x1 x2 a y1 y2')
f1= x1**2+25*x2**2
g1,g2=sp.diff(f1,x1),sp.diff(f1,x2)
#print(g1.subs({x1:1}),g2)
y1=x1-a*g1
y2=x2-a*g2
f1= y1**2+25*y2**2
ga=sp.diff(f1,a)
#print(ga)
x1num=2
x2num=2
iter = 0
t = 1
while t > 1e-6:
    w=sp.solve(ga,a)[0]
    amin=w.subs({x1:x1num,x2:x2num})
    d1=-g1.subs({x1:x1num,x2:x2num})
    d2=-g2.subs({x1:x1num,x2:x2num})
    x=[x1num,x2num]
    x1num=x1num+amin*d1
    x2num=x2num+amin*d2
    y=[float(x1num),float(x2num)]
    t = fun(x,y)
    iter = iter+1
print(y,iter)