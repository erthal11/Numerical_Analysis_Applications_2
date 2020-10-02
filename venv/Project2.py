import sympy as sym
from mpmath import *

e = mp.e

def newtonsdivideddifference(x, x0, fx0, x1, fx1, x2, fx2, x3, fx3):

    def orderfunc(xa, fxa, xb, fxb):
        return (fxb - fxa)/(xb-xa)

    fx0x1 = orderfunc(x0, fx0, x1, fx1)
    fx1x2 = orderfunc(x1, fx1, x2, fx2)
    fx2x3 = orderfunc(x2, fx2, x3, fx3)

    fx0x1x2 = (fx1x2 - fx0x1)/(x2-x0)
    fx1x2x3 = (fx2x3 - fx1x2)/(x3-x1)

    fx0x1x2x3 = (fx1x2x3 - fx0x1x2)/(x3 - x0)

    def order1(x): return fx0 + fx0x1*(x-x0)
    def order2(x): return order1(x) + fx0x1x2 * (x-x0) * (x-x1)
    def order3(x): return order2(x) + fx0x1x2x3 * (x-x0) * (x-x1) * (x-x2)

    return order3(x)


print("Newton's Divided Difference: ", newtonsdivideddifference(8.4, 8.1, 16.94410, 8.3, 17.56492, 8.6, 18.50515, 8.7, 18.82091))
print()

def hermiteinterpolation(x,x0,fx0, fpx0, x1, fx1, fpx1):

    z0 = x0
    z1 = x0
    z2 = x1
    z3 = x1

    fz0 = fx0
    fz1 = fx0
    fz2 = fx1
    fz3 = fx1

    def orderfunc(xa, fxa, xb, fxb):
        return (fxb - fxa)/(xb-xa)

    fx0x1 = fpx0
    fx1x2 = orderfunc(z1, fz1, z2, fz2)
    fx2x3 = fpx1

    fx0x1x2 = (fx1x2 - fx0x1)/(z2-z0)
    fx1x2x3 = (fx2x3 - fx1x2)/(z3-z1)

    fx0x1x2x3 = (fx1x2x3 - fx0x1x2)/(z3 - z0)

    def order1(x): return fz0 + fx0x1*(x-z0)
    def order2(x): return order1(x) + fx0x1x2 * (x-z0) * (x-z1)
    def order3(x): return order2(x) + fx0x1x2x3 * (x-z0) * (x-z1) * (x-z2)

    return order3(x)


print("Hermite Interpolation at t=10: ", hermiteinterpolation(10,8,623,74,13,993,72))
# divided difference table on written portion
print()


def funcb(x):
    return (pow(e,2*x)*sin(3*x))
    #integral: -14.21397713


# for compositeTrapezoidalRule
def summation(lower,upper,func,xs):
    summation = 0
    i = lower
    while i<=upper:
        summation += func(xs[i])
        i += 1
    return summation


def compositeTrapezoidalRule(a,b,func,n):

    h = (b-a)/n

    xs = []
    i = a
    while i<=b:
        xs.append(i)
        i = i + h

    return ((h/2)*(func(a) + 2*(summation(1,n-1,func,xs)) + func(b)))



# for compositeSimpsonsRule and compositeMidpointRule
def summation2(lower,upper,func,xs):
    summation = 0
    i = lower
    while i<=upper:
        summation += func(xs[2*i])
        i += 1
    return summation

# for compositeSimpsonsRule
def summation3(lower,upper,func,xs):
    summation = 0
    i = lower
    while i<=upper:
        if len(xs)>= (2*i):
            summation += func(xs[2*i-1])
        i += 1
    return summation


def compositeSimpsonsRule(a,b,func,n):

    h = (b-a)/n

    xs = []
    i = a
    while i <= b:
        xs.append(i)
        i = i + h

    return ((h/3)*(func(a) + 2*(summation2(1,(n/2)-1,func,xs)) + 4*(summation3(1,n/2,func,xs)) + func(b)))



def compositeMidpointRule(a,b,func,n):

    h = (b-a)/(n+2)

    xs = []
    i = a + h
    while i <= b:
        xs.append(i)
        i = i + h

    return (2*h*(summation2(0,n/2,func,xs)))




trap4 = compositeTrapezoidalRule(0,2,funcb,4)
simp4 = compositeSimpsonsRule(0,2,funcb,4)
mid4 = compositeMidpointRule(0,2,funcb,4)

trap8 = compositeTrapezoidalRule(0,2,funcb,8)
simp8 = compositeSimpsonsRule(0,2,funcb,8)
mid8 = compositeMidpointRule(0,2,funcb,8)

trap16 = compositeTrapezoidalRule(0,2,funcb,16)
simp16 = compositeSimpsonsRule(0,2,funcb,16)
mid16 = compositeMidpointRule(0,2,funcb,16)

exact = -14.21397713

print ("For n = 4:")
print ("Composite Trapezoidal Rule gives:", trap4, " , with absolute error: ", abs(trap4-exact))
print ("Composite Simpson's Rule gives:", simp4, " , with absolute error: ", abs(simp4-exact))
print ("Composite Midpoint Rule gives:", mid4, " , with absolute error: ", abs(mid4-exact))
print()

print ("For n = 8:")
print ("Composite Trapezoidal Rule gives:", trap8, " , with absolute error: ", abs(trap8-exact))
print("Composite Simpson's Rule gives:", simp8, " , with absolute error: ", abs(simp8-exact))
print("Composite Midpoint Rule gives:", mid8, " , with absolute error: ", abs(mid8-exact))
print()

print ("For n = 16:")
print ("Composite Trapezoidal Rule gives:", trap16, " , with absolute error: ", abs(trap16-exact))
print("Composite Simpson's Rule gives:", simp16, " , with absolute error: ", abs(simp16-exact))
print("Composite Midpoint Rule gives:", mid16, " , with absolute error: ", abs(mid16-exact))
print()
