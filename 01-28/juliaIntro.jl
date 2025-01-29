# using Pkg
# Pkg.activate(".")
# using Plots

using Printf
# using LinearAlgebra

## Strings
a = "this is a string"

## Dictionaries
# key => value
A = Dict( 1 => "value1", "2"=>rand(), :hamburger => t-> cos(t) )

# cos(1)
A[:hamburger](1)

## Strucs
mutable struct myCustomStruct
    prop1
    prop2::String
    prop3::Int
end

B = myCustomStruct(pi, "This is a String", 2)

## Tuple
C = (1,2,"a")


## functions

function myfunc(x)
    y = x+1
    return y
end


function myfunc(x, z)
    y = x+z
    return (y,z^2)
end

function myfunc(n::Int)
    return sin(n)
end

myfunc(1.)
d,e = myfunc(1,2)

# inline function
f(x) = cos(2*x) + 1

## extra keywords and default values
global w = 3
function mycalc(x,y,z; tol = 1e-6, α = 1)

    result = α*(x+y+z) + w
    
    if result > tol
        return result
    else
        return NaN
    end

end
mycalc(1,2,3, α = 1e-8)
mycalc(1,2,3)
w

function inplacechange!(A)

    # A[1] = A[1] + A[2]
    A[1] += A[2] 

end;
A = [1,2,3]
inplacechange!(A)
A

1+1

## Memory Scope
# for loops and simple functions

# for loops:
x = 0
for i = 1:4
    x =  i
end

## copies of arrays
A
B = A # makes a ref, not a copy
B[1] = 0

C = copy(A)
C[1] = -100

# list Comprehension
# for i = 1:n
#     x[i] = some work
# end
x = [ i^2*j for i in 1:4, j in 1:3 ]


## while loops
iter = 0
n = 3
while iter < n
    iter += 1
    @show iter
end


## what's up with the dot?
A = rand(3,2)
t = LinRange(0,20,100); # like linspace in Matlab

x = cos.(2*t)

A.^2/2

f.(A)
