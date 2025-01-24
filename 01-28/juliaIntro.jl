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


## Memory Scope
# for loops and simple functions

# for loops:


# list Comprehension


## while loops


## Modules (groupings of functions)
