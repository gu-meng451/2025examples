using LinearAlgebra

A = [2 1 -1
    -3 -1 2
    -2 1 2.0]
b = [8; -11; -3.0]

## Augment the system [A b I]
B = hcat(A, b, I)

## Gaussian Elimination
B[3, :] = B[3, :] + B[1, :]
B[2, :] = B[2, :] + 3 / 2 * B[1, :]
B[3, :] = B[3, :] - 4 * B[2, :]
B[3, :] = -B[3, :]
B[2, :] = B[2, :] - 0.5 * B[3, :]
B[1, :] = B[1, :] + B[3, :]
B[2, :] = B[2, :] * 2
B[1, :] = B[1, :] - B[2, :]
B[1, :] = B[1, :] / 2

## extract the column that held b, it should be x
x = B[:,4]
Ainv = B[:,5:7]

## Let's check if it is a solution to Ax=b
norm( A*x - b )

## We can use the built-in solver to get x as well
x = A\b

## Let's check what's going on with Ainv
Ainv*A - I |> norm
A*Ainv - I |> norm

## We don't typically want Ainv, but should use the \ operator instead

## using Gaussian Elimination to build LU:
A = [2 1 -1
    -3 -1 2
    -2 1 2.0]

## Augment the system [A b I]
U = copy(A)
L = diagm([1,1,1.])

## Gaussian Elimination
# 1: (2,1)
U[2, :] += 3/2*U[1, :]
L[2, 1] = -3/2

# 2: (3,1)
U[3,:] += +1*U[1,:]
L[3,1] = -1

# 3: (3,2)
U[3,:] += -4*U[2,:]
L[3, 2] = 4

# check to make sure we did that right:
L*U - A |> norm

# Now that it's factored, we can solve (LU*x = b)
# first solve L*y = b (letting U*x = y).  This is the forward solve
y = L\b
# next do a back solve to get the final result
x = U\y

# check to make sure we solved the problem
A*x - b |> norm

## Built in LU
F = lu(A)

F.L*F.U - F.P*A |> norm
