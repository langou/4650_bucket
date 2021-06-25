##############################################################################################
#
# This work is licensed under a Creative Commons Attribution 4.0 International License
# https://creativecommons.org/licenses/by/4.0/
# Copyright (c) 2021, Julien Langou. All rights reserved.
#
##############################################################################################
##############################################################################################
#
# lu_with_partial_pivoting
#
# friendly useage
#
#     P, L, U = lu_with_partial_pivoting( A )
#
#        A:     ( input )  matrix to be `PA=LU`-factored,
#               A is not changed on exit
#        P,L,U: ( output ) matrices such that (1) PA=LU, (2) U is upper triangular, 
#               (3) P is a permutation matrix and (4) L is unit lower triangular
#               matrix with elements below diagonal of absolute value less than 
#               or equal to 1.
#
#  advanced useage
#
#     perm = lu_with_partial_pivoting( A, inplace = True )
#
#        A:    ( input/ouput ) in input, matrix to be `PA=LU`-factored. In output, 
#              L and U stacked on top of the other in array A
#        perm: ( ouptut ) permutation vector representating the row permutations
#
def lu_with_partial_pivoting( A_, inplace = False ):
#
  import numpy as np
  import copy
#    
  if ( inplace ):
    A = A_
  else: 
    A = copy.deepcopy(A_)
#
# begin main code
#  
  n = A.shape[0]
  
  perm = np.arange(n)

  for k in range(0,n-1):

    i = k + np.argmax( abs( A[k:n,k] ) )
    A[[k,i],:] = A[[i,k],:]
    perm[[k,i]] = perm[[i,k]]

    if ( A[k,k] == 0. ): print("oops, matrix is singular\n"); break

    for i in range(k+1,n):
      A[i,k] = A[i,k] / A[k,k]
      for j in range(k+1,n):
        A[i,j] = A[i,j] - A[k,j] * A[i,k]
#        
# end main code   
#
  if inplace:
    return perm
  else:
    L = np.tril(A,-1)+np.eye(n)
    U = np.triu(A)
    P = np.eye(n)[perm,:]
    return P, L, U
#
##############################################################################################
#
# backward_substitution
#
#
# perform backward substitution Ax = b where A is upper triangular
#
# Note that we do not check whether A is upper triangular, it is assumed to be,
# Only the upper part of A is referenced. The zeros in the lower part are assumed.
# You can store whatever you want in the lower part of A.
#
# useage: 
# => for the ``standard interface``, do
#           x = backward_substitution( A, b )
#
# => for the ``in-place`` interface, do
#           backward_substitution( A, b, inplace = True )
#
# input:  A and b: represent the upper triangular linear system of equations 
#                  Ax=b we want to solve
# output:          solution of x of Ax=b, 
#    
def backward_substitution( A, b_, inplace = False ):
#
  import numpy as np
  import copy
#    
  if ( inplace ):
    b = b_
  else: 
    b = copy.deepcopy(b_)
#    
  n = A.shape[0]
  for i in range(n-1,-1,-1):
    for j in range(i+1,n):
      b[i] = b[i] - A[i,j] * b[j]
    b[i] = b[i] / A[i,i]
#    
  if not inplace:
    return b
#
##############################################################################################
#
# forward__substitution
#
# This work is licensed under a Creative Commons Attribution 4.0 International License
# https://creativecommons.org/licenses/by/4.0/
# Copyright (c) 2021, Julien Langou. All rights reserved.
#
# useage: 
# => for the ``standard interface``, do
#           x = forward__substitution___x( A, b )
#
# => for the ``in-place`` interface, do
#           forward__substitution___x( A, b, inplace = True )
#
# the flag ``unit`` is to say whether there are one on the diagonal of A.
# If unit is ``True`` then 1 are assumed on the diagonal of A. We do not reference the diagonal
# and you can store whatever in it.
#
# Note that we do not check whether A is lower triangular, it is assumed to be.
# Only the lower part of A is referenced. The zeros in the upper part are assumed.
# You can store whatever you want in the upper part of A.
#
def forward__substitution( A, b_, inplace = False, unit = False ):
#
  import numpy as np
  import copy
#
  if ( inplace ):
    b = b_
  else: 
    b = copy.deepcopy(b_)
#    
  n = A.shape[0]
  for i in range(0,n,):
    for j in range(0,i):
      b[i] = b[i] - A[i,j] * b[j]
    if not unit:
      b[i] = b[i] / A[i,i]
#    
  if not inplace:
    return b
