# This work is licensed under a Creative Commons Attribution 4.0 International License
# https://creativecommons.org/licenses/by/4.0/
# Copyright (c) 2021, Julien Langou. All rights reserved.
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
#
def lu_with_partial_pivoting( A_, inplace = False ):
#
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
