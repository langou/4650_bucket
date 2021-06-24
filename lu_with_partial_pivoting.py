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
#import copy
#
def lu_with_partial_pivoting( A_, inplace = False ):
  return 1
