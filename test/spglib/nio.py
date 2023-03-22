import numpy as np
import spglib

# lattice constant
a = 4.19

# structural information ( [lattice constants], [coordination], [atomic numbers] )
# Caution: it should be a tuple
nio = \
    ( [ [ a, 0, 0 ], [ 0, a, 0 ], [ 0, 0, a ] ],
      [ [ 0,   0,   0   ],  # Ni
        [ 0,   0.5, 0.5 ],  # Ni
        [ 0.5, 0,   0.5 ],  # Ni
        [ 0.5, 0.5, 0   ],  # Ni
        [ 0.5, 0,   0   ],  # O
        [ 0,   0.5, 0   ],  # O
        [ 0,   0,   0.5 ],  # O
        [ 0.5, 0.5, 0.5 ]   # O
      ],
      [28,] * 4 + [8,] * 4
    )

# just space group 
print( spglib.get_spacegroup( nio ) )
#Fm-3m (225)

# primitive cell
print( spglib.find_primitive( nio ) )
#(array([[0.   , 2.095, 2.095],
#       [2.095, 0.   , 2.095],
#       [2.095, 2.095, 0.   ]]), array([[0. , 0. , 0. ],
#       [0.5, 0.5, 0.5]]), array([28,  8], dtype=int32))
