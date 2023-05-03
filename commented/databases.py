# -*- coding: utf-8 -*-
r"""
Access functions to online databases for coding theory
"""



def self_orthogonal_binary_codes(n, k, b=2, parent=None, BC=None, equal=False,
    in_test=None):
    """
    Returns a Python iterator which generates a complete set of
    representatives of all permutation equivalence classes of
    self-orthogonal binary linear codes of length in ``[1..n]`` and
    dimension in ``[1..k]``.

    INPUT:

    -  ``n`` - Integer, maximal length

    -  ``k`` - Integer, maximal dimension

    -  ``b`` - Integer, requires that the generators all have weight divisible
       by ``b`` (if ``b=2``, all self-orthogonal codes are generated, and if
       ``b=4``, all doubly even codes are generated). Must be an even positive
       integer.

    -  ``parent`` - Used in recursion (default: ``None``)

    -  ``BC`` - Used in recursion (default: ``None``)

    -  ``equal`` - If ``True`` generates only [n, k] codes (default: ``False``)

    -  ``in_test`` - Used in recursion (default: ``None``)

    EXAMPLES:

    Generate all self-orthogonal codes of length up to 7 and dimension up
    to 3::

        sage: for B in codes.databases.self_orthogonal_binary_codes(7,3):
        ....:    print(B)
        [2, 1] linear code over GF(2)
        [4, 2] linear code over GF(2)
        [6, 3] linear code over GF(2)
        [4, 1] linear code over GF(2)
        [6, 2] linear code over GF(2)
        [6, 2] linear code over GF(2)
        [7, 3] linear code over GF(2)
        [6, 1] linear code over GF(2)

    Generate all doubly-even codes of length up to 7 and dimension up
    to 3::

        sage: for B in codes.databases.self_orthogonal_binary_codes(7,3,4):
        ....:    print(B); print(B.generator_matrix())
        [4, 1] linear code over GF(2)
        [1 1 1 1]
        [6, 2] linear code over GF(2)
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]
        [7, 3] linear code over GF(2)
        [1 0 1 1 0 1 0]
        [0 1 0 1 1 1 0]
        [0 0 1 0 1 1 1]

    Generate all doubly-even codes of length up to 7 and dimension up
    to 2::

        sage: for B in codes.databases.self_orthogonal_binary_codes(7,2,4):
        ....:    print(B); print(B.generator_matrix())
        [4, 1] linear code over GF(2)
        [1 1 1 1]
        [6, 2] linear code over GF(2)
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]

    Generate all self-orthogonal codes of length equal to 8 and
    dimension equal to 4::

        sage: for B in codes.databases.self_orthogonal_binary_codes(8, 4, equal=True):
        ....:     print(B); print(B.generator_matrix())
        [8, 4] linear code over GF(2)
        [1 0 0 1 0 0 0 0]
        [0 1 0 0 1 0 0 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 0 0 0 1 1]
        [8, 4] linear code over GF(2)
        [1 0 0 1 1 0 1 0]
        [0 1 0 1 1 1 0 0]
        [0 0 1 0 1 1 1 0]
        [0 0 0 1 0 1 1 1]

    Since all the codes will be self-orthogonal, b must be divisible by
    2::

        sage: list(codes.databases.self_orthogonal_binary_codes(8, 4, 1, equal=True))
        Traceback (most recent call last):
        ...
        ValueError: b (1) must be a positive even integer.
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.matrix.constructor import Matrix

    d=int(b)#this is 4 for doubly-even codes
    if d!=b or d%2==1 or d <= 0:
        raise ValueError("b (%s) must be a positive even integer."%b)
    from .linear_code import LinearCode
    from .binary_code import BinaryCode, BinaryCodeClassifier
    
    # first optimization: separate function assuming results of these
    #                     if statements don't change between recursions
    
    if k < 1 or n < 2:
        return
    
    # these if-else blocks overwrite the input in_test function
    if equal:
        in_test = lambda M: (M.ncols() - M.nrows()) <= (n-k)
        # test for recursion: see if parent code has size 
        # similar to size of codes we want
        out_test = lambda C: (C.dimension() == k) and (C.length() == n)
        # test for output: we only output codes that have the specified
        # dimensions that we want
    else:
        # maybe add:
        # if in_test != None:
        in_test = lambda M: True
        out_test = lambda C: True
    if BC is None:
        BC = BinaryCodeClassifier()
    if parent is None:
        # initial function call before recursion
        for j in range(d, n+1, d):
            M = Matrix(FiniteField(2), [[1]*j])
            # initial parents are 1 dimensional
            # one parent for each multiple of d (4 for doubly even)
            if in_test(M):
                for N in self_orthogonal_binary_codes(n, k, d, M, BC, in_test=in_test):
                    if out_test(N):
                        yield N
    else: # recursion
        C = LinearCode(parent)
        
        if out_test(C):
            yield C
        if k == parent.nrows():
            return
        # exit condition: the parent code has maximal dimension
        for nn in range(parent.ncols()+1, n+1):
            # iterate over codes with degree greater than the parent
            # and leq given degree
            if in_test(parent): # we can move this if statement outside of the loop
                for child in BC.generate_children(BinaryCode(parent), nn, d):
                    for N in self_orthogonal_binary_codes(n, k, d, child, BC, in_test=in_test):
                        if out_test(N): 
                            # this test might be applied to some codes
                            # that have already passed it
                            yield N

