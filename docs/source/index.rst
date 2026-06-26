burnsidelib
===========

``burnsidelib`` is a small, single-file SageMath library for computations in
Burnside rings of finite groups.

It is designed for a classroom workflow: students copy ``burnsidelib.py`` into
the same folder as their Sage worksheet or script and load it directly.

Quick start
-----------

::

    sage: load("burnsidelib.py")
    sage: C2 = gap("CyclicGroup(2)")
    sage: A = BurnsideRing(C2)
    sage: A.basis_names()
    ('[C2/1]', '[C2/C2]')
    sage: 3*A.gen(0) - 4*A.one()
    3*[C2/1] - 4*[C2/C2]

The main object is :class:`burnsidelib.BurnsideRing`, whose elements support
marks, restriction, transfer, norms, and several power operations.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   tutorial
   api
