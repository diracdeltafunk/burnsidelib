Tutorial
========

.. currentmodule:: burnsidelib

Creating a Burnside ring
------------------------

Start Sage in the folder containing ``burnsidelib.py`` and load the file::

    sage: load("burnsidelib.py")

Create a group with GAP, then pass it to :class:`BurnsideRing`::

    sage: C2 = gap("CyclicGroup(2)")
    sage: A = BurnsideRing(C2)
    sage: A
    Burnside ring of C2

The basis is indexed by conjugacy classes of subgroups.  The displayed names
are transitive orbits ``G/H``::

    sage: A.basis_names()
    ('[C2/1]', '[C2/C2]')
    sage: A.gens()
    ([C2/1], [C2/C2])

Creating elements
-----------------

Elements can be built with ordinary Sage arithmetic::

    sage: x = 3*A.gen(0) - 4*A.one()
    sage: x
    3*[C2/1] - 4*[C2/C2]

They can also be built from coefficient vectors in the orbit basis::

    sage: A.from_coefficients([3, -4])
    3*[C2/1] - 4*[C2/C2]

The coefficient vector and the mark vector are both available::

    sage: x.coefficients()
    (3, -4)
    sage: x.marks()
    (2, -4)

The table of marks can be inspected directly::

    sage: A.table_of_marks()
    [2 0]
    [1 1]

Restriction, transfer, and norm
-------------------------------

Use GAP subgroups for restriction.  For example, let ``e`` be the trivial
subgroup of ``C2``::

    sage: e = C2.TrivialSubgroup()
    sage: A.gen(0).restrict(e)
    2*[1/1]

The short alias ``res`` is also available.

Transfers and norms go from a subgroup to a larger group::

    sage: Ae = BurnsideRing(e)
    sage: Ae.one().transfer(C2)
    [C2/1]
    sage: Ae.one().norm(C2)
    [C2/C2]

The aliases ``tr`` and ``nm`` are also available.

Matrices
--------

Restriction and transfer matrices can be returned in the orbit-basis
coordinates or in mark coordinates::

    sage: A.restriction_matrix(e)
    [2]
    [1]
    sage: A.restriction_matrix(e, basis="marks")
    [1]
    [0]

Power operations
----------------

The library includes symmetric powers, exterior powers, and Siebeneicher
powers of virtual Burnside ring elements::

    sage: v = A.from_coefficients([-1, 2])
    sage: v.symmetric_power(2)
    -[C2/1] + 2*[C2/C2]
    sage: v.exterior_power(2)
    0
    sage: v.siebeneicher_power(2)
    -[C2/1] + 2*[C2/C2]

Pickling
--------

If a ring needs to be pickled, construct it from a GAP expression string::

    sage: B = BurnsideRing.from_gap("CyclicGroup(2)")

Rings built from live GAP objects still work normally, but they cannot be
pickled because the live GAP object itself cannot be reconstructed later.
