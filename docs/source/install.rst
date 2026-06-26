Installation
============

For students
------------

There is no package installation step.

1. Download ``burnsidelib.py`` from the repository.
2. Put it in the same folder as your Sage notebook or ``.sage`` file.
3. Add this near the top of the worksheet or script::

       load("burnsidelib.py")

After that, the class ``BurnsideRing`` is available in the worksheet.

In a Sage script, Python-style imports also work when ``burnsidelib.py`` is in
the same directory::

    from burnsidelib import *

Requirements
------------

``burnsidelib`` expects to run inside SageMath.  It uses Sage's interfaces to
GAP, matrices, vectors, power series, and combinatorial free modules.

Building these docs locally
---------------------------

From the repository root, run::

    sage -pip install -r docs/requirements.txt
    sage -python -m sphinx -b html docs/source docs/_build/html

Then open ``docs/_build/html/index.html`` in a browser.
