# burnsidelib

`burnsidelib` is a small single-file SageMath library for computations with
Burnside rings of finite groups.

## Installation

There is no package to install.

1. Put `burnsidelib.py` in the same folder as your Sage notebook or `.sage`
   file.
2. Add this near the top:

   ```python
   load("burnsidelib.py")
   ```

Then create Burnside rings from GAP groups:

```python
C2 = gap("CyclicGroup(2)")
A = BurnsideRing(C2)
x = 3*A.gen(0) - 4*A.one()
x.marks()
```

## Documentation

The documentation source lives in `docs/source`.  To build it locally:

```sh
sage -pip install -r docs/requirements.txt
sage -python -m sphinx -b html docs/source docs/_build/html
```

The GitHub Actions workflow in `.github/workflows/docs.yml` builds the docs
with Sage and deploys them to GitHub Pages on pushes to `master` or `main`.
In the repository settings, set Pages to deploy from GitHub Actions.

## Tests

Run the smoke tests with:

```sh
sage burnsidelib_test.sage
```
