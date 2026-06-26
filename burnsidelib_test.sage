from burnsidelib import *

import pickle


def assert_raises(exc_type, func):
    try:
        func()
    except exc_type:
        return
    except Exception as err:
        raise AssertionError(f"expected {exc_type.__name__}, got {type(err).__name__}") from err
    raise AssertionError(f"expected {exc_type.__name__}")


def assert_pointwise_marks_multiply(A):
    for x in A.gens():
        for y in A.gens():
            assert (x * y).marks() == vector(ZZ, [a * b for a, b in zip(x.marks(), y.marks())])



# Compatibility with Sage generator shorthand from the original smoke test.
A.<t> = BurnsideRing(gap("CyclicGroup(2)"))
assert (3 * t - 4).marks() == vector(ZZ, [2, -4])


# Core C2 behavior.
C2 = gap("CyclicGroup(2)")
e = C2.TrivialSubgroup()
A = BurnsideRing(C2)
Ae = BurnsideRing(e)

assert str(A.zero()) == "0"
assert A.one() == A.gen(1)
assert A.gens() == A.orbit_basis()
assert A.ngens() == 2
assert A.basis_names() == ("[C2/1]", "[C2/C2]")
assert A.gen_names() == A.basis_names()
assert str(A.subgroup_representatives()[0].StructureDescription()) == "1"

x = A.from_coefficients([3, -4])
assert x.coefficients() == vector(ZZ, [3, -4])
assert x.coeffs() == vector(ZZ, [3, -4])
assert x.list() == [3, -4]
assert x.marks() == vector(ZZ, [2, -4])
assert A.from_vec([3, -4]) == x
assert A.from_marks([2, -4]) == x
assert A.gen(0) * A.gen(0) == A.from_coefficients([2, 0])
assert_pointwise_marks_multiply(A)

coeffs = x.coefficients()
coeffs[0] = 999
assert x.coefficients() == vector(ZZ, [3, -4])

tom = A.table_of_marks()
tom[0, 0] = 999
assert A.table_of_marks()[0, 0] == 2

assert_raises(ValueError, lambda: A.from_coefficients([1]))
assert_raises(ValueError, lambda: A.from_coefficients([1, 2, 3]))
assert_raises(ValueError, lambda: A.from_marks([1, 0]))
assert_raises(ValueError, lambda: A.restriction_matrix(e, basis="ghost"))
assert_raises(ValueError, lambda: A.from_coefficients([1 / 2, 0]))


# Restriction, transfer, and norm for an actual subgroup inclusion e <= C2.
assert A.gen(0).restrict(e) == Ae.from_coefficients([2])
assert A.one().res(e) == Ae.one()
assert Ae.one().transfer(C2) == A.gen(0)
assert Ae.one().tr(C2) == A.gen(0)
assert Ae.from_coefficients([4]).norm(C2) == A.from_coefficients([6, 4])
assert Ae.from_coefficients([4]).nm(C2) == A.from_coefficients([6, 4])
assert Ae.from_coefficients([-1]).norm(C2) == A.from_coefficients([1, -1])
assert Ae.zero().norm(C2) == A.zero()
assert A.restriction_matrix(e) == matrix(ZZ, [[2], [1]])
assert A.restriction_matrix(e, basis="marks") == matrix(ZZ, [[1], [0]])
assert Ae.transfer_matrix(C2) == matrix(ZZ, [[1, 0]])
assert Ae.transfer_matrix(C2, basis="marks") == matrix(ZZ, [[2, 0]])


# A slightly larger group: S3 restricted to its standard S2 subgroup.
S3 = gap("SymmetricGroup(3)")
S2 = gap("SymmetricGroup(2)")
AS3 = BurnsideRing(S3)
res_sum = sum(AS3.gens()).restrict(S2)
assert res_sum.parent() is BurnsideRing(S2)
assert res_sum.coefficients() == vector(ZZ, [5, 2])
assert res_sum.marks() == vector(ZZ, [12, 2])
assert_pointwise_marks_multiply(AS3)


# Duplicate abstract subgroup names are disambiguated by basis index.
Q8 = gap("QuaternionGroup(8)")
AQ8 = BurnsideRing(Q8)
assert len(set(AQ8.basis_names())) == len(AQ8.basis_names())
assert "[2: Q8/C4]" in AQ8.basis_names()
assert "[3: Q8/C4]" in AQ8.basis_names()
assert "[4: Q8/C4]" in AQ8.basis_names()


# Power operations on virtual Burnside ring elements.
v = A.from_coefficients([-1, 2])
assert v.symmetric_power(0) == A.one()
assert v.exterior_power(0) == A.one()
assert v.siebeneicher_power(0) == A.one()
assert v.symmetric_power(1) == v
assert v.exterior_power(1) == v
assert v.siebeneicher_power(1) == v
assert v.symmetric_power(2) == v
assert v.exterior_power(2) == A.zero()
assert v.siebeneicher_power(2) == v
assert_raises(ValueError, lambda: v.symmetric_power(-1))

for k in range(1, 5):
    total = A.zero()
    for i in range(k + 1):
        total += (-1) ** i * v.siebeneicher_power(i) * v.symmetric_power(k - i)
    assert total == A.zero()


# Pickling is supported for rings rebuilt from a GAP expression.
B = BurnsideRing.from_gap("CyclicGroup(2)")
assert B is BurnsideRing.from_gap("CyclicGroup(2)")
assert pickle.loads(pickle.dumps(B)) is B
assert pickle.loads(pickle.dumps(B.gen(0))) == B.gen(0)
TestSuite(B).run()

print("burnsidelib tests passed")
