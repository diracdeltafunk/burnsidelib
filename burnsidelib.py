from collections import Counter

from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.interfaces.gap import gap
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ


def _rebuild_burnside_ring_from_gap(gap_expression):
    return BurnsideRing.from_gap(gap_expression)


def _gap_bool(value):
    try:
        return bool(value.sage())
    except AttributeError:
        return bool(value)


def _integer_vector(values, length, name):
    try:
        entries = list(values)
    except TypeError as err:
        raise TypeError(f"{name} must be an iterable of integers") from err
    if len(entries) != length:
        raise ValueError(f"{name} must have length {length}; got length {len(entries)}")
    try:
        return vector(ZZ, entries)
    except (TypeError, ValueError) as err:
        raise ValueError(f"{name} must contain only integers") from err


def _nonnegative_integer(k, name):
    try:
        k = ZZ(k)
    except (TypeError, ValueError) as err:
        raise ValueError(f"{name} must be a nonnegative integer") from err
    if k < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return k


class BurnsideRingElement(CombinatorialFreeModule.Element):
    def coefficients(self):
        parent = self.parent()
        coeffs = self.monomial_coefficients()
        return vector(
            ZZ, [coeffs.get(i, ZZ.zero()) for i in range(parent._num_cc_subgroups)]
        )

    def coeffs(self):
        return self.coefficients()

    def list(self):
        return list(self.coefficients())

    def marks(self):
        return self.coefficients() * self.parent()._tommat

    def restrict(self, H):
        A_G = self.parent()
        A_H = BurnsideRing(H)
        return A_H.from_coefficients(self.coefficients() * A_G.restriction_matrix(H))

    def res(self, H):
        return self.restrict(H)

    def transfer(self, G):
        A_H = self.parent()
        A_G = BurnsideRing(G)
        return A_G.from_coefficients(self.coefficients() * A_H.transfer_matrix(G))

    def tr(self, G):
        return self.transfer(G)

    def norm(self, G):
        A_H = self.parent()
        A_G = BurnsideRing(G)
        A_H._check_subgroup_of(G)

        H = A_H.group()
        marks_in = self.marks()
        marks_out = vector(ZZ, [1] * A_G._num_cc_subgroups)

        for i, K in enumerate(A_G._cc_reps):
            for g, _ in G.DoubleCosetRepsAndSizes(K, H):
                subgroup = gap.Intersection(H, K.ConjugateGroup(g))
                j = A_H._subgroup_index(subgroup)
                marks_out[i] *= marks_in[j]
        return A_G.from_marks(marks_out)

    def nm(self, G):
        return self.norm(G)

    def symmetric_power(self, k):
        return self._power_operation(k, "symmetric")

    def exterior_power(self, k):
        return self._power_operation(k, "exterior")

    def siebeneicher_power(self, k):
        return self._power_operation(k, "siebeneicher")

    def _power_operation(self, k, operation):
        k = _nonnegative_integer(k, "k")
        A_G = self.parent()
        R = PowerSeriesRing(QQ, "t", default_prec=int(k) + 1)
        t = R.gen()
        marks_out = vector(QQ, A_G._num_cc_subgroups)

        for i, H in enumerate(A_G._cc_reps):
            restricted = self.restrict(H)
            A_H = restricted.parent()
            series = R.one()
            for j, coeff in enumerate(restricted.coefficients()):
                orbit_size = ZZ(A_H._tommat[j, 0])
                if coeff == 0:
                    continue
                if operation == "symmetric":
                    series *= (1 - t**orbit_size) ** (-coeff)
                elif operation == "exterior":
                    series *= (1 + t**orbit_size) ** coeff
                elif operation == "siebeneicher":
                    series *= (1 - (-t) ** orbit_size) ** coeff
                elif operation == "dual exterior":
                    series *= (1 + (-t) ** orbit_size) ** (-coeff)
                else:
                    raise ValueError(f"unknown power operation: {operation}")
            marks_out[i] = series.padded_list(int(k) + 1)[int(k)]

        return A_G.from_marks(marks_out)


class BurnsideRing(CombinatorialFreeModule):
    Element = BurnsideRingElement

    @staticmethod
    def __classcall_private__(cls, group, names=None, gap_expression=None):
        if gap_expression is not None:
            gap_expression = str(gap_expression)
            return super().__classcall__(
                cls, gap_expression, names=None, gap_expression=gap_expression
            )
        return super().__classcall__(cls, group, names=None, gap_expression=None)

    def __init__(self, group, names=None, gap_expression=None):
        if gap_expression is not None:
            group = gap(gap_expression)
        self._group = group
        self._gap_expression = gap_expression
        self._tom = group.TableOfMarks()
        self._tommat = matrix(ZZ, self._tom.MatTom())
        self._tommat_inv = self._tommat.inverse()
        self._num_cc_subgroups = self._tommat.nrows()
        self._index_of_unit = next(
            i for i in range(self._num_cc_subgroups) if self._tommat[i, 0] == 1
        )
        self._cc_reps = tuple(
            self._tom.RepresentativeTom(n + 1) for n in range(self._num_cc_subgroups)
        )
        self._basis_names = self._make_basis_names()
        self._multiplication_table = self._make_multiplication_table()

        super().__init__(
            ZZ,
            range(self._num_cc_subgroups),
            element_class=BurnsideRingElement,
            category=AlgebrasWithBasis(ZZ).Commutative(),
            prefix="b",
        )

    @classmethod
    def from_gap(cls, gap_expression):
        return cls(gap_expression, gap_expression=gap_expression)

    def __reduce__(self):
        if self._gap_expression is None:
            raise TypeError(
                "BurnsideRing objects built from live GAP groups are not pickleable; "
                "construct with BurnsideRing.from_gap(...) to make pickling possible"
            )
        return (_rebuild_burnside_ring_from_gap, (self._gap_expression,))

    def _make_basis_names(self):
        group_name = str(self._group.StructureDescription())
        subgroup_names = [str(H.StructureDescription()) for H in self._cc_reps]
        counts = Counter(subgroup_names)
        names = []
        for i, subgroup_name in enumerate(subgroup_names):
            if counts[subgroup_name] == 1:
                names.append(f"[{group_name}/{subgroup_name}]")
            else:
                names.append(f"[{i}: {group_name}/{subgroup_name}]")
        return tuple(names)

    def _make_multiplication_table(self):
        table = {}
        for i in range(self._num_cc_subgroups):
            for j in range(i, self._num_cc_subgroups):
                product_marks = vector(
                    ZZ,
                    [
                        self._tommat[i, k] * self._tommat[j, k]
                        for k in range(self._num_cc_subgroups)
                    ],
                )
                table[(i, j)] = self._coefficients_from_marks(product_marks)
        return table

    def _coefficients_from_marks(self, marks):
        coeffs = vector(QQ, marks) * self._tommat_inv
        try:
            return vector(ZZ, coeffs)
        except (TypeError, ValueError) as err:
            raise ValueError(
                "the given marks are not marks of a Burnside ring element for this group"
            ) from err

    def _subgroup_index(self, subgroup):
        for i, representative in enumerate(self._cc_reps):
            if _gap_bool(self._group.IsConjugate(representative, subgroup)):
                return i
        raise ValueError(
            "subgroup is not conjugate to any table-of-marks representative"
        )

    def _check_subgroup_of(self, G):
        if not _gap_bool(gap.IsSubgroup(G, self._group)):
            raise ValueError("the source group must be a subgroup of the target group")

    def _check_contains_subgroup(self, H):
        if not _gap_bool(gap.IsSubgroup(self._group, H)):
            raise ValueError("H must be a subgroup of the group of this Burnside ring")

    def _repr_(self):
        return f"Burnside ring of {self._group.StructureDescription()}"

    def _repr_short(self):
        return f"A({self._group.StructureDescription()})"

    def _repr_term(self, i):
        return self._basis_names[i]

    def one_basis(self):
        return self._index_of_unit

    def product_on_basis(self, i, j):
        if j < i:
            i, j = j, i
        return self.from_coefficients(self._multiplication_table[(i, j)])

    def group(self):
        return self._group

    def table_of_marks(self):
        return matrix(ZZ, self._tommat)

    def subgroup_representatives(self):
        return tuple(self._cc_reps)

    def basis_names(self):
        return tuple(self._basis_names)

    def gen_names(self):
        return self.basis_names()

    def orbit_basis(self):
        return tuple(self.monomial(i) for i in range(self._num_cc_subgroups))

    def gens(self):
        return self.orbit_basis()

    def gen(self, i=0):
        return self.orbit_basis()[i]

    def ngens(self):
        return self._num_cc_subgroups

    def _first_ngens(self, n):
        return self.orbit_basis()[:n]

    def from_coefficients(self, v):
        coeffs = _integer_vector(v, self._num_cc_subgroups, "coefficients")
        return self._from_dict(
            {i: coeff for i, coeff in enumerate(coeffs) if coeff != 0}, coerce=False
        )

    def from_vec(self, v):
        return self.from_coefficients(v)

    def from_marks(self, v):
        marks = _integer_vector(v, self._num_cc_subgroups, "marks")
        return self.from_coefficients(self._coefficients_from_marks(marks))

    def restriction_matrix(self, H, basis="coefficients"):
        self._check_contains_subgroup(H)
        A_H = BurnsideRing(H)
        marks_matrix = matrix(
            ZZ,
            self._num_cc_subgroups,
            A_H._num_cc_subgroups,
            lambda i, j: ZZ(
                _gap_bool(self._group.IsConjugate(self._cc_reps[i], A_H._cc_reps[j]))
            ),
        )
        if basis == "marks":
            return marks_matrix
        if basis == "coefficients":
            return matrix(ZZ, self._tommat * marks_matrix * A_H._tommat_inv)
        raise ValueError('basis must be either "coefficients" or "marks"')

    def transfer_matrix(self, G, basis="coefficients"):
        self._check_subgroup_of(G)
        A_G = BurnsideRing(G)
        coefficient_matrix = matrix(
            ZZ,
            self._num_cc_subgroups,
            A_G._num_cc_subgroups,
            lambda i, j: ZZ(
                _gap_bool(G.IsConjugate(self._cc_reps[i], A_G._cc_reps[j]))
            ),
        )
        if basis == "coefficients":
            return coefficient_matrix
        if basis == "marks":
            return matrix(ZZ, self._tommat_inv * coefficient_matrix * A_G._tommat)
        raise ValueError('basis must be either "coefficients" or "marks"')
