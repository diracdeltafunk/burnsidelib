API reference
=============

.. currentmodule:: burnsidelib

Burnside rings
--------------

.. autoclass:: BurnsideRing
   :members:
      from_gap,
      one_basis,
      product_on_basis,
      group,
      table_of_marks,
      subgroup_representatives,
      basis_names,
      gen_names,
      orbit_basis,
      gens,
      gen,
      ngens,
      from_coefficients,
      from_vec,
      from_marks,
      restriction_matrix,
      transfer_matrix

Elements
--------

.. autoclass:: BurnsideRingElement
   :members:
      coefficients,
      coeffs,
      list,
      marks,
      restrict,
      res,
      transfer,
      tr,
      norm,
      nm,
      symmetric_power,
      exterior_power,
      dual_exterior_power,
      siebeneicher_power,
      adams_operation,
      adams_operation_by_marks,
      symmetric_adams_operation,
      exterior_adams_operation,
      dual_exterior_adams_operation,
      siebeneicher_adams_operation
