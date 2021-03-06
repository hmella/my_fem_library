// Copyright (C) 2017 Chris Richardson and Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <petscksp.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscdm.h>

// pybind11 casters for PETSc/petsc4py objects

#include <petsc4py/petsc4py.h>

// Import petsc4py on demand
#define VERIFY_PETSC4PY(func)   \
  if (!func)                    \
  {                             \
    if (import_petsc4py() != 0) \
    {                           \
      std::cout << "ERROR: could not import petsc4py!" << std::endl; \
      throw std::runtime_error("Error when importing petsc4py");     \
    }                           \
  }

// Macro for casting between dolfin and petsc4py objects
#define PETSC_CASTER_MACRO(TYPE, NAME)          \
  template <> class type_caster<_p_##TYPE>      \
    {                                           \
    public:                                     \
      PYBIND11_TYPE_CASTER(TYPE, _(#NAME));     \
      bool load(handle src, bool)               \
      {                                         \
        VERIFY_PETSC4PY(PyPetsc##TYPE##_Get);   \
        if (PyObject_TypeCheck(src.ptr(), &PyPetsc##TYPE##_Type) == 0)  \
          return false;                                                 \
        value = PyPetsc##TYPE##_Get(src.ptr());                         \
        return true;                                                    \
      }                                                                 \
                                                                        \
      static handle cast(TYPE src, pybind11::return_value_policy policy, handle parent) \
      {                                                                 \
        VERIFY_PETSC4PY(PyPetsc##TYPE##_New);                           \
        return pybind11::handle(PyPetsc##TYPE##_New(src));              \
      }                                                                 \
                                                                        \
      operator TYPE()                                                   \
      { return value; }                                                 \
    }



namespace pybind11
{
  namespace detail
  {
    PETSC_CASTER_MACRO(Vec, vec);
    PETSC_CASTER_MACRO(Mat, mat);
    PETSC_CASTER_MACRO(KSP, ksp);
    PETSC_CASTER_MACRO(DM, DM);
  }
}

#undef PETSC_CASTER_MACRO
