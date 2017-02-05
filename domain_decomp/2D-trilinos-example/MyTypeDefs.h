#ifndef MYTYPEDEFS_H
#define MYTYPEDEFS_H
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
typedef Tpetra::CrsMatrix<>      matrix_type;
typedef Tpetra::MultiVector<>    vector_type;
typedef Tpetra::MultiVector<int> int_vector_type;
typedef Kokkos::View<double **, Kokkos::Serial> view_type;
typedef Tpetra::Map<> map_type;
#endif
