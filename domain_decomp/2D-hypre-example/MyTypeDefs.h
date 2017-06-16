#ifndef MYTYPEDEFS_H
#define MYTYPEDEFS_H
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix_Helpers.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
typedef double scalar_type;
typedef Tpetra::CrsMatrix<scalar_type>                    matrix_type;
typedef Tpetra::Experimental::BlockCrsMatrix<scalar_type> block_matrix_type;
typedef Tpetra::MultiVector<scalar_type>                  vector_type;
typedef Tpetra::Vector<scalar_type>                       single_vector_type;
typedef Tpetra::MultiVector<int>                          int_vector_type;
typedef Kokkos::View<double **, Kokkos::Serial> view_type;
typedef Tpetra::Map<> map_type;
typedef Tpetra::CrsGraph<> graph_type;
#endif
