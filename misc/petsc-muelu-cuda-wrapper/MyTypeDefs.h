#ifndef MYTYPEDEFS_H
#define MYTYPEDEFS_H
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
typedef double                                            scalar_type;
typedef Tpetra::CrsMatrix<scalar_type>                    matrix_type;
typedef Tpetra::MultiVector<scalar_type>                  vec_type;
typedef Tpetra::Map<>             map_type;
#endif
