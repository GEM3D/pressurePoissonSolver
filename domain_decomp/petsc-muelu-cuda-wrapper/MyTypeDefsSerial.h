#ifndef MYTYPEDEFS_H
#define MYTYPEDEFS_H
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
typedef Kokkos::Compat::KokkosSerialWrapperNode node_type;
typedef double                                            scalar_type;
typedef Tpetra::CrsMatrix<scalar_type,int,int,node_type
>                    matrix_type;
typedef Tpetra::MultiVector<scalar_type,int,int,node_type
>                  vec_type;
typedef Tpetra::Map<int,int,node_type
>             map_type;
typedef Tpetra::Operator<scalar_type,int,int,node_type
> op_type;
#endif
