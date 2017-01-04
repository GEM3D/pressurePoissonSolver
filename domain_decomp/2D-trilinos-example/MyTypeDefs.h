#ifndef MYTYPEDEFS_H
#define MYTYPEDEFS_H
#include <Amesos2.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
typedef Epetra_CrsMatrix   matrix_type;
typedef Epetra_MultiVector vector_type;
typedef Amesos2::Solver<matrix_type, vector_type> solver_type;
typedef Epetra_Map map_type;
#endif
