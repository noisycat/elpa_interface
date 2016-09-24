/* Charles Lena - 2016 */
#ifndef __EXTERN_UTILITY_H__
#define __EXTERN_UTILITY_H__
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif 
	int section_matrix_(int* na, double* a, int* na_rows, int* nblk, int* my_prow, int* my_pcol, int* np_rows, int* np_cols, MPI_Fint* my_mpi_comm_world, double* input_matrix);
	int gather_matrix_(int* na, double* z, int* na_rows, int* na_cols, int* nblk, int* my_prow, int* my_pcol, int* np_rows, int* np_cols, MPI_Fint* my_mpi_comm_world, MPI_Fint* my_blacs_ctxt, double* input_matrix);

#ifdef __cplusplus
}
#endif 
#endif
