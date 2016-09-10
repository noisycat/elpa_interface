/* This file describes a somewhat ghetto ELPA interface from C++ */
#ifndef __ELPA_INTERFACE_HPP__
#define __ELPA_INTERFACE_HPP__

#include <vector>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <iomanip>

#include "mycomplex.h"

#ifdef DEBUG
extern "C" { void plus_(double * A, int * N, double* B); }
#endif
extern "C" {
	/* raw copy and past from elpa header because of C99 / C++ disagreement over complex stypes in GNU */
	int elpa_get_communicators(int mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);
	int get_elpa_communicators(int mpi_comm_world, int my_prow, int my_pcol, int *mpi_comm_rows, int *mpi_comm_cols);
	int elpa_solve_evp_real_1stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols);
	int elpa_solve_evp_complex_1stage(int na, int nev, double _Complex *a, int lda, double *ev, double _Complex *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols);
	int elpa_solve_evp_real_2stage(int na, int nev, double *a, int lda, double *ev, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int mpi_comm_all, int THIS_REAL_ELPA_KERNEL_API, int useQR);
	int elpa_solve_evp_complex_2stage(int na, int nev, double _Complex *a, int lda, double *ev, double _Complex *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int mpi_comm_all, int THIS_COMPLEX_ELPA_KERNEL_API);
	int elpa_solve_tridi(int na, int nev, double *d, double *e, double *q, int ldq, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
	int elpa_mult_at_b_real(char uplo_a, char uplo_c, int na, int ncb, double *a, int lda, double *b, int ldb, int nlbk, int mpi_comm_rows, int mpi_comm_cols, double *c, int ldc);
	int elpa_mult_ah_b_complex(char uplo_a, char uplo_c, int na, int ncb, double _Complex *a, int lda, double _Complex *b, int ldb, int nblk, int mpi_comm_rows, int mpi_comm_cols, double _Complex *c, int ldc);
	int elpa_invert_trm_real(int na, double *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
	int elpa_invert_trm_complex(int na, double _Complex *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
	int elpa_cholesky_real(int na, double *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);
	int elpa_cholesky_complex(int na, double _Complex *a, int lda, int nblk, int matrixCols, int mpi_comm_rows, int mpi_comm_cols, int wantDebug);


	/* A reminder - all of these are the recorded FORTRAN interfaces that emulate Scalapack interfaces
	   so name mangling exists and is a very real thing */
	void get_elpa_row_col_comms_(MPI_Fint * MPI_COMM_GLOBAL, int * MY_PROW, int * MY_PCOL, MPI_Fint * MPI_COMM_ROWS, MPI_Fint * MPI_COMM_COLS);
	/* ELPA1 solvers */
	void solve_evp_real_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	/* ELPA2 solvers */
	void solve_evp_real_2stage_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_2stage_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);

	// pretty sure this works?
	int numroc_(int, int, int, int);

	//
	void solve_full_(MPI_Fint *, double *, int* , double *, int*);
}

using std::vector;

template<float> class ELPA_Interface;
template<std::complex> class ELPA_Interface;

template<typename T> class ELPA_Interface {
	// Keep in mind, dealing with A * v = B * w * v
	/* All ELPA use comes in essentially 3 steps:
	 * - Set up the BLACS distribution
	 * - Create the ELPA specific communicators for rows and columns
	 * - Call the solver for ELPA1 or ELPA2
	 */
	private:
		int N, M; // actual matrix dimension
		int nblk; // nblk: Blocking factor in block cyclic distribution
		int na; // na: size of system
		int nev; // nev: number of eigenvalues to calculate
		int mypid, nprocs;
	public:
		// Constructor 
		ELPA_Interface() {};
		void StartBLACS() {};
		void StopBLACS() {};
		void Transfer(vector< vector<T> > &A, T* a, int nblk = 16){};
		void Broadcast(vector< vector<T> > &A, T* a, int nblk = 16){};
		void Gather(vector< vector<T> > &A, T* a, int nblk = 16){};

		// Associate - transfer matrix
		void Associate(vector< vector<T> > &A, double* a, int nblk = 16) 
		{ 
			this->N = A.size();
			this->M = A[0].size();
			std::cout << this->N << ' ' << this->M << '\n';
			assert(N == M);
			//-------------------------------------------------------------------------------
			// Please set system size parameters below!
			// na:   System size
			// nev:  Number of eigenvectors to be calculated
			// nblk: Blocking factor in block cyclic distribution
			//-------------------------------------------------------------------------------
			this->nblk = nblk;
			this->na = this->N;
			this->nev = this->N;

			for(int i = 0; i < N; i++) {
				for(int j = 0; j < M; j++) {
					a[M*i+j] = (A[i])[j];
				}
			}
		};
		// New interface
		void Solve(double* A, int N_, MPI_Comm* the_comm, double* eigvals_, int nblk = 16) {
			int N = N_;
			int myid;
			MPI_Fint the_comm_f = MPI_Comm_c2f(*the_comm);
			//fprintf(stderr,"%d %p %d %p %d\n", the_comm_f, A, N, eigvals_, nblk );
			solve_full_(&the_comm_f, A, &N, eigvals_, &nblk);
#if defined(SOLVE_PRINT_EXPLODE)
			MPI_Comm_rank(MPI_COMM_WORLD,&myid);
			if(myid==0) {
				for(int i = 0; i < N; i++) {
					for(int j = 0; j < N; j++) {
						printf(" %e",A[i*N + j]);
					}
					printf("\n");
				}
			}
			for(int i = 0; i < N; i++) {
				printf("%d %e\n",i,eigvals_[i]);
			}
#endif
		};

		// Destructor - ELPA_Interface object goes byebye
		~ELPA_Interface() {  };
#if defined(INTERFACE_UNNECSSARY_CRUFT)
		int setSolve(MPI_Comm * mpi_comm_rows, MPI_Comm * mpi_comm_cols,MPI_Comm * the_comm){
			this->ew.resize(na);
			int success = solve_evp_real_2stage_(this->na, this->nev, this->a.data(), this->na_rows, this->ev.data(), this->ew.data(), this->na_rows, this->nblk, this->mpi_comm_rows, this->mpi_comm_cols, this->the_comm);
			return success;
		};

		// Estimate - dumps to stdout an estimated time of execution for a matrix, optionally kills if above limit
		double Estimate() { return pow(sqrt(this->A.size()),3);};
		double Estimate(double timelimit) { };

		// Create ELPA specific layout - required to solve
		void CreateELPAComms(MPI_Comm* the_comm, int* my_prow,int* my_pcol,MPI_Comm* mpi_comm_rows, MPI_Comm* mpi_comm_cols) { 
			// All ELPA routines need MPI communicators for communicating within
			// rows or columns of processes, these are set in get_elpa_row_col_comms.
			MPI_Fint my_blacs_ctxt = MPI_Comm_c2f(*the_comm);
			MPI_Fint mpi_comm_rows_f;
			MPI_Fint mpi_comm_cols_f;

			get_elpa_row_col_comms_(&my_blacs_ctxt, my_prow, my_pcol, &mpi_comm_rows_f, &mpi_comm_cols_f);
			*mpi_comm_rows = MPI_Comm_f2c(mpi_comm_rows_f);
			*mpi_comm_cols = MPI_Comm_f2c(mpi_comm_cols_f);
		} ;

		// Regather
		void Gather() {};

		// Solve from start to finish
		void Solve2(vector< vector<T> > &A, MPI_Comm* the_comm, int myid, int nprocs, int required_mpi_thread_level, int provided_mpi_thread_level){};
		double Residual(int myid, double err){ };
		double EigenOrth(vector<T> tmp, double* z, double & err, int myid){};
#endif
};
#endif
