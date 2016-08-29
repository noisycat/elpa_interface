/* This file describes a somewhat ghetto ELPA interface from C++ */
#ifndef __ELPA_INTERFACE_HPP__
#define __ELPA_INTERFACE_HPP__

#include <vector>
#include "mpi.h"
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
	#include <elpa/elpa.h>

	/* A reminder - all of these are the recorded FORTRAN interfaces that emulate Scalapack interfaces
	   so name mangling exists and is a very real thing */
	void get_elpa_row_col_comms_(MPI_Fint * MPI_COMM_GLOBAL, int * MY_PROW, int * MY_PCOL, MPI_Fint * MPI_COMM_ROWS, MPI_Fint * MPI_COMM_COLS);
	/* ELPA1 solvers */
	void solve_evp_real_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	/* ELPA2 solvers */
	void solve_evp_real_2stage_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_2stage_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);

	// haven't checked these
	void Cblacs_gridinit(int* icontxt, char *order, int* nprow, int *npcol);
	void Cblacs_gridinfo(int* icontxt, int *nprow, int *npcol, int *myprow, int *mypcol);
	void Cblacs_gridexit(int* icontxt);

	// pretty sure this works?
	int numroc_(int, int, int, int);

	//
	void solve_full_(MPI_Fint *, double *, int* , double *, int*);
}

using std::vector;

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
	public:
		// Constructor 
		ELPA_Interface() { };

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
		void newSolve(){
			int argc, char** argv;
			int myid;
			int nprocs;
#ifndef WITH_MPI
			int MPI_COMM_WORLD;
#endif
			int na, nev, nblk;

			int status;

			int np_cols, np_rows, np_colsStart;

			int my_blacs_ctxt, nprow, npcol, my_prow, my_pcol;

			int mpierr;

			int my_mpi_comm_world;
			int mpi_comm_rows, mpi_comm_cols;

			int info, *sc_desc;

			int na_rows, na_cols;
			double startVal;

			complex double *a, *z, *as, *tmp1, *tmp2;

			double *ev, *xr;

			int *iseed;

			int success;

			int THIS_COMPLEX_ELPA_KERNEL_API;
#ifdef WITH_MPI
			MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
			MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#else
			nprocs = 1;
			myid =0;
			MPI_COMM_WORLD=1;
#endif
			na = 1000;
			nev = 500;
			nblk = 16;

			if (myid == 0) {
				printf("This is the c version of an ELPA test-programm\n");
				printf("\n");
				printf("It will call the 1stage ELPA complex solver for a matrix\n");
				printf("of matrix size %d. It will compute %d eigenvalues\n",na,nev);
				printf("and uses a blocksize of %d\n",nblk);
				printf("\n");
				printf("This is an example program with much less functionality\n");
				printf("as it's Fortran counterpart. It's only purpose is to show how \n");
				printf("to evoke ELPA1 from a c programm\n");
				printf("\n");
			}

			status = 0;

			startVal = sqrt((double) nprocs);
			np_colsStart = (int) round(startVal);
			for (np_cols=np_colsStart;np_cols>1;np_cols--){
				if (nprocs %np_cols ==0){
					break;
				}
			}

			np_rows = nprocs/np_cols;

			if (myid == 0) {
				printf("\n");
				printf("Number of processor rows %d, cols %d, total %d \n",np_rows,np_cols,nprocs);
			}

			/* set up blacs */
			/* convert communicators before */
#ifdef WITH_MPI
			my_mpi_comm_world = MPI_Comm_c2f(MPI_COMM_WORLD);
#else
			my_mpi_comm_world = 1;
#endif
			set_up_blacsgrid_from_fortran(my_mpi_comm_world, &my_blacs_ctxt, &np_rows, &np_cols, &nprow, &npcol, &my_prow, &my_pcol);

			if (myid == 0) {
				printf("\n");
				printf("Past BLACS_Gridinfo...\n");
				printf("\n");
			}

			/* get the ELPA row and col communicators. */
			/* These are NOT usable in C without calling the MPI_Comm_f2c function on them !! */
#ifdef WITH_MPI
			my_mpi_comm_world = MPI_Comm_c2f(MPI_COMM_WORLD);
#endif
			mpierr = get_elpa_communicators(my_mpi_comm_world, my_prow, my_pcol, &mpi_comm_rows, &mpi_comm_cols);

			if (myid == 0) {
				printf("\n");
				printf("Past split communicator setup for rows and columns...\n");
				printf("\n");
			}

			sc_desc = malloc(9*sizeof(int));

			set_up_blacs_descriptor_from_fortran(na, nblk, my_prow, my_pcol, np_rows, np_cols, &na_rows, &na_cols, sc_desc, my_blacs_ctxt, &info);

			if (myid == 0) {
				printf("\n");
				printf("Past scalapack descriptor setup...\n");
				printf("\n");
			}

			/* allocate the matrices needed for elpa */
			if (myid == 0) {
				printf("\n");
				printf("Allocating matrices with na_rows=%d and na_cols=%d\n",na_rows, na_cols);
				printf("\n");
			}

			a  = malloc(na_rows*na_cols*sizeof(complex double));
			z  = malloc(na_rows*na_cols*sizeof(complex double));
			as = malloc(na_rows*na_cols*sizeof(complex double));

			xr = malloc(na_rows*na_cols*sizeof(double));


			ev = malloc(na*sizeof(double));

			tmp1  = malloc(na_rows*na_cols*sizeof(complex double));
			tmp2 = malloc(na_rows*na_cols*sizeof(complex double));

			iseed = malloc(4096*sizeof(int));

			prepare_matrix_complex_from_fortran(na, myid, na_rows, na_cols, sc_desc, iseed, xr, a, z, as);

			free(xr);

			if (myid == 0) {
				printf("\n");
				printf("Entering ELPA 2stage complex solver\n");
				printf("\n");
			}
#ifdef WITH_MPI
			mpierr = MPI_Barrier(MPI_COMM_WORLD);
#endif
			THIS_COMPLEX_ELPA_KERNEL_API = ELPA2_COMPLEX_KERNEL_GENERIC;
			success = elpa_solve_evp_complex_2stage(na, nev, a, na_rows, ev, z, na_rows, nblk, na_cols, mpi_comm_rows, mpi_comm_cols, my_mpi_comm_world, THIS_COMPLEX_ELPA_KERNEL_API);

			if (success != 1) {
				printf("error in ELPA solve \n");
#ifdef WITH_MPI
				mpierr = MPI_Abort(MPI_COMM_WORLD, 99);
#endif
			}


			if (myid == 0) {
				printf("\n");
				printf("2stage ELPA complex solver complete\n");
				printf("\n");
			}

			/* check the results */
			status = check_correctness_complex_from_fortran(na, nev, na_rows, na_cols, as, z, ev, sc_desc, myid, tmp1, tmp2);

			if (status !=0){
				printf("The computed EVs are not correct !\n");
			}
			if (status ==0){
				if (myid == 0) {
					printf("All ok!\n");
				}
			}

			free(sc_desc);
			free(a);
			free(z);
			free(as);

			free(tmp1);
			free(tmp2);
		}
		double Residual(int myid, double err){ };
		double EigenOrth(vector<T> tmp, double* z, double & err, int myid){};
#endif
};
#endif
