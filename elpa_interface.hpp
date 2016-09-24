/* This file describes a somewhat ghetto ELPA interface from C++ */
#ifndef __ELPA_INTERFACE_HPP__
#define __ELPA_INTERFACE_HPP__

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include "extern_utility.h"
#include <mpi.h>

#include <elpa/elpa.h>

extern "C" {

	// pretty sure this works?
	int numroc_(int* na, int* nblk, int* my_proc_coord, const int* Offset_coord, int* num_coord);
	int descinit_(int* sc_desc, int* N, int* M, int* nblkX, int* nblkY, const int* offset_Z, const int* offset_Y, MPI_Fint* my_blacs_ctxt, int* na_rows, int* info);
	int blacs_gridinit_(int* my_blacs_ctxt, char[1], int* np_rows, int* np_cols);
	int blacs_gridinfo_(int* my_blacs_ctxt, int* nprow, int* npcol, int* my_prow, int* my_pcol);
	int blacs_gridexit_(int* my_blacs_ctxt);

	// original fake distributed H~ way
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
		// Required for the fortran provided BLACS / Scalapack interfaces 
		const int zero_int = 0;
		const double zero_dbl = 0.0;
		const float zero_flt = 0.f;
		const int one_int = 1;
		const double one_dbl = 1.0;
		const float one_flt = 1.f;

		int N, M; // actual matrix dimension
		int nblk; // nblk: Blocking factor in block cyclic distribution
		int na; // na: size of system
		int nev; // nev: number of eigenvalues to calculate

		double* A_mm; // personal block of distributed matrix A
		double* Q_mm; // personal block of distributed eigenvectors Q
		double* A; // distributed matrix A
		double* ev; // eigenvalues

		int myid, nprocs;
		// Base MPI / BLACS Communicators
		//MPI_Comm my_mpi_comm_world, my_blacs_ctxt;


		MPI_Fint my_mpi_comm_world, my_blacs_ctxt;

		int np_rows, np_cols;

		MPI_Fint nprow, npcol, my_prow, my_pcol;

		// ELPA Communicators
		MPI_Fint mpi_comm_rows, mpi_comm_cols;

		//
		int info;
		int* sc_desc;
 
		//
		int na_rows, na_cols;

		char order;

		int useQr;

		std::ostream& outfile;
	public:
		// Constructor 
#ifdef WITH_MPI
		ELPA_Interface(const MPI_Comm _comm=MPI_COMM_WORLD) : 
#else
		ELPA_Interface(const int _comm=1) : 
#endif
			N(0), M(0), nblk(16), na(0), nev(0), 
			A_mm(NULL), Q_mm(NULL), A(NULL), ev(NULL),
			np_rows(0), np_cols(0), my_prow(0), my_pcol(0), 
			order('C'), 
			myid(0), nprocs(1), my_mpi_comm_world(1), my_blacs_ctxt(1), 
			info(0), sc_desc(NULL), 
			useQr(0), outfile(std::cout) {
#ifdef WITH_MPI
			MPI_Comm_size(_comm, &this->nprocs);
			MPI_Comm_rank(_comm, &this->myid);
			this->my_mpi_comm_world = MPI_Comm_c2f(_comm);
#endif
		};
		// Destructor - ELPA_Interface object goes byebye
		~ELPA_Interface()
		{
			if(sc_desc) free(sc_desc);
		}

		int StartBLACS() {
#ifdef WITH_MPI
			this->my_blacs_ctxt =  MPI_Comm_c2f(MPI_COMM_WORLD);
			double startVal = sqrt((double) nprocs);
			int np_colsStart = (int) round(startVal);
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

			blacs_gridinit_(&this->my_blacs_ctxt, "C", &this->np_rows, &this->np_cols);
			blacs_gridinfo_(&this->my_blacs_ctxt, &this->nprow, &this->npcol, &this->my_prow, &this->my_pcol);
#endif
			if (myid == 0) {
				printf("\n\nPast BLACS_Gridinfo...\n");
			}
		};

		void StopBLACS() {
			blacs_gridexit_(&this->my_blacs_ctxt);
		};

		int GetELPAComms()
		{

			int mpierr = get_elpa_communicators(my_mpi_comm_world, my_prow, my_pcol, &mpi_comm_rows, &mpi_comm_cols);
			if (myid == 0) {
				printf("\nPast split elpa communicator setup for rows and columns...\n");
				if (mpierr == 0) {
					printf("..But there was an error");
				}
			}
			return mpierr;
		};

		std::pair<int,int> MatrixLocalSize(int _na, 
				int _nblk_row, int _nblk_col, 
				int _my_prow, int _my_pcol, 
				int _np_rows, int _np_cols){
			na_rows = numroc_(&_na, &_nblk_row, &_my_prow, &zero_int, &_np_rows);
			na_cols = numroc_(&_na, &_nblk_col, &_my_pcol, &zero_int, &_np_cols);
			if (myid == 0) printf("\nPast numroc_ calculations for rows and columns...\n");
			return std::pair<int,int>(na_rows,na_cols);
		}
		std::pair<int,int> MatrixLocalSize() {
			return MatrixLocalSize(this->na,
				this->nblk, this->nblk, 
				this->my_prow, this->my_pcol, 
				this->np_rows, this->np_cols);
		}

		int ScaLAPACKDesc(){
			if(! sc_desc) sc_desc = (int*) malloc(9*sizeof(int));
			descinit_( sc_desc, &na, &na, &nblk, &nblk,
				       	&zero_int, &zero_int, &my_blacs_ctxt, &na_rows, &info );
			if (myid == 0) printf("\nPast scalapack descriptor init...\n");
		};
		int SolveEVP2Stage(int elpaKernel = ELPA2_REAL_KERNEL_GENERIC){
			if (myid == 0) printf("\nEntering ELPA 2stage real solver\n");
			assert(na_rows);
			assert(na_cols);
			assert(nblk);
			assert(na);
			assert(nev);
			assert(A_mm);
			assert(Q_mm);
			assert(ev);

			int success = elpa_solve_evp_real_2stage(na, nev, A_mm, na_rows, ev, Q_mm, na_rows, nblk, na_cols, mpi_comm_rows, mpi_comm_cols, my_mpi_comm_world, elpaKernel, useQr);
			if (success != 1) {
				printf("%d - error in ELPA solve \n",this->myid);
#ifdef WITH_MPI
				int mpierr = MPI_Abort(MPI_COMM_WORLD, 99);
#endif
			}
			if (myid == 0) printf("\nExiting ELPA 2stage real solver\n");
			return 0;
		};



		void Transfer(vector< vector<T> > &A, T* a, int nblk = 16){};
		void Broadcast(vector< vector<T> > &A, T* a, int nblk = 16){};
		void Gather(vector< vector<T> > &A, T* a, int nblk = 16){};

		// Take a portion of a matrix - as described by numroc_ - and use that instead of a full matrix
		void AssociatePartial(vector< vector<T> > &A, T* a, int nblk = 16) {};

		// Associate - transfer matrix
		void Associate(vector< vector<T> > &A, double* a, int nblk = 16) 
		{ 
			this->N = A.size();
			this->M = A[0].size();
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
		int MatrixLocalSetupOriginal(T* Ain){
			this->A = Ain;
			return this->MatrixLocalSetupOriginal();
		};
		int MatrixLocalSetupOriginal(){
			assert(A != NULL);
			this->A_mm = new double[na_rows * na_cols];
			this->Q_mm = new double[na_rows * na_cols];
			section_matrix_(&na, A_mm, &na_rows, &nblk, &my_prow, &my_pcol, &np_rows, &np_cols, &my_mpi_comm_world, A);
			return 0;
		};
		int MatrixLocalGather(){
			gather_matrix_(&na, Q_mm, &na_rows, &na_cols, &nblk, &my_prow, &my_pcol, &np_rows, &np_cols, &my_mpi_comm_world, &my_blacs_ctxt, A);
		}
		int EigenDecompOriginal(double* WTHW, int n_WTHW, double* eigvals, int n_ev){
			/* This calls the pieces for the eigen decomposition of H~: H~ W = W D --> W,D
			 * It uses the following parameters:
			 *
			 */
			assert(this->na == n_WTHW);
			assert(WTHW != NULL);
			assert(eigvals != NULL);
			this->A = WTHW;
			this->ev = eigvals;
			this->nev = n_ev;

			this->StartBLACS();
			this->GetELPAComms();
			this->ScaLAPACKDesc();
			this->MatrixLocalSize();
			this->MatrixLocalSetupOriginal();
			this->SolveEVP2Stage();
			this->MatrixLocalGather();
			this->StopBLACS();
		};
#if defined(UNIT_TESTS)
		void Test_Solve(vector<vector<double> >& A, double* a, int N, double* eigvals, int nev, MPI_Comm the_comm, char filename[] = "eigvals_Solve.txt")
		{
			for(size_t i = 0; i < (size_t)N; a[i++]=0.0);
			for(size_t i = 0; i < (size_t)nev; eigvals[i++]=0.0);
			this->Associate(A,a);
			this->Solve(a,N,&the_comm,eigvals);

			if(myid==0) {
				FILE* comparison = fopen(filename,"w");
				for(int i = 0; i < N; i++) fprintf(comparison,"%d %e\n",i,eigvals[i]);
				fclose(comparison);
			}
		};
		void Test_EigenDecompOriginal(vector<vector<double> >& A, double* a, int N, double* eigvals, int nev, MPI_Comm the_comm, char filename[] = "eigvals_EigenDecompOriginal.txt")
		{
			for(size_t i = 0; i < (size_t)N; a[i++]=0.0);
			for(size_t i = 0; i < (size_t)nev; eigvals[i++]=0.0);
			this->Associate(A,a);
			this->EigenDecompOriginal(a,N,eigvals,nev);

			if(myid==0) {
				FILE* comparison = fopen(filename,"w");
				for(int i = 0; i < N; i++) fprintf(comparison,"%d %e\n",i,eigvals[i]);
				fclose(comparison);
			}
		};
#endif 
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
