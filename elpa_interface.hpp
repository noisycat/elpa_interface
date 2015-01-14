/* This file describes a somewhat ghetto ELPA interface from C++ */
#include <vector>
#include "mycomplex.h"
#include "mpi.h"
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <iomanip>

extern "C" {
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
	void solve_provided_(MPI_Fint * op_comm, double** input_matrix, int* N, int* M, double** z, int *z_rows, int *z_cols, double ** ev);
	void solve_full_(MPI_Fint *, double **, int* , double **, int *, int * , double **);
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
		int np_rows, np_cols; 
 		int N, M; // actual matrix dimension
		int na_rows, na_cols; // local subsection of matrix dimension
		int nblk; // nblk: Blocking factor in block cyclic distribution
		int na; // na: size of system
		int nev; // nev: number of eigenvalues to calculate
		int nprow;
		int npcol;
		MPI_Comm the_comm; // The comm used for this entire operation
		vector<T> a; // where we coalesce a vector< vector< T > > before passing it to its doom (overwritten during routine)
		vector<T> as; // store original a for checking against later - can remove this
		vector<T> ev; // where we store the answered eigenvectors
		vector<T> ew; // where we store the answered eigenvalues
		vector<int> sc_desc; //scalapack descriptor?

	public:
		// Constructor 
		ELPA_Interface() {
			this->sc_desc.resize(9);
		};
		// Associate - transfer matrix
		void Associate(vector< vector<T> > &A, int nblk = 16) 
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

			this->a.empty(); // shouldn't be necessary, but we're gonna be safe
			this->a.resize(this->N * this->N);
			for(int i = 0; i < N; i++) {
				for(int j = 0; j < M; j++) {
					this->a[M*i+j] = (A[i])[j];
					std::cout << ' ' << std::setw(3) << (A[i])[j];
				}
				std::cout << std::endl;
			}
		};
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
		void Solve2(vector< vector<T> > &A, MPI_Comm* the_comm, int myid, int nprocs, int required_mpi_thread_level, int provided_mpi_thread_level) {

			if(myid==0) this->Associate(A);

			double *data = this->a.data();
			double *eigvecs, *eigvals;
			int     eigvecs_rows,  eigvecs_cols;
			MPI_Fint the_comm_f = MPI_Comm_c2f(*the_comm);

			//solve_provided_(&the_comm_f, &data, &this->N, &this->N, &eigvecs, &eigvecs_rows, &eigvecs_cols,  &eigvals);
			//solve_full_(&the_comm_f, &data, &this->N);
			solve_full_(&the_comm_f, &data, &this->N, &eigvecs, &eigvecs_rows, &eigvecs_cols,  &eigvals);
		};
		void Solve(vector< vector<T> > &A, MPI_Comm* the_comm, int myid, int nprocs, int required_mpi_thread_level, int provided_mpi_thread_level) {

			if(myid==0) this->Associate(A);

			double *data = this->a.data();
			double *eigvecs, *eigvals;
			int     eigvecs_rows,  eigvecs_cols;
			MPI_Fint the_comm_f = MPI_Comm_c2f(*the_comm);

			solve_provided_(&the_comm_f, &data, &this->N, &this->N, &eigvecs, &eigvecs_rows, &eigvecs_cols,  &eigvals);
		};
		double Residual(){
			vector<T> tmp1;
			tmp1.resize(this->na);
			// tmp1 =  A * Z
			//call pdgemm('N','N',na,nev,na,1.d0,as,1,1,sc_desc, z,1,1,sc_desc,0.d0,tmp1,1,1,sc_desc);

			//double* tmp2 = new double [na_rows * na_cols];
			vector<T> tmp2;

			// tmp2 = Zi*EVi
			tmp2.assign(this->ew.begin(), this->ew.end());
			for(int i=0; i < nev; i++){
				//pdscal(this->na,this->ev(i),tmp2,1,i,sc_desc,1);
			}

			//  tmp1 = A*Zi - Zi*EVi
			//tmp1(:,:) =  tmp1(:,:) - tmp2(:,:);
			for(int i = 0; i < na; i++) tmp1.data()[i] = tmp1.data()[i] - tmp2.data()[i];

			// Get maximum norm of columns of tmp1
			errmax = 0.0;
			for(int i=0; i < nev; i++){
				err = 0.0;
				//pdnrm2(na,err,tmp1,1,i,sc_desc,1);
				errmax = max<double>(errmax, err);
			}

			// Get maximum error norm over all processors
			err = errmax;
			mpierr = MPI_Allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD);
			if(myid==0) std::cout << std::endl;
			if(myid==0) std::cout << "Error Residual     :" << errmax << std::endl;

			if (errmax > 5e-12) {
				status = 1;
			}
		}
		double EigenOrth(){
			vector<T> tmp1;
			vector<T> tmp2;
			tmp1.resize(this->na);
			// tmp1 = Z**T * Z
			tmp1.assign(tmp.size(),0.0);
			pdgemm('T','N',nev,nev,na,1.0,z,1,1,sc_desc, z,1,1,sc_desc,0.0,tmp1,1,1,sc_desc);


			// tmp2 = Zi*EVi
			tmp2.assign(this->ew.begin(), this->ew.end());
			// Initialize tmp2 to unit matrix
			tmp2.assign(tmp.size(),0.0);
			pdlaset('A',nev,nev,0.0,1.0,tmp2,1,1,sc_desc);

			// tmp1 = Z**T * Z - Unit Matrix
			//tmp1(:,:) =  tmp1(:,:) - tmp2(:,:);

			// Get maximum error (max abs value in tmp1)
			err = maxval(abs(tmp1));
			call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr);
			if(myid==0) std::cout << "Error Orthogonality:" << errmax << std::endl;

			if (errmax > 5e-12) {
				status = 1;
			}
		};
		// Destructor - ELPA_Interface object goes byebye
		~ELPA_Interface() {  };
};

