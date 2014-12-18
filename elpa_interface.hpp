/* This file describes a somewhat ghetto ELPA interface from C++ */
#include <vector>
#include "mycomplex.h"

extern "C" {
	/* A reminder - all of these are the recorded FORTRAN interfaces that emulate Scalapack interfaces
	   so name mangling exists and is a very real thing */
	void get_elpa_row_col_comms_(int MPI_COMM_GLOBAL, int MY_PROW, int MY_PCOL, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	/* ELPA1 solvers */
	void solve_evp_real_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	/* ELPA2 solvers */
	void solve_evp_real_2stage_(int NA, int NEV, double* A, int LDA, double* EV, double* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
	void solve_evp_complex_2stage_(int NA, int NEV, mycomplex* A, int LDA, double* EV, mycomplex* Q, int LDQ, int NBLK, int MPI_COMM_ROWS, int MPI_COMM_COLS);
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
		int N, M;
		vector<int> comm_rows, comm_cols; // I assume these are going to be used coming from our friendly functions uptop
		vector<T> A; // where we coalesce a vector< vector< T > > before passing it to its doom (overwritten during routine)
		vector<T> eigvects; // where we store the answered eigenvectors
		vector<T> eigvals; // where we store the answered eigenvalues
	public:
		// Constructor 
		ELPA_Interface() { };
		// Associate - transfer matrix
		void Associate(vector< vector<T> >) { };
		// Estimate - dumps to stdout an estimated time of execution for a matrix, optionally kills if above limit
		double Estimate() { };
		double Estimate(double timelimit) { };
		// Create BLACS layout - required to solve
		void CreateBLACS() { };
		// Create ELPA specific layout - required to solve
		void CreateELPAComms() { } ;
		// Solve - required to solve
		void Solve() { };
		// Regather
		// Destructor - ELPA_Interface object goes byebye
		~ELPA_Interface() {  };
};
