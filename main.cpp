#include <vector>
#include <iostream>
#include "elpa_interface.hpp"
#include "test_functions.hpp"
#include "mpi.h"
#include <omp.h>

using std::vector;

int main(int argc, char* argv[])
{
	int required = MPI_THREAD_MULTIPLE;
	int provided = 0;
	int mpierr = 0;
	mpierr = MPI_Init_thread(&argc,&argv,required,&provided);
	/* get N, M values from commandline for tests */
	int N = 34;
	int M = N;

	/* MPI TASK value */
	int myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols;
	int i, my_blacs_ctxt, sc_desc[9], info, nprow, npcol;

	MPI_Comm the_comm = MPI_COMM_WORLD;
	mpierr = MPI_Comm_rank(the_comm, &myid);
	mpierr = MPI_Comm_size(the_comm, &nprocs);

	Test_Greeting(myid, nprocs);
	MPI_Barrier(the_comm);
	/* test 'matrix' */
	vector< vector<double> > A;

	if (myid == 0) {
		std::cout << "Threaded version of test program" << std::endl;
		std::cout << "Using " << omp_get_max_threads() << " threads" << std::endl << std::endl;
	}


	/* make sure that Matrix is symmetric - ELPA requires users to ensure themselves it is */
#if defined(TEST1)
	Test_FillA_Set1(A,N,M,myid);
#elif defined(TEST2)
	Test_FillA_Set2(A,N,M,myid);
#elif defined(TEST3)
	Test_FillA_Set3(A,N,M,myid);
#elif defined(TEST4)
	Test_FillA_Set4(A,N,M,myid);
#else
	Test_FillA_CommsTest(A,N,M,myid);
#endif

	/* ELPA_Interface */
	ELPA_Interface<double> elpa;

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!
	/* !!!! elpa.Solve(A) !!!!! */
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!
	double* eigvals = new double[N];
	elpa.Associate(A);

	elpa.Solve(elpa.data(),N,&the_comm,eigvals);

	if(myid==0) {
		FILE* comparison = fopen("comparison.txt","w");
		for(int i = 0; i < N; i++) fprintf(comparison,"%d %e\n",i,eigvals[i]);
		fclose(comparison);
	}

	delete [] eigvals;
	MPI_Finalize();
	return 0;
}
