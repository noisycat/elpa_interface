#include <vector>
#include <iostream>
#include "elpa_interface.hpp"
#include "test_functions.hpp"

using std::vector;

int main()
{
	/* get N, M values from commandline for tests */
	int N = 500;
	int M = 500;

	/* MPI TASK value */
	int myid = 0;

	/* test 'matrix' */
	vector< vector<double> > A;

	/* make sure that Matrix is symmetric - ELPA requires users to ensure themselves it is */
	Test_FillA(A,N,M,myid);

	/* ELPA_Interface */
	ELPA_Interface<double> elpa;

	/* Associate Interface with Matrix */
	elpa.Associate(&A);

	/* ELPA - BLACS Comms Gen */
	/* ELPA Comms Gen */
	/* ELPA Solve */

	return 0;
}
