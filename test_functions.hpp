#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
void Test_FillA(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		A.push_back(row);
		double sum = 0.0;
		for( int j = 0; j < M; j++) {
			const double my_random = rand();
			row.push_back( my_random );
		}
		Z.push_back(row);
	}

	// make A symmetric - same way they did it in ELPA test: A = Z + Z.Transpose()
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			A[i][j] = Z[i][j] + Z[j][i];
		}
	}
}
void Test_FillA_Set1(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		A.push_back(row);
		double sum = 0.0;
		for( int j = 0; j < M; j++) {
			const double my_random = rand();
			row.push_back( my_random );
		}
		Z.push_back(row);
	}

	// make A symmetric - same way they did it in ELPA test: A = Z + Z.Transpose()
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			A[i][j] = Z[i][j] + Z[j][i];
		}
	}
}
void Test_FillA_Set2(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		A.push_back(row);
		double sum = 0.0;
		for( int j = 0; j < M; j++) {
			const double my_random = rand();
			row.push_back( my_random );
		}
		Z.push_back(row);
	}

	// make A symmetric - same way they did it in ELPA test: A = Z + Z.Transpose()
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			A[i][j] = Z[i][j] + Z[j][i];
		}
	}
}
void Test_FillA_Set3(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		A.push_back(row);
		double sum = 0.0;
		for( int j = 0; j < M; j++) {
			const double my_random = rand();
			row.push_back( my_random );
		}
		Z.push_back(row);
	}

	// make A symmetric - same way they did it in ELPA test: A = Z + Z.Transpose()
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			A[i][j] = Z[i][j] + Z[j][i];
		}
	}
}
void Test_FillA_Set4(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( 0 );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		for( int j = 0; j < M; j++) {
			double my_random = rand() /((float) RAND_MAX);
			row.push_back( my_random );
		}
		A.push_back(row);
		Z.push_back(row);
	}

	// make A symmetric - same way they did it in ELPA test: A = Z + Z.Transpose()
	for(int i = 0; i < N; i++) {
		for(int j = 0; j < M; j++) {
			A[i][j] = Z[i][j] + Z[j][i];
		}
	}
}
void Test_FillA_CommsTest(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	A.clear();

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.reserve(M);
		double sum = 0.0;
		for( int j = 0; j < M; j++) {
			const double my_random = j + M*i;
			row.push_back( my_random );
			std::cout << ' ' << std::setw(3) << my_random;
		}
		A.push_back(row);
		std::cout << std::endl;
	}
}
int Test_Greeting(int myid, int nprocs)
{
	printf("%d of %d\n",myid,nprocs);
	return 0;
}
