#include <stdlib.h>
void Test_FillA(vector< vector< double > > &A, const int N, const int M, const int my_seed)
{
	// A is the N x M matrix to fill out
	
	vector< vector<double> > Z;

	// Init random num generator
	srand( my_seed );

	// give Z vectors for rows (I assume this is how Grady did his matrix)
	
	for(int i = 0; i < N; i++) {
		vector<double> row;
		row.resize(M);
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
