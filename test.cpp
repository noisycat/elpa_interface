#include <omp.h>
#include <mpi.h>

int main(int argc, char** argv){
	int required = MPI_THREAD_MULTIPLE;
	int provided;
	MPI_Init_thread(&argc,&argv,required,&provided);
	int myid;
	int ntasks;
	MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#pragma omp parallel
{
	int nthreads = omp_get_num_threads();
	if(myid == 0) printf("%d %d\n",ntasks,nthreads);
}
	MPI_Finalize();
}
