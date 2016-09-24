/*     This file is part of ELPA. */
/*  */
/*     The ELPA library was originally created by the ELPA consortium, */
/*     consisting of the following organizations: */
/*  */
/*     - Max Planck Computing and Data Facility (MPCDF), formerly known as */
/*       Rechenzentrum Garching der Max-Planck-Gesellschaft (RZG), */
/*     - Bergische Universität Wuppertal, Lehrstuhl für angewandte */
/*       Informatik, */
/*     - Technische Universität München, Lehrstuhl für Informatik mit */
/*       Schwerpunkt Wissenschaftliches Rechnen , */
/*     - Fritz-Haber-Institut, Berlin, Abt. Theorie, */
/*     - Max-Plack-Institut für Mathematik in den Naturwissenschaften, */
/*       Leipzig, Abt. Komplexe Strukutren in Biologie und Kognition, */
/*       and */
/*     - IBM Deutschland GmbH */
/*  */
/*  */
/*     More information can be found here: */
/*     http://elpa.mpcdf.mpg.de/ */
/*  */
/*     ELPA is free software: you can redistribute it and/or modify */
/*     it under the terms of the version 3 of the license of the */
/*     GNU Lesser General Public License as published by the Free */
/*     Software Foundation. */
/*  */
/*     ELPA is distributed in the hope that it will be useful, */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/*     GNU Lesser General Public License for more details. */
/*  */
/*     You should have received a copy of the GNU Lesser General Public License */
/*     along with ELPA.  If not, see <http://www.gnu.org/licenses/> */
/*  */
/*     ELPA reflects a substantial effort on the part of the original */
/*     ELPA consortium, and we ask you to respect the spirit of the */
/*     license that we chose: i.e., please contribute any changes you */
/*     may have back to the original ELPA library distribution, and keep */
/*     any derivatives of ELPA under the same license that we chose for */
/*     the original distribution, the GNU Lesser General Public License. */
/*  */
/*  */

#include <stdio.h>
#include <stdlib.h>
#ifdef WITH_MPI
#include <mpi.h>
#endif
#include <math.h>

#include <elpa/elpa.h>

void 	pdtran_ (int *M, int *N, double *ALPHA, double *A, int *IA, int *JA, int *DESCA, double *BETA, double *C, int *IC, int *JC, int *DESCC);

int main(int argc, char** argv) {
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

   double *a, *z, *as, *ev, *tmp1, *tmp2;

   int iseed;

   int success;

   int useQr, THIS_REAL_ELPA_KERNEL_API;
#ifdef WITH_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#else
   nprocs = 1;
   myid=0;
   MPI_COMM_WORLD=1;
#endif
   na = 1000;
   nev = 500;
   nblk = 16;

   if (myid == 0) {
     printf("This is the c version of an ELPA test-programm\n");
     printf("\n");
     printf("It will call the 1stage ELPA real solver for an\n");
     printf("of matrix size %d. It will compute %d eigenvalues\n",na,nev);
     printf("and uses a blocksize of %d\n",nblk);
     printf("\n");
     printf("This is an example program with much less functionality\n");
     printf("as it's Fortran counterpart. It's only purpose is to show how \n");
     printf("to evoke ELPA1 from a c programm\n");
     printf("\n");

   }
   double* input_matrix = (double*) calloc(na*na,sizeof(double));
   double* test_matrix = (double*) calloc(na*na,sizeof(double));
   // might be a problem here
   srand(0);
   size_t c_i, c_j;
   for(c_i = 0; c_i <  (size_t)na; c_i++){
	   for(c_j = c_i; c_j <  (size_t)na; c_j++){
		   double tmp = rand() / (double)RAND_MAX;
		   input_matrix[c_i*na + c_j] = tmp;
		   input_matrix[c_j*na + c_i] += tmp;
	   }
   }
   for(c_i = 0; c_i <  (size_t)(na*na); c_i++)
	   test_matrix[c_i] = input_matrix[c_i];

   char jobz = 'V';
   char uplo = 'U';
   int lda = na;
   int lwork = na*na;
   double* eigs = (double*) calloc(na,sizeof(double));
   double* work = (double*) calloc(na*na,sizeof(double));
   dsyev_(&jobz, &uplo, &na, test_matrix, &lda, eigs, work, &lwork, &info);
   if (myid==0) printf("| Accurate solve complete work(1)-%lf info-%d\n",work[0],info);
   if(0 && myid == 0){
	   FILE* w_output = fopen("eigs_output.txt","w");
	   for(c_i = 0; c_i < (size_t) nev; c_i++)
	   {
		   fprintf(stdout,"%03d %lf\n",c_i,eigs[c_i]);
		   fprintf(w_output,"%03d %lf\n",c_i,eigs[c_i]);
	   }
	   fclose(w_output);
   }

   if (myid==0) printf("| Random matrix block has been set up. (only processor 0 confirms this step)\n");

   double double_one=1.0;
   int int_one=1;

   if (myid==0) printf("| Random matrix has been symmetrized.\n");

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
   my_blacs_ctxt =  MPI_Comm_c2f(MPI_COMM_WORLD);
   const char order = 'C';
   const int zero_int = 0;
#else
  my_mpi_comm_world = 1;
#endif
   blacs_gridinit_(&my_blacs_ctxt, "C", &np_rows, &np_cols);
   blacs_gridinfo_(&my_blacs_ctxt, &nprow, &npcol, &my_prow, &my_pcol);

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
   na_rows = numroc_(&na, &nblk, &my_prow, &zero_int, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero_int, &np_cols);

   sc_desc = (int*) malloc(9*sizeof(int));

   descinit_( sc_desc, &na, &na, &nblk, &nblk, &zero_int, &zero_int, &my_blacs_ctxt, &na_rows, &info );

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

   a  = malloc(na_rows*na_cols*sizeof(double));
   z  = malloc(na_rows*na_cols*sizeof(double));
   as = malloc(na_rows*na_cols*sizeof(double));


   ev = malloc(na*sizeof(double));

   tmp1 = malloc(na_rows*na_cols*sizeof(double));
   tmp2 = malloc(na_rows*na_cols*sizeof(double));

   section_matrix_(&na, a, &na_rows, &nblk, &my_prow, &my_pcol, &np_rows, &np_cols, &my_mpi_comm_world, input_matrix);

   // backup our part of the matrix
   for(c_i=0; c_i < (size_t)(na_rows*na_cols); c_i++)
   	z[c_i] = a[c_i];

   if (myid == 0) {
     printf("\n");
     printf("Entering ELPA 2stage real solver\n");
     printf("\n");
   }
#ifdef WITH_MPI
   mpierr = MPI_Barrier(MPI_COMM_WORLD);
#endif
   useQr = 0;
   THIS_REAL_ELPA_KERNEL_API = ELPA2_REAL_KERNEL_GENERIC;

   success = elpa_solve_evp_real_2stage(na, nev, a, na_rows, ev, z, na_rows, nblk, na_cols, mpi_comm_rows, mpi_comm_cols, my_mpi_comm_world, THIS_REAL_ELPA_KERNEL_API, useQr);

   if (success != 1) {
     printf("error in ELPA solve \n");
#ifdef WITH_MPI
     mpierr = MPI_Abort(MPI_COMM_WORLD, 99);
#endif
   }
   if(myid == 0){
	   FILE* ev_output = fopen("ev_ouput.txt","w");
	   for(c_i = 0; c_i < (size_t) nev; c_i++)
	   {
		   fprintf(ev_output,"%03d %lf %lf %lg\n",c_i,ev[c_i],eigs[c_i],ev[c_i]-eigs[c_i]);
	   }
	   fclose(ev_output);
   }


   if (myid == 0) {
     printf("\n");
     printf("2stage ELPA real solver complete\n");
     printf("\n");
   }

   /* check the results */
   gather_matrix_(&na, z, &na_rows, &na_cols, &nblk, &my_prow, &my_pcol, &np_rows, &np_cols, &my_mpi_comm_world, &my_blacs_ctxt, input_matrix);
   status = 0;

   if (myid==0) printf("| input matrix has been rewritten.\n");

   if (status !=0){
     printf("The computed EVs are not correct !\n");
   }
   if (status ==0){
     if (myid ==0) {
       printf("All ok!\n");
     }
   }

   free(sc_desc);
   free(a);
   free(z);
   free(as);

   free(tmp1);
   free(tmp2);
#ifdef WITH_MPI
   blacs_gridexit_(&my_blacs_ctxt);
   MPI_Finalize();
#endif
   return 0;
}
