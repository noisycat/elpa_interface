!#include "config-f90.h"
subroutine solve_full(op_comm, input_matrix, N, ev_ptr_c, nblk)

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program makes use of ELPA for a C++-callable routine.
! MVP
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program demonstrates the use of the ELPA module
! together with standard scalapack routines
! 
! Copyright of the original code rests with the authors inside the ELPA
! consortium. The copyright of any additional modifications shall rest
! with their original authors, but shall adhere to the licensing terms
! distributed along with the original code in the file "COPYING".
!
!-------------------------------------------------------------------------------

   use elpa1
   use elpa2
#ifdef WITH_OPENMP
   use test_util
#endif

#ifdef HAVE_ISO_FORTRAN_ENV
  use iso_fortran_env, only : error_unit
#endif
  use iso_c_binding

   implicit none
   include 'mpif.h'
   ! -- beginning of input header
   integer, intent(in) :: op_comm
   integer(c_int), intent(in) :: N
   real(c_double), dimension(1:N,1:N), intent(inout) :: input_matrix 
   real(c_double), dimension(1:N), intent(inout) :: ev_ptr_c

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:   System size
   ! nev:  Number of eigenvectors to be calculated
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer, intent(in) :: nblk
   integer na, nev

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol

   integer, external :: numroc

   real*8 err, errmax
   real*8, allocatable :: a(:,:)

#if defined(ELPA_INTERNAL_DEBUG)
   real*8, allocatable :: tmp1(:,:), tmp2(:,:), as(:,:)
#endif

   real*8, allocatable, target :: z(:,:), ev(:)

   integer :: iseed(4096) ! Random seed, size should be sufficient for every generator
   integer :: STATUS
#ifdef WITH_OPENMP
   integer :: omp_get_max_threads,  required_mpi_thread_level, provided_mpi_thread_level
#endif
   logical :: write_to_file
   !-------------------------------------------------------------------------------
   !  Parse command line argumnents, if given
   character*16 arg1
   character*16 arg2
   character*16 arg3
   character*16 arg4
   character*13 filename

#ifndef HAVE_ISO_FORTRAN_ENV
  integer, parameter   :: error_unit = 6
#endif

  logical :: success

   write_to_file = .false.
   success = .true.

   arg4 = "output"
   na = N
   nev = N

   !-------------------------------------------------------------------------------
   !  MPI Initialization

   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

   STATUS = 0
#ifdef WITH_OPENMP
   if (myid .eq. 0) then
      print *,"Threaded version of test program"
      print *,"Using ",omp_get_max_threads()," threads"
      print *," "
   endif
#endif

   if (myid .eq. 0) then
      print *," "
      print *,"This ELPA2 is build with"
#ifdef WITH_REAL_AVX_BLOCK2_KERNEL
      print *,"AVX optimized kernel (2 blocking) for real matrices"
#endif
#ifdef WITH_REAL_AVX_BLOCK4_KERNEL
      print *,"AVX optimized kernel (4 blocking) for real matrices"
#endif
#ifdef WITH_REAL_AVX_BLOCK6_KERNEL
      print *,"AVX optimized kernel (6 blocking) for real matrices"
#endif

#ifdef WITH_REAL_GENERIC_KERNEL
     print *,"GENERIC kernel for real matrices"
#endif
#ifdef WITH_REAL_GENERIC_SIMPLE_KERNEL
     print *,"GENERIC SIMPLE kernel for real matrices"
#endif
#ifdef WITH_REAL_SSE_KERNEL
     print *,"SSE ASSEMBLER kernel for real matrices"
#endif
#ifdef WITH_REAL_BGP_KERNEL
     print *,"BGP kernel for real matrices"
#endif
#ifdef WITH_REAL_BGQ_KERNEL
     print *,"BGQ kernel for real matrices"
#endif
   endif
   if (arg4 .eq. "output") then 
      write_to_file = .true.
      if (myid .eq. 0) print *,"Writing output files"
   endif
   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs
   ! and decrement it until a divisor is found.

   do np_cols = NINT(SQRT(REAL(nprocs))),2,-1
      if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_rows = nprocs/np_cols

   if(myid==0) then
      print *
      print '(a)','Standard eigenvalue problem - REAL version'
      print *
      print '(3(a,i0))','Matrix size=',na,', Number of eigenvectors=',nev,', Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !
   ! The BLACS context is only necessary for using Scalapack.
   !
   ! For ELPA, the MPI communicators along rows/cols are sufficient,
   ! and the grid setup may be done in an arbitrary way as long as it is
   ! consistent (i.e. 0<=my_prow<np_rows, 0<=my_pcol<np_cols and every
   ! process has a unique (my_prow,my_pcol) pair).

   my_blacs_ctxt = mpi_comm_world
   call BLACS_Gridinit( my_blacs_ctxt, 'C', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   if (myid==0) then
     print '(a)','| Past BLACS_Gridinfo.'
   end if

   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in get_elpa_row_col_comms.

   mpierr = get_elpa_communicators(mpi_comm_world, my_prow, my_pcol, &
                                   mpi_comm_rows, mpi_comm_cols)


   if (myid==0) then
     print '(a)','| Past split communicator setup for rows and columns.'
   end if

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )

   if (myid==0) then
     print '(a)','| Past scalapack descriptor setup.'
   end if

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem

   allocate(a (na_rows,na_cols))
   allocate(z (na_rows,na_cols))

#if defined(ELPA_INTERNAL_DEBUG)
   allocate(as(na_rows,na_cols))
#endif

   allocate(ev(na))

   ! For getting a symmetric test matrix A we get a random matrix Z
   ! and calculate A = Z + Z**T

   ! We want different random numbers on every process
   ! (otherways A might get rank deficient):
#if 0
   iseed(:) = myid
   call RANDOM_SEED(put=iseed)

   call RANDOM_NUMBER(z)

   a(:,:) = z(:,:)

   if (myid==0) then
     print '(a)','| Random matrix block has been set up. (only processor 0 confirms this step)'
   end if

   call pdtran(na, na,  1.d0, z, 1, 1, sc_desc, 1.d0, a, 1, 1, sc_desc) ! A = A + Z**T

   if (myid==0) then
     print '(a)','| Random matrix has been symmetrized.'
   end if

#else

   call section_matrix(na, a, na_rows, nblk, my_prow, my_pcol, np_rows, &
   np_cols, op_comm, input_matrix)
   if (myid==0) then
     print '(a)','| input matrix has been assigned.'
   end if

#ifdef DEBUG
   write(filename,"(A6,I3.3,A4)") "output",myid,".txt"
   open(12,file=filename,status="unknown")
   write(12,"(20E12.5)") (input_matrix(i,:),i=1,20)
   close(12)
#endif

#endif
   ! Save original matrix A for later accuracy checks

#if defined(ELPA_INTERNAL_DEBUG)
   as = a
#endif

   ! set print flag in elpa1
   elpa_print_times = .true.

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors

   if (myid==0) then
     print '(a)','| Entering two-stage ELPA solver ... '
     print *
   end if

#if defined(LOCKING_TIMING)
   call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
#endif
   success = solve_evp_real_2stage(na, nev, a, na_rows, ev, z, na_rows,  nblk, na_cols, &
                              mpi_comm_rows, mpi_comm_cols, mpi_comm_world)

   if (.not.(success)) then
      write(error_unit,*) "solve_evp_real_2stage produced an error! Aborting..."
      call MPI_ABORT(mpi_comm_world, mpierr)
   endif

   if (myid==0) then
     print '(a)','| Two-step ELPA solver complete.'
     print *
   end if

   if(myid == 0) print *,'Time transform to tridi :',time_evp_fwd
   if(myid == 0) print *,'Time solve tridi        :',time_evp_solve
   if(myid == 0) print *,'Time transform back EVs :',time_evp_back
   if(myid == 0) print *,'Total time (sum above)  :',time_evp_back+time_evp_solve+time_evp_fwd
      
#if defined(ELPA_DEBUG)
   if(write_to_file) then
      if (myid == 0) then
         open(17,file="EVs_real2_out.txt",form='formatted',status='unknown')
         do i=1,na
            write(17,*) i,ev(i)
         enddo
         close(17)
      endif
   endif
#endif
   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)
   call gather_matrix(na, z, na_rows, na_cols, nblk, my_prow, my_pcol, np_rows, &
   np_cols, op_comm, my_blacs_ctxt, input_matrix)
   if (myid==0) then
     print '(a)','| input matrix has been rewritten.'
   end if
#if defined(ELPA_INTERNAL_DEBUG)

   call BLACS_Gridinit( my_blacs_ctxt, 'C', 1, 1)
   allocate(tmp1(na_rows,na_cols))

   ! 1. Residual (maximum of || A*Zi - Zi*EVi ||)

   ! tmp1 =  A * Z
   call pdgemm('N','N',na,nev,na,1.d0,as,1,1,sc_desc, &
           z,1,1,sc_desc,0.d0,tmp1,1,1,sc_desc)

   deallocate(as)
   allocate(tmp2(na_rows,na_cols))

   ! tmp2 = Zi*EVi
   tmp2(:,:) = z(:,:)
   do i=1,nev
      call pdscal(na,ev(i),tmp2,1,i,sc_desc,1)
   enddo

   !  tmp1 = A*Zi - Zi*EVi
   tmp1(:,:) =  tmp1(:,:) - tmp2(:,:)

   ! Get maximum norm of columns of tmp1
   errmax = 0
   do i=1,nev
      err = 0
      call pdnrm2(na,err,tmp1,1,i,sc_desc,1)
      errmax = max(errmax, err)
   enddo

   ! Get maximum error norm over all processors
   err = errmax
   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
   if(myid==0) print *
   if(myid==0) print *,'Error Residual     :',errmax

   if (errmax .gt. 5e-12) then
      status = 1
   endif

   ! 2. Eigenvector orthogonality

   ! tmp1 = Z**T * Z
   tmp1 = 0
   call pdgemm('T','N',nev,nev,na,1.d0,z,1,1,sc_desc, &
           z,1,1,sc_desc,0.d0,tmp1,1,1,sc_desc)
   ! Initialize tmp2 to unit matrix
   tmp2 = 0
   call pdlaset('A',nev,nev,0.d0,1.d0,tmp2,1,1,sc_desc)

   ! tmp1 = Z**T * Z - Unit Matrix
   tmp1(:,:) =  tmp1(:,:) - tmp2(:,:)

   ! Get maximum error (max abs value in tmp1)
   err = maxval(abs(tmp1))
   call mpi_allreduce(err,errmax,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,mpierr)
   if(myid==0) print *,'Error Orthogonality:',errmax
   
   if (errmax .gt. 5e-12) then
      status = 1
   endif

   deallocate(tmp1)
   deallocate(tmp2)
#endif
   ev_ptr_c(:) = ev(:)
   deallocate(a)
   deallocate(z)
   deallocate(ev)
   call blacs_gridexit(my_blacs_ctxt)
end subroutine

!-----------------------------------------------------------------------
