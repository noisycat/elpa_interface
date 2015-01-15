#include "config-f90.h"
subroutine solve_full(op_comm, input_matrix, N, &
z_ptr_c, na_rows_out, na_cols_out, ev_ptr_c)

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

   use ELPA1
   use ELPA2
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
   double precision, dimension(1:N,1:N), intent(inout) :: input_matrix 
   type(c_ptr)  :: z_ptr_c, ev_ptr_c
   integer, intent(out) :: na_rows_out, na_cols_out

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:   System size
   ! nev:  Number of eigenvectors to be calculated
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer :: nblk
   integer na, nev

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol

   integer, external :: numroc

   real*8 err, errmax
   real*8, allocatable :: a(:,:), z(:,:), tmp1(:,:), tmp2(:,:), as(:,:), ev(:)

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

#ifndef HAVE_ISO_FORTRAN_ENV
  integer, parameter   :: error_unit = 6
#endif

  logical :: success

   write_to_file = .false.
   success = .true.

   nblk = 16
#if 1
   na = N
   nev = N
#else
   na = 20
   nev = 20
#endif
#if 0
   if (COMMAND_ARGUMENT_COUNT() == 3) then
      call GET_COMMAND_ARGUMENT(1, arg1)
      call GET_COMMAND_ARGUMENT(2, arg2)
      call GET_COMMAND_ARGUMENT(3, arg3)

      read(arg1, *) na
      read(arg2, *) nev
      read(arg3, *) nblk
   endif

   if (COMMAND_ARGUMENT_COUNT() == 4) then
      call GET_COMMAND_ARGUMENT(1, arg1)
      call GET_COMMAND_ARGUMENT(2, arg2)
      call GET_COMMAND_ARGUMENT(3, arg3)
      call GET_COMMAND_ARGUMENT(4, arg4)
      read(arg1, *) na
      read(arg2, *) nev
      read(arg3, *) nblk
      
   endif
#endif
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

   call get_elpa_row_col_comms(mpi_comm_world, my_prow, my_pcol, &
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
   allocate(as(na_rows,na_cols))

   allocate(ev(na))

   ! For getting a symmetric test matrix A we get a random matrix Z
   ! and calculate A = Z + Z**T

   ! We want different random numbers on every process
   ! (otherways A might get rank deficient):

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

   ! Save original matrix A for later accuracy checks

   as = a

   ! set print flag in elpa1
   elpa_print_times = .true.

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors

   if (myid==0) then
     print '(a)','| Entering two-stage ELPA solver ... '
     print *
   end if

   call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
   success = solve_evp_real_2stage(na, nev, a, na_rows, ev, z, na_rows, nblk, &
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
      


   if(write_to_file) then
      if (myid == 0) then
         open(17,file="EVs_real2_out.txt",form='formatted',status='new')
         do i=1,na
            write(17,*) i,ev(i)
         enddo
         close(17)
      endif
   endif
   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)

   deallocate(a)
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

   deallocate(z)
   deallocate(tmp1)
   deallocate(tmp2)
   deallocate(ev)
   call blacs_gridexit(my_blacs_ctxt)
end

!-------------------------------------------------------------------------------

subroutine solve_provided(op_comm, input_matrix, N, M, z_ptr_c, &
                na_rows, na_cols, ev_ptr_c)
!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program makes use of ELPA for a C++-callable routine.
! MVP
!-------------------------------------------------------------------------------

   use ELPA1
   use test_util
   use iso_c_binding

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer, parameter :: nblk = 16

   !-------------------------------------------------------------------------------
   ! INPUT:
   !
   !
   !-------------------------------------------------------------------------------
   integer, intent(in) :: op_comm
   integer, intent(in) :: N, M
   double precision, dimension(1:N,1:M), intent(inout) :: input_matrix 
   type(C_PTR)  :: z_ptr_c, ev_ptr_c
   integer, intent(out) :: na_rows, na_cols
   !-------------------------------------------------------------------------------
   !  Local Variables

   integer na, nev

   integer np_rows, np_cols 

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol, lenarg

   integer, external :: numroc

   real*8 err, errmax
   real*8, allocatable :: a(:,:), tmp1(:,:), tmp2(:,:), as(:,:)
   real*8, allocatable, save :: z(:,:), ev(:)


   character*256 filename
   integer :: omp_get_max_threads
   logical :: success
   !-------------------------------------------------------------------------------
   !  MPI Initialization

   call mpi_comm_rank(op_comm,myid,mpierr)
   call mpi_comm_size(op_comm,nprocs,mpierr)

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

   my_blacs_ctxt = op_comm
   call BLACS_Gridinit( my_blacs_ctxt, 'C', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in get_elpa_row_col_comms.

   call get_elpa_row_col_comms(op_comm, my_prow, my_pcol, &
                               mpi_comm_rows, mpi_comm_cols)

!---------------------------------------------------------------------
! since we're not reading matrix size, we're passing it in, we should
! know what na is before a reading step
!---------------------------------------------------------------------
   na = N
   if(myid==0) print *,'Matrix size: ',na

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
   print *, 'id ', ' na_rows', ' na_cols', myid, na_rows, na_cols

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices 

   allocate(a (na_rows,na_cols))
   allocate(z (na_rows,na_cols))
   allocate(as(na_rows,na_cols))
   allocate(ev(na))

   !-------------------------------------------------------------------------------
   ! Read matrix

   call bcast_matrix(na, a, ubound(a,1), nblk, my_prow, my_pcol, np_rows, np_cols, op_comm, input_matrix)

   nev = na ! all eigenvaules

   ! Save original matrix A for later accuracy checks

   as = a

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors

   success = solve_evp_real(na, nev, a, na_rows, ev, z, na_rows, nblk, &
                       mpi_comm_rows, mpi_comm_cols)

   if(myid == 0) print *,'Time tridiag_real :',time_evp_fwd
   if(myid == 0) print *,'Time solve_tridi  :',time_evp_solve
   if(myid == 0) print *,'Time trans_ev_real:',time_evp_back

   if(myid == 0) then
      do i=1,nev
         print '(i6,g25.15)',i,ev(i)
      enddo
   endif

   z_ptr_c = c_loc(z)
   ev_ptr_c = c_loc(ev)
   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)

   deallocate(a)
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
   
   deallocate(tmp1)
   deallocate(tmp2)
   call blacs_gridexit(my_blacs_ctxt)

end subroutine solve_provided

!----------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------
subroutine bcast_matrix(na, a, lda, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm, a_global)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: na, lda, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm
   real*8, intent(out) :: a(lda, *)
   real*8, intent(in) :: a_global(na,na)

   integer i, j, lr, lc, myid, mpierr
   integer, allocatable :: l_row(:), l_col(:)

   real*8, allocatable :: col(:)

   ! allocate and set index arrays

   allocate(l_row(na))
   allocate(l_col(na))

   ! Mapping of global rows/cols to local

   l_row(:) = 0
   l_col(:) = 0

   lr = 0 ! local row counter
   lc = 0 ! local column counter

   do i = 1, na

     if( MOD((i-1)/nblk,np_rows) == my_prow) then
       ! row i is on local processor
       lr = lr+1
       l_row(i) = lr
     endif

     if( MOD((i-1)/nblk,np_cols) == my_pcol) then
       ! column i is on local processor
       lc = lc+1
       l_col(i) = lc
     endif

   enddo

   call mpi_comm_rank(op_comm,myid,mpierr)
   allocate(col(na))

   do i=1,na
      if(myid==0) col(1:i) = a_global(1:i,i)
      call mpi_bcast(col,i,MPI_REAL8,0,op_comm,mpierr)
      if(l_col(i) > 0) then
         do j=1,i
            if(l_row(j)>0) a(l_row(j),l_col(i)) = col(j)
         enddo
      endif
      if(l_row(i) > 0) then
         do j=1,i-1
            if(l_col(j)>0) a(l_row(i),l_col(j)) = col(j)
         enddo
      endif
   enddo

   deallocate(l_row, l_col, col)

end subroutine bcast_matrix

!-------------------------------------------------------------------------------
subroutine read_matrix(iunit, na, a, lda, nblk, my_prow, my_pcol, np_rows, np_cols)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: iunit, na, lda, nblk, my_prow, my_pcol, np_rows, np_cols
   real*8, intent(out) :: a(lda, *)

   integer i, j, lr, lc, myid, mpierr
   integer, allocatable :: l_row(:), l_col(:)

   real*8, allocatable :: col(:)

   ! allocate and set index arrays

   allocate(l_row(na))
   allocate(l_col(na))

   ! Mapping of global rows/cols to local

   l_row(:) = 0
   l_col(:) = 0

   lr = 0 ! local row counter
   lc = 0 ! local column counter

   do i = 1, na

     if( MOD((i-1)/nblk,np_rows) == my_prow) then
       ! row i is on local processor
       lr = lr+1
       l_row(i) = lr
     endif

     if( MOD((i-1)/nblk,np_cols) == my_pcol) then
       ! column i is on local processor
       lc = lc+1
       l_col(i) = lc
     endif

   enddo

   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   allocate(col(na))

   do i=1,na
      if(myid==0) read(iunit,*) col(1:i)
      call mpi_bcast(col,i,MPI_REAL8,0,MPI_COMM_WORLD,mpierr)
      if(l_col(i) > 0) then
         do j=1,i
            if(l_row(j)>0) a(l_row(j),l_col(i)) = col(j)
         enddo
      endif
      if(l_row(i) > 0) then
         do j=1,i-1
            if(l_col(j)>0) a(l_row(i),l_col(j)) = col(j)
         enddo
      endif
   enddo

   deallocate(l_row, l_col, col)

end subroutine read_matrix

subroutine solve_provided_full(op_comm, input_matrix, N, M)
!-------------------------------------------------------------------------------
! Standard eigenvalue problem - REAL version
!
! This program makes use of ELPA for a C++-callable routine.
! MVP
!-------------------------------------------------------------------------------

   use ELPA1
   use test_util

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer, parameter :: nblk = 16

   !-------------------------------------------------------------------------------
   ! INPUT:
   !
   !
   !-------------------------------------------------------------------------------
   integer, intent(in) :: op_comm
   integer, intent(in) :: N, M
   double precision, dimension(1:N,1:M), intent(inout) :: input_matrix 

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer na, nev

   integer np_rows, np_cols, na_rows, na_cols

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol, lenarg

   integer, external :: numroc

   real*8 err, errmax
   real*8, allocatable :: a(:,:), z(:,:), tmp1(:,:), tmp2(:,:), as(:,:), ev(:)

   character*256 filename
   integer :: omp_get_max_threads
   logical :: success
   !-------------------------------------------------------------------------------
   !  MPI Initialization

   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

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

   my_blacs_ctxt = op_comm
   call BLACS_Gridinit( my_blacs_ctxt, 'C', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in get_elpa_row_col_comms.

   call get_elpa_row_col_comms(op_comm, my_prow, my_pcol, &
                               mpi_comm_rows, mpi_comm_cols)

#if 0
   ! Read matrix size
   if(myid==0) read(10,*) na
   call mpi_bcast(na, 1, mpi_integer, 0, mpi_comm_world, mpierr)

   ! Quick check for plausibility
   if(na<=0 .or. na>10000000) then
      if(myid==0) print *,'Illegal value for matrix size: ',na
      call mpi_finalize(mpierr)
      stop
   endif
#else
!---------------------------------------------------------------------
! since we're not reading matrix size, we're passing it in, we should
! know what na is before a reading step
!---------------------------------------------------------------------
   na = N
#endif
   if(myid==0) print *,'Matrix size: ',na

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( sc_desc, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices 

   allocate(a (na_rows,na_cols))
   allocate(z (na_rows,na_cols))
   allocate(as(na_rows,na_cols))
   allocate(ev(na))

   !-------------------------------------------------------------------------------
   ! Read matrix

   call bcast_matrix(na, a, ubound(a,1), nblk, my_prow, my_pcol, np_rows, np_cols,op_comm,input_matrix)

   nev = na ! all eigenvaules

   ! Save original matrix A for later accuracy checks

   as = a

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors

   call mpi_barrier(mpi_comm_world, mpierr) ! for correct timings only
   success = solve_evp_real(na, nev, a, na_rows, ev, z, na_rows, nblk, &
                       mpi_comm_rows, mpi_comm_cols)

   if(myid == 0) print *,'Time tridiag_real :',time_evp_fwd
   if(myid == 0) print *,'Time solve_tridi  :',time_evp_solve
   if(myid == 0) print *,'Time trans_ev_real:',time_evp_back

   if(myid == 0) then
      do i=1,nev
         print '(i6,g25.15)',i,ev(i)
      enddo
   endif

   !-------------------------------------------------------------------------------
   ! Test correctness of result (using plain scalapack routines)

   deallocate(a)
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
   
   deallocate(z)
   deallocate(tmp1)
   deallocate(tmp2)
   deallocate(ev)
   call blacs_gridexit(my_blacs_ctxt)
   call mpi_finalize(mpierr)

end subroutine solve_provided_full
