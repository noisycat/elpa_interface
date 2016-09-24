subroutine gather_matrix(na, a, lda, na_cols, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm, ctxt, a_global)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: na, lda, na_cols, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm, ctxt
   real*8, intent(in) :: a(lda, na_cols)
   real*8, intent(out) :: a_global(na,na)

   integer i, j, lr, lc, myid, mpierr
   integer steps, starti, endi
   integer, allocatable :: l_row(:), l_col(:)
#if defined(PRINT_GATHER_DEBUG)
   character*16 :: filename
#endif

   ! allocate and set index arrays

   allocate(l_row(na))
   allocate(l_col(na))

   ! Mapping of global rows/cols to local

   l_row(:) = 0
   l_col(:) = 0

   lr = 0 ! local row counter
   lc = 0 ! local column counter

   call mpi_comm_rank(op_comm,myid,mpierr)
   do i = 1, na

     if( MOD((i-1)/nblk,np_rows) == my_prow) then
       ! row i is on local processor
       lr = lr+1
       l_row(i) = lr
#if defined(PRINT_GATHER_DEBUG2)
       write(*,*) myid,' has row ',i,' on ',lr
#endif
     endif

     if( MOD((i-1)/nblk,np_cols) == my_pcol) then
       ! column i is on local processor
       lc = lc+1
       l_col(i) = lc
#if defined(PRINT_GATHER_DEBUG2)
       write(*,*) myid,' has col ',i,' on ',lc
#endif
     endif

   enddo

   a_global = 0

   do i=1,na
      if(l_col(i) > 0) then
         do j=1,i
            if(l_row(j)>0) a_global(j,i) = a(l_row(j),l_col(i))
         enddo
      endif
      if(l_row(i) > 0) then
         do j=1,i-1
            if(l_col(j)>0) a_global(i,j) = a(l_row(i),l_col(j))
         enddo
      endif
   enddo

#if defined(PRINT_GATHER_DEBUG)
   write(filename,"(A6,I3.3,A4)") "gather",myid,".txt"
   open(12,file=filename,status="unknown")
   write(12,"(34E26.13E2)") (a_global(i,:),i=1,34)
   close(12)
#endif

#if 0
   ! global reduction that stitches together the a_global
   ! Doesn't work for large enough matrices (na >> 15k on edison) because of 
   ! limitations on BLACS buffer space
   call DGSUM2D(ctxt,'A',' ',na,na,a_global,na,-1,-1,-1)
#elif 1
   steps = 8
   do i=1,steps
      starti = (i-1)*(na/steps)+1
      endi = i*(na/steps)
      if(myid==0) write(*,*) i,starti,endi
      call DGSUM2D(ctxt,'A',' ',na,endi-starti+1,a_global(:,starti:endi),na,-1,-1,-1)
   enddo
   if (.not. (endi == na)) then
      starti=endi+1
      endi=na
      if(myid==0) write(*,*) i,starti,endi
      call DGSUM2D(ctxt,'A',' ',na,endi-starti+1,a_global(:,starti:endi),na,-1,-1,-1)
   endif
#endif

#if defined(PRINT_GATHER_DEBUG)
   write(filename,"(A6,I3.3,A4)") "collect",myid,".txt"
   open(12,file=filename,status="unknown")
   write(12,"(34E26.13E2)") (a_global(i,:),i=1,PRINT_GATHER_DEBUG_N)
   close(12)
#endif

   deallocate(l_row, l_col)

end subroutine gather_matrix

!----------------------------------------------------------------------------------
subroutine section_matrix(na, a, lda, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm, a_global)

   implicit none
   include 'mpif.h'

   integer, intent(in) :: na, lda, nblk, my_prow, my_pcol, np_rows, np_cols, op_comm
   real*8, intent(out) :: a(lda, *)
   real*8, intent(in) :: a_global(na,na)

   integer i, j, lr, lc, myid, mpierr
   integer, allocatable :: l_row(:), l_col(:)

   real*8, allocatable :: col(:)

#if defined(PRINT_SECTION_DEBUG)
   character*16 :: filename
#endif

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
      col(1:i) = a_global(1:i,i)
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

#if defined(PRINT_SECTION_DEBUG)
   write(filename,"(A6,I3.3,A4)") "sectio",myid,".txt"
   open(12,file=filename,status="unknown")
   write(12,"(34E26.13E2)") (a_global(i,:),i=1,34)
   close(12)
#endif

   deallocate(l_row, l_col, col)

end subroutine section_matrix
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
   character*16 :: filename

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

