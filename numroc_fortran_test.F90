PROGRAM NUMROCTEST
implicit none
include 'mpif.h'
integer :: N, NB
integer :: MYROW, MYCOL
integer :: NPROW, NPCOL
integer :: na_row, na_col
integer :: mpierr
integer :: op_comm, myid, nprocs
integer :: step
integer, external :: numroc

call MPI_INIT(mpierr)
op_comm = MPI_COMM_WORLD
call mpi_comm_rank(op_comm,myid,mpierr)
call mpi_comm_size(op_comm,nprocs,mpierr)

print *, 'id of nprocs ', myid, nprocs

NPROW = 2
NPCOL = 2

call BLACS_Gridinit(MPI_COMM_WORLD, "C", NPROW, NPCOL)
call BLACS_Gridinfo(MPI_COMM_WORLD, NPROW, NPCOL, MYROW, MYCOL)

print *, 'id ', ' myrow', ' mycol', myid, myrow, mycol
NB = 2
N = 16
do step = 1,5
N = N * 2
na_row = numroc(N,NB,MYROW,0,NPROW)
na_col = numroc(N,NB,MYCOL,0,NPCOL)

print *, 'id ', ' na_rows', ' na_cols', myid, na_row, na_col
enddo

call MPI_FINALIZE(mpierr)


END PROGRAM NUMROCTEST
