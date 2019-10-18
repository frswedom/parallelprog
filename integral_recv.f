       program integral_recv
       implicit none
       include 'mpif.h'

       integer size, rank, ierr, i, n, my_n, status(MPI_STATUS_SIZE)
       double precision sum, gsum, a, b, time_start, time_end
       double precision al, bl, x, f, isum, step
       double precision time_1, time_2, time_total
       character (len=32) :: arg

       Parameter (a=0.0, b=3.14)

       f(x)= sin(x)

       call MPI_INIT(ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
       call MPi_COMM_RanK(MPI_COMM_WORLD, rank, ierr)
       
       if (rank.eq.0) time_1 = MPI_WTime()

       al = a + (b - a) * rank / size
       bl = al + (b - a) / size

       if (iargc().eq.0) then
            open(12, FILE='config')
            READ(12, *) n
            close(12)
       else
            call getarg(1, arg)
            read (arg, '(i12)') n
       endif

       time_start = MPI_WTime()

       my_n = n / size
       x = al
       step = (bl - al) / my_n
       sum = f(x) / 2
       do i = 1, my_n
            sum = sum + f(x)
            x = x + step
       end do

       sum = (sum + f(x) / 2) * step

       time_end = MPI_WTime()
       
       if(size.eq.1) then 
           open(13, FILE='output.recv.f.txt')
       else
           open(13, FILE='output.recv.f.txt', position='append')
       end if
       
       if (rank.ne.0) then
           call MPI_Send(sum, 1, MPI_DOUBLE_PRECISION, 0,
     $ 0, MPI_COMM_WORLD, ierr)
       else
           gsum = sum

           do i=1,size-1
               call MPI_RECV(isum, 1, MPI_DOUBLE_PRECISION,
     $ MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, status, ierr)
               gsum = gsum + isum
           end do
           
           time_2 = MPI_WTime()
           time_total = (time_2 - time_1) * 1000

           write (* , *) 'Integral = ', gsum
           write (* , *) 'Total time = ', time_total, 'ms'
           write(13, *) size, ' ', time_total
           close(13)
           
       end if

       call MPi_finaliZE(ierr)
       stop
       end
