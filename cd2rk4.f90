! This program calculates the solution of 1D wave eqn 
! using CD2-Rk4 method using MPI

      MODULE VARIABLES
      IMPLICIT NONE
      INCLUDE "mpif.h"

      INTEGER :: rank,root,ierr
      INTEGER, PARAMETER :: max=1000
      INTEGER :: status(MPI_STATUS_SIZE)
      INTEGER :: j,jg,t1,t1max,jm1,jp1,n,npro,local_n,jmin,jmax
      DOUBLE PRECISION :: a,b,h,x(max),alpha,x0,k,kh,Nc,c,dt,t,local_a,local_b
      DOUBLE PRECISION :: u0(max),un(max),unp1(max),usol(max),local_x(max)
      DOUBLE PRECISION :: u1(max),u2(max),u3(max),data_tra(max)
      DOUBLE PRECISION :: temp0(max),temp1(max),temp2(max),temp3(max)

      END MODULE VARIABLES

      PROGRAM WAVE_PROPAGATION
      USE VARIABLES
      IMPLICIT NONE

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npro,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

      Nc = 0.10d0 ; kh = 0.50d0
      c = 1.0d0 
      t1max = 400 
      x0 = 20.0d0
       
!     calculate initial solution 
      a=0.0d0 ; b=40.0d0; n=800 
      h = (b-a)/(dble(n))
      k=kh/h
      
!     define global coordinates
      do jg = 1,n+1
         x(jg) = (jg-1)*h
         alpha = 1.0d0*(x(jg)-x0)**2.0d0   
         u0(jg) = exp(-alpha)*cos(k*x(jg))
      end do        
      
!      write the initial condition
       open(5,file='ini.dat')
       write(5,*) 'Variables = X,U'
       do jg = 1,n+1
          write(5,*) x(jg),u0(jg)
       end do
       close(5)       
       
!     define local coordinates
      local_n = n/npro
      local_a = a+rank*local_n*h
      local_b = local_a+local_n*h
        
      do j=1,local_n+1
         local_x(j)=local_a+(j-1)*h
         alpha = 1.0d0*(local_x(j)-x0)**2.0d0
         un(j) = exp(-alpha)*cos(k*local_x(j))
      end do
      
      jmin=2 ; jmax=local_n+1

!     time loop starts
      dt = Nc*h/c
   
      DO t1=1,t1max
         t = t1*dt

!     RK-4 1st stage
         DO j=1,jmax+1
            data_tra(j)=un(j)
         ENDDO
         
         CALL DATA_TRANSFER

         DO j=1,jmax+1
            un(j)=data_tra(j)
         ENDDO

         DO j=jmin,jmax
            temp0(j) = un(j-1)-un(j+1)      
            u1(j) = un(j) + 0.5d0*Nc*temp0(j)
         ENDDO 
 
!     RK-4 2nd stage

         DO j=1,jmax+1
            data_tra(j)=u1(j)
         ENDDO

         CALL DATA_TRANSFER

         DO j=1,jmax+1
            u1(j)=data_tra(j)
         ENDDO

         DO j=jmin,jmax
            temp1(j) = u1(j-1)-u1(j+1)
            u2(j) = un(j) + 0.5d0*Nc*temp1(j)
         ENDDO
 
!      RK-4 3rd stage
 
         DO j=1,jmax+1
            data_tra(j)=u2(j)
         ENDDO

         CALL DATA_TRANSFER

         DO j=1,jmax+1
            u2(j)=data_tra(j)
         ENDDO
 
         DO j=jmin,jmax
            temp2(j) = u2(j-1)-u2(j+1)
            u3(j) = un(j) + Nc*temp2(j)
         ENDDO
 
!      RK-4 4th stage

         DO j=1,jmax+1
            data_tra(j)=u3(j)
         ENDDO

         CALL DATA_TRANSFER

         DO j=1,jmax+1
            u3(j)=data_tra(j)
         ENDDO

         DO j=jmin,jmax
            temp3(j)=u3(j-1)-u3(j+1)
            unp1(j) = un(j) + (1.0d0/6.0d0)*Nc* &
                 (temp0(j)+2.0d0*(temp1(j)+temp2(j))+temp3(j))
         ENDDO
         
         DO j=1,local_n+1
            un(j)=unp1(j)
         ENDDO

         IF(rank.EQ.0) WRITE(*,*) 'timecount=',t1,'time=',t

         IF(t.EQ.(t1max*dt)) THEN
            root = 0
            CALL MPI_GATHER(unp1(jmin:jmax),local_n,MPI_double_precision, &
                 usol,local_n,MPI_double_precision,root,MPI_Comm_World,ierr)
            IF(rank.EQ.0) then
              OPEN(10,FILE='soln.dat')
              WRITE(10,*) 'Variables = X,U'
              DO j=2,n+1
                 WRITE(10,*) x(j), usol(j)+t
              ENDDO
              CLOSE(10)
            ENDIF
         ENDIF
 
      ENDDO !  time loop  !

      CALL MPI_FINALIZE(ierr)
      
      END PROGRAM WAVE_PROPAGATION

!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------

      SUBROUTINE DATA_TRANSFER
      use VARIABLES 
      IMPLICIT NONE
      integer :: source,dest,tag,count,sendtag,recvtag
 
       IF(rank.eq.0) THEN
         sendtag=0 ; recvtag=1
         source=rank+1  ; dest=rank+1
         Call MPI_Send(data_tra(jmax),1,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(data_tra(jmax+1),1,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       ELSE IF (rank.EQ.npro-1) THEN
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(data_tra(jmin-1),1,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(data_tra(jmin),1,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
       ELSE
         sendtag=1   ;  recvtag=0
         source=rank-1     ;  dest=rank-1
         Call MPI_Recv(data_tra(jmin-1),1,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
         Call MPI_Send(data_tra(jmin),1,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)

         sendtag=0   ;  recvtag=1
         source=rank+1     ;  dest=rank+1
         Call MPI_Send(data_tra(jmax),1,MPI_double_precision,dest, &
               sendtag,MPI_Comm_World,ierr)
         Call MPI_Recv(data_tra(jmax+1),1,MPI_double_precision,source, &
               recvtag,MPI_Comm_World,status,ierr)
       END IF
 
      ENDSUBROUTINE DATA_TRANSFER           

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
