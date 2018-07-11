PROGRAM SWE

  implicit none

  REAL,PARAMETER                             ::  Ho = 1000.                  ! equilibrium height (m)
  REAL,PARAMETER                             ::  A = 15.                     ! amplitude
  REAL,PARAMETER                             ::  Dx = 100000.                ! space step (m)
  REAL,PARAMETER                             ::  Dt = 600.                   ! time step
  REAL,PARAMETER                             ::  sigma = 5*Dx                ! width of the gaussian
  REAL,PARAMETER                             ::  L = 62*Dx                   ! total lenght of domain
  REAL,PARAMETER                             ::  Xm = L/2                    ! middle point
  REAL,PARAMETER                             ::  f = 1.26e-4                 ! Coriolis parameter
  REAL,PARAMETER                             ::  g = 10.                     ! Gravity acceleration

  INTEGER,PARAMETER                          ::  N1 = 1                      ! 1st time level for printout
  INTEGER,PARAMETER                          ::  N2 = 10                     ! 2nd time level for printout
  INTEGER,PARAMETER                          ::  N3 = 19                     ! 3rd time level for printout
  INTEGER,PARAMETER                          ::  N4 = 28                     ! 4th time level for printout
  INTEGER,PARAMETER                          ::  N5 = 37                     ! 5th time level for printout
  INTEGER,PARAMETER                          ::  N6 = 61                     ! 6th time level for printout
  INTEGER,PARAMETER                          ::  NMAX = 61                   ! total number of time steps
  INTEGER,PARAMETER                          ::  JM = 63                     ! total number of space steps
  INTEGER,PARAMETER                          ::  NM = 3                      ! number of time steps saved at any time (number of columns per matrix)
  INTEGER,PARAMETER                          ::  NP = 6                      ! number of time steps saved to plot 

  INTEGER                                    ::  j,N                         ! counters in space and time
  INTEGER                                    ::  isave                       ! auxiliary integer
  INTEGER                                    ::  iold                        ! previous time step (n-1)
  INTEGER                                    ::  icur                        ! current time step (n)
  INTEGER                                    ::  inew                        ! future time step (n+1)
  INTEGER                                    ::  res                         ! 
  REAL                                       ::  h_bound                     ! boundary value (x=0,L)

  REAL,DIMENSION(JM,NM)                      ::  u                           ! local velocity (x)
  REAL,DIMENSION(JM,NM)                      ::  v                           ! meridional velocity (y)
  REAL,DIMENSION(JM,NM)                      ::  h                           ! 
  REAL,DIMENSION(JM,NP)                      ::  h_plot                      !
  REAL,DIMENSION(JM)                         ::  X                           ! height
  REAL,DIMENSION(JM)                         ::  alpha                       ! relaxation parameter
  REAL,DIMENSION(JM)                         ::  vg_plot                     ! meridional geostrophic velocity for print
  REAL,DIMENSION(JM)                         ::  v_plot                      ! meridional velocity for print


  OPEN(UNIT=10,FILE='SWE.txt',FORM="FORMATTED",IOSTAT=res)
  OPEN(UNIT=11,FILE='vg.txt',FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
     PRINT *, "Error in opening output file, status: ", res
     STOP					! Stop the program
  END IF


  DO j=1,JM                                ! fill up the height matrix
     X(j) = (j-1)*Dx
  END DO


  h_bound = Ho + A*exp(-(X(j)-Xm/sigma)**2)             ! boundary value for point h. on notes


  DO j=7,JM-6                                 ! fill up the alpha array
     alpha(j) = 0                             ! points from 7 to JM-6 (inner domain)
  END DO                                      ! (j(7) and j(JM-6) not included into inner domain)


  alpha(1)=1.                 ! points from 1 to 6
  alpha(2)=0.69               ! left buffer domain
  alpha(3)=0.44
  alpha(4)=0.25
  alpha(5)=0.11
  alpha(6)=0.03

  alpha(JM-5)=0.03            ! points from JM-5 to JM 
  alpha(JM-4)=0.11            ! right buffer domain
  alpha(JM-3)=0.25
  alpha(JM-2)=0.44
  alpha(JM-1)=0.69
  alpha(JM)=1.


  isave = 0                   ! Initialize ...
  iold = 1
  icur = 2
  inew = 3

!======================= n=1 ===========================

  DO j=1,JM                                        ! IC
     h(j,icur)=Ho + A*exp(-((X(j)-Xm)/sigma)**2)
     u(j,icur)=0
     v(j,icur)=0
  END DO

  DO j=2,JM                                              ! change IC for v as per point j. in notes
     v(j,icur)=(g/(f*2*Dx))*(h(j+1,icur)-h(j-1,icur))
  END DO

  DO j=1,JM                                       ! FRS on IC
     u(j,icur) = (1-alpha(j))*u(j,icur)
     v(j,icur) = (1-alpha(j))*v(j,icur)
     h(j,icur) = (1-alpha(j))*h(j,icur)+alpha(j)*h_bound
  END DO

  DO j=1,JM                                        ! Save values of IC for print
     h_plot(j,1)=h(j,icur)
  END DO

!======================= n=2 ======================== 
!! FTCS step

  DO j=2,JM-1
     h(j,inew)=h(j,icur)-u(j,icur)*(Dt/(2*Dx))*(h(j+1,icur)-h(j-1,icur))-h(j,icur)*(Dt/(2*Dx))*(u(j+1,icur)-u(j-1,icur))
     u(j,inew)=u(j,icur)+Dt*f*v(j,icur)-(Dt/(2*Dx))*u(j,icur)*(u(j+1,icur)-u(j-1,icur))-g*(Dt/(2*Dx))*(h(j+1,icur)-h(j-1,icur))
     v(j,inew)=v(j,icur)-Dt*f*u(j,icur)-(Dt/(2*Dx))*u(j,icur)*(v(j+1,icur)-v(j-1,icur))
  END DO


  u(1,inew) = 0                          !! BCs
  u(JM,inew) = 0
  v(1,inew) = 0
  v(JM,inew) = 0
  h(1,inew) = h_bound
  h(JM,inew) = h_bound


  DO j=1,JM                              ! swop the values
     h(j,iold)=h(j,icur)
     h(j,icur)=h(j,inew)
     u(j,iold)=u(j,icur)
     u(j,icur)=u(j,inew)
     v(j,iold)=v(j,icur)
     v(j,icur)=v(j,inew)
  END DO


  DO j=1,JM                                       ! FRS on FTCS n=2
     u(j,icur) = (1-alpha(j))*u(j,icur)
     v(j,icur) = (1-alpha(j))*v(j,icur)
     h(j,icur) = (1-alpha(j))*h(j,icur)+alpha(j)*h_bound
  END DO


!!====================== n=3 upwards =======================
!! Leapfrog (CTCS)

  DO N=3,NMAX
     
     DO j=2,JM-1
        h(j,inew)=h(j,iold)-(Dt/Dx)*u(j,icur)*(h(j+1,icur)-h(j-1,icur))-(Dt/Dx)*h(j,icur)*(u(j+1,icur)-u(j-1,icur))
        u(j,inew)=u(j,iold)+2*Dt*f*v(j,icur)-(Dt/Dx)*u(j,icur)*(u(j+1,icur)-u(j-1,icur))-g*(Dt/Dx)*(h(j+1,icur)-h(j-1,icur))
        v(j,inew)=v(j,iold)-2*Dt*f*u(j,icur)-(Dt/Dx)*u(j,icur)*(v(j+1,icur)-v(j-1,icur))
     END DO

     u(1,inew) = 0                          !! BCs
     u(JM,inew) = 0
     v(1,inew) = 0
     v(JM,inew) = 0
     h(1,inew) = h_bound
     h(JM,inew) = h_bound

     DO j=1,JM                                       ! FRS on CTCS n=3 upwards
        u(j,inew) = (1-alpha(j))*u(j,inew)
        v(j,inew) = (1-alpha(j))*v(j,inew)
        h(j,inew) = (1-alpha(j))*h(j,inew)+alpha(j)*h_bound
     END DO

 
     IF(N.EQ.N2) THEN                !! Ploting of relevant time steps
        DO j=1,JM
           h_plot(j,2)=h(j,inew)
        END DO
     ELSEIF(N.EQ.N3) THEN
        DO j=1,JM
           h_plot(j,3)=h(j,inew)
        END DO
     ELSEIF(N.EQ.N4) THEN
        DO j=1,JM
           h_plot(j,4)=h(j,inew)
        END DO
     ELSEIF(N.EQ.N5) THEN
        DO j=1,JM
           h_plot(j,5)=h(j,inew)
        END DO
     ELSEIF(N.EQ.N6) THEN
        DO j=1,JM
           h_plot(j,6)=h(j,inew)
        END DO
     ENDIF


     IF(N.EQ.N5) THEN                                             ! as per point i. in notes
        DO j=2,JM-1
           vg_plot(j)=(g/(2*Dx*f))*(h(j+1,inew)-h(j-1,inew))
        END DO
        vg_plot(1)=v(1,inew)
        vg_plot(JM)=v(JM,inew)
        DO j=1,JM
           v_plot(j)=v(j,inew)
        END DO
     END IF

     IF(N.EQ.NMAX) THEN
        EXIT
     ELSE
        isave = iold
        iold = icur
        icur = inew
        inew = isave
     END IF

  END DO

  DO j=1,JM
     WRITE(UNIT=10,FMT='(1X,F10.4,2X,6F12.4)',IOSTAT=res) &
          X(j)/1000,h_plot(j,1),h_plot(j,2),h_plot(j,3),h_plot(j,4),h_plot(j,5),h_plot(j,6)
     WRITE(UNIT=11,FMT='(1X,F10.4,2X,4F12.4)',IOSTAT=res) &
          X(j)/1000,v_plot(j),vg_plot(j)
     IF(res /= 0) THEN                 ! An error has occured
        PRINT *, "Error in writing file, status: ", res
        EXIT                            ! Exit the loop
     END IF
  END DO

  CLOSE(UNIT=10)                      !Close the output file
  CLOSE(UNIT=11)                      !Close the output file

  STOP  'after successful run'        !End program in an orderly fashion
  

END PROGRAM SWE
