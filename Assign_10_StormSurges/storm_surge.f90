PROGRAM storm_surge

  implicit none

  !
  REAL,PARAMETER                ::  g = 10.
  REAL,PARAMETER                ::  TAUsx = 0.
  REAL,PARAMETER                ::  TAUsy = 0.1
  REAL,PARAMETER                ::  dens = 1000.
  REAL,PARAMETER                ::  f = 0.0001
  REAL,PARAMETER                ::  Ho = 300.
  REAL,PARAMETER                ::  R = 0.0024
  REAL,PARAMETER                ::  Lr = sqrt(g*Ho)/f
  REAL,PARAMETER                ::  Uek = TAUsy/(dens*f)
  REAL,PARAMETER                ::  L = 10.*Lr
  INTEGER,PARAMETER             ::  JM = 101
  REAL,PARAMETER                ::  Dt = 600.
  REAL,PARAMETER                ::  Dx = L/(JM-1)
  INTEGER,PARAMETER             ::  NMAX = 630
  INTEGER,PARAMETER             ::  case = 2

  INTEGER                       ::  n,j
  INTEGER                       ::  res
  REAL                          ::  x
  REAL                          ::  t

  REAL,DIMENSION(JM)            ::  u
  REAL,DIMENSION(JM)            ::  v
  REAL,DIMENSION(JM)            ::  h
  REAL,DIMENSION(JM,NMAX+1)     ::  up
  REAL,DIMENSION(JM,NMAX+1)     ::  vp
  REAL,DIMENSION(JM,NMAX+1)     ::  hp


  OPEN(UNIT=11,FILE='up.txt',FORM="FORMATTED",IOSTAT=res)
  OPEN(UNIT=12,FILE='vp.txt',FORM="FORMATTED",IOSTAT=res)
  OPEN(UNIT=13,FILE='hp.txt',FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
     PRINT *, "Error in opening output file, status: ", res
     STOP					! Stop the program
  END IF

  !IC
  DO j=1,JM
     u(j)=0
     v(j)=0
     h(j)=Ho
     up(j,1)=u(j)
     vp(j,1)=v(j)
     hp(j,1)=h(j)
  END DO

  IF (case == 1) THEN
     DO n=2,NMAX+1
        
        DO j=1,JM
           v(j) = (1-R*Dt/Ho)*v(j) - f*Dt*u(j) + TAUsy*Dt/dens
        END DO

        DO j=2,JM-1
           u(j) = (1-R*Dt/Ho)*u(j) + f*Dt*v(j) - g*Ho*(Dt/Dx)*(h(j+1)-h(j))
        END DO

        u(JM)=0                                                                  !BC at coast (x=0)
        u(1)=(1-R*Dt/Ho)*u(1) + f*Dt*v(1) -g*Ho*(Dt/Dx)*(h(2)-h(1))                     !BC at open ocean (x=L)

        DO j=2,JM
           h(j) = h(j) - (Dt/Dx)*(u(j)-u(j-1))
        END DO

        h(1)=Ho                                                                  !BC at open ocean (x=L)

        DO j=1,JM
           up(j,n)=u(j)
           vp(j,n)=v(j)
           hp(j,n)=h(j)
        END DO

     END DO

  ELSE
     DO N=1,NMAX-1
        DO j=1,JM
           x=-1*L+(j-1)*Dx
           t=(n-1)*Dt

           CALL anal(u(j),v(j),h(j),Uek,Lr,Ho,f,x,t)
                up(j,n)=u(j)
                vp(j,n)=v(j)
                hp(j,n)=h(j)
        END DO
     END DO

  END IF

!!$  DO j=1,JM
     WRITE(UNIT=11,FMT='(1X,106F12.4)',IOSTAT=res),((up(j,n),n=1,NMAX+1,6),j=1,JM)
     WRITE(UNIT=12,FMT='(1X,106F12.4)',IOSTAT=res),((vp(j,n),n=1,NMAX+1,6),j=1,JM)
     WRITE(UNIT=13,FMT='(1X,106F12.4)',IOSTAT=res),((hp(j,n),n=1,NMAX+1,6),j=1,JM)
     IF(res /= 0) THEN                 ! An error has occured
        PRINT *, "Error in writing file, status: ", res
!!$        EXIT                            ! Exit the loop
     END IF
!!$  END DO

  CLOSE(UNIT=11)                      !Close the output file
  CLOSE(UNIT=12)                      !Close the output file
  CLOSE(UNIT=13)                      !Close the output file

  STOP  'after successful run'        !End program in an orderly fashion

END PROGRAM storm_surge


SUBROUTINE anal(uu,vv,hh,Uek,Lr,Ho,f,x,t)

implicit none

REAL                    ::  uu
REAL                    ::  vv
REAL                    ::  hh

REAL                    ::  Uek
REAL                    ::  f
REAL                    ::  Lr
REAL                    ::  Ho
REAL                    ::  x 
REAL                    ::  t

uu=Uek*(1-exp(x/Lr))
vv=f*t*Uek*exp(x/Lr)
hh=Ho*(1+t*Uek*exp(x/Lr)/(Lr*Ho))

END SUBROUTINE anal
