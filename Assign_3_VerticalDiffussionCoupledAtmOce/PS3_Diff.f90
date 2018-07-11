
PROGRAM Diff

implicit none


!define parameters and constants

INTEGER                ::   J                           !counter along z_axis  
INTEGER                ::   NTIME                       !counter along t_axis
INTEGER,PARAMETER      ::   JMAXa=29                    !maximun # of grid points along z_axis for atm
INTEGER,PARAMETER      ::   JMAXo=302                   !maximun # of grid points along z_axis for oce
INTEGER,PARAMETER      ::   JMAXt= JMAXo + JMAXa -1     !maximun # of grid points along z_axis for both atm and oce
INTEGER,PARAMETER      ::   COLS=2                      !columns to store INEW and ICUR Th values
INTEGER,PARAMETER      ::   PRINTS=6                    !columns for Th values to plot
INTEGER,PARAMETER      ::   JJa=JMAXa-1                 !Index # for last wet grid point along z-axis for atm (=JMAXa-1)
INTEGER,PARAMETER      ::   JJo=JMAXo-1                 !Index # for last wet grid point along z-axis for oce (=JMAXo-1) 
INTEGER,PARAMETER      ::   NTMAX=500001                !Maximum # of time levels NTMAX = 1201 
INTEGER,PARAMETER      ::   olun=11                     !Logical Unit Number (lun) for referencing the file holding the results (outfile)
INTEGER	               ::   N1	                        !1st time level for printout
INTEGER	               ::   N2	                        !2nd time level for printout
INTEGER	               ::   N3	                        !3rd time level for printout
INTEGER		       ::   N4	                        !4th time level for printout
INTEGER		       ::   N5	                        !5th time level for printout
INTEGER 	       ::   ISAVE                       !Auxiliary integer
INTEGER		       ::   INEW                        !New time step (n+1)
INTEGER	               ::   ICUR                        !Current time step (n)
INTEGER                ::   res                         !A status variable to hold the result of file operations
REAL                   ::   D=30.0                      !depth OML
REAL                   ::   H=270.0                     !height ABL
REAL                   ::   Th_btm=10.0                 !fixed temp at bottom OML (oC)
REAL                   ::   Th_top=0.0                  !fixed temp at top ABL (oC)
REAL                   ::   ka=30.0                     !Atm mixing coefficient (m².s⁻¹)
REAL                   ::   ko=0.003                    !Oce mixing coefficient (m².s⁻¹)
REAL                   ::   K_a=0.45
REAL                   ::   K_o=0.45
REAL                   ::   A                           !=(ko/ka)*(Dza/Dzo)


!define variables
REAL,DIMENSION(JMAXt)            :: Z
REAL,DIMENSION(JMAXt,COLS)       :: Th
REAL,DIMENSION(JMAXt,PRINTS)     :: ThP

!define time-space intervals

REAL                   ::   Dta        !time step for atm
REAL                   ::   Dto        !time step for oce
REAL                   ::   Dza        !increment along the z_axis for atm
REAL                   ::   Dzo        !increment along the z_axis for oce

!Declare character strings to hold the name of the file containing the results
CHARACTER(LEN=80)      ::   outfile    

!compute space increments Dz_a and Dz_o and location of vertical levels
   Dza=H/JJa
   Dzo=D/JJo

!compute time steps 
   Dta = K_a*Dza*Dza/ka	! Atmospheric application	
   Dto = K_o*Dzo*Dzo/ko	! Oceanic application
  
   A=(ko/ka)*(Dza/Dzo)

   outfile="diffusion.txt"

! Open the outfile and logfile and test the success
  OPEN(UNIT=olun,FILE=outfile,FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
    PRINT *, "Error in opening output file, status: ", res
    STOP					! Stop the program
  END IF       

!fill up height/depth vector

    DO J=1,JMAXo                     ! Ocean
      Z(J) = - D + (J-1)*Dzo
    ENDDO

!  Z(J_a=1)=Z(JMAX_o)                ! interfase

    DO J=JMAXo+1,JMAXt               ! Atmosphere
      Z(J) = ((J-1)-(JMAXo)+1)*Dza
    ENDDO  

!Initialize time counters and time
     ISAVE=0
     ICUR=1
     INEW=2
     NTIME=0

!Specify time levels for printout 

     !for K_o=K_a=0.45
     N1=50
     N2=500
     N3=5000
     N4=50000
     N5=500000

     !for K_o=K_a=0.55
     !N1=10
     !N2=20
     !N3=30
     !N4=40
     !N5=50

!initial condition
     
    DO J=1,JMAXt-1
      Th(J,ICUR)=Th_btm
      
      Th(J,INEW) = 0.0         !ensure that there are zeroes in the new temperature initially

      ThP(J,1)=Th(J,ICUR)      !extract the initial values for plotting
    ENDDO

      Th(JMAXt,ICUR)=Th_top

      ThP(JMAXt,1)=Th(JMAXt,ICUR)

!Calculate the interior temperature at the next time level
 
 DO NTIME=1,NTMAX

    DO J=2,JJo
      Th(J,INEW) = Th(J,ICUR) + K_o*(Th(J-1,ICUR) - 2.*Th(J,ICUR) + Th(J+1,ICUR))
    ENDDO
    
    DO J=JMAXo+1,JMAXt-1
      Th(J,INEW) = Th(J,ICUR) + K_a*(Th(J-1,ICUR) - 2.*Th(J,ICUR) + Th(J+1,ICUR))
    ENDDO

      Th(JMAXo,INEW)=(A*Th(JMAXo-2,INEW)+Th(JMAXo+2,INEW))/(A+1)     !interfase Fo=Fa

!boundary conditions

     Th(1,INEW)=Th_btm
     Th(JMAXt,INEW)=Th_top    !Th(JMAXt,INEW)=Th(JMAXt-1,INEW)  for the no flux case

!Plot time step
    IF(NTIME.EQ.N1) THEN
      DO J=1,JMAXt
        ThP(J,2)=Th(J,INEW)
      END DO
    ELSEIF(NTIME.EQ.N2) THEN
      DO J=1,JMAXt
        ThP(J,3)=Th(J,INEW)
      END DO
    ELSEIF(NTIME.EQ.N3) THEN
      DO J=1,JMAXt
        ThP(J,4)=Th(J,INEW)
      END DO
    ELSEIF(NTIME.EQ.N4) THEN
      DO J=1,JMAXt
        ThP(J,5)=Th(J,INEW)
      END DO
    ELSEIF(NTIME.EQ.N5) THEN
      DO J=1,JMAXt
        ThP(J,6)=Th(J,INEW)
      END DO
    ENDIF

!More time steps
    IF(NTIME.GE.NTMAX) THEN                   ! Jump out of the loop
      EXIT
    ELSE                                      ! Swap old and new time steps
      ISAVE=INEW
      INEW=ICUR
      ICUR=ISAVE
    ENDIF

 ENDDO

!// Write results to outfile and test whether the operation was successful
  DO J=1,JMAXt
    WRITE(UNIT=olun,FMT='(1X,F10.4,2X,6F12.4)',IOSTAT=res) &
      Z(J),ThP(J,1),ThP(J,2),ThP(J,3),ThP(J,4),ThP(J,5),ThP(J,6)
    IF(res /= 0) THEN                 ! An error has occured
      PRINT *, "Error in writing file, status: ", res
      EXIT                            ! Exit the loop
    END IF
  END DO 

!// Close the output file
  CLOSE(UNIT=olun)

!// End program in an orderly fashion
  STOP  'after successful run'

END PROGRAM Diff
