!/////////////////////////////////////////////////////////////////////////////////////////
!//                                                                                     //
!// This is the solution to Problem set 2 of GEF4510 FALL 2015                          //
!//-------------------------------------------------------------------------------------//
!// Revision history:									//
!//   - 2013-09-13: 	Reprogrammed by Lars Petter Røed: MET Norway based on earlier 	//
!//			Fortran 77 versions						//
!//   - 2013-09-24: 	Completed						 	//
!//   - 2015-09-15: 	Added revision history, and changed the explanation 	 	//
!//                                                                                     //
!//-------------------------------------------------------------------------------------//
!// Purpose:                                                                            //
!//                                                                                     //
!// To solve the diffusion equation in the form                                         //
!//                                                                                     //
!//	d(theta)/dt=k*d[d(theta)/dz]/dz							//
!//                                                                                     //
!// for an atmospheric application (A) for which theta (herafter th) is the potential	//
!// temperature, 0<z<D, and for an oceanic application (O) for which -D<z<0. D is	//
!// height in the atmospheric application and depth in the oceanic application.		//
!//                                                                                     //
!//-------------------------------------------------------------------------------------//
!// Method:                                                                             //
!//                                                                                     //
!// An explicit forward in time, centered in space scheme is used, that is,		//
!//                                                                                     //
!//	th_{j,n+1}=th_{j,n}+K*(th_{j+1,n}-2th_{j,n}+th_{j-1,n})				//
!//                                                                                     //
!// (latex notation) where								//
!//                                                                                     //
!//	j = 2(1)JMAX-1 and n = 0(1)NTMAX.						//
!//                                                                                     //
!// K is given by									//
!//                                                                                     //
!//	K=(kappa*Dt)/(2*Dz),								//
!//                                                                                     //
!// where Dt is the time step, Dz is the space increment, and kappa is the mixing	//
!// coefficient. The total number of time levels are specified in NTMAX, and maximum	//
!// number of space points is set to JMAX = 27.						//
!//                                                                                     //
!// Regarding the atmospheric application we are asked to compute the temperature for	//
!// two values of the constant K, namely K = 0.45 and K = 0.55 given the boundary       //
!// conditions:                                                                         //   
!//                                                                                     //
!//	th_{1,n} = 0, th_{JMAX,n} = 0 ; n = 0(1)N					//
!//                                                                                     //
!// and initial condition:                                                              //
!//                                                                                     //
!//       th_{j,0} = th0*sin(pi*z/D)                                                    //
!//                                                                                     //
!// where th0=10 degrees Celsius.                                                       // 
!//                                                                                     //
!// Regarding the oceanic application we are asked to compute th, with K=0.45, given    // 
!// the boundary conditions:                                                            // 
!//                                                                                     //
!//       th_{1,n} = 0,                                                                 //
!//       th_{JMAX,n} = th0*n*Dt/tc  for n*Dt < tc,                                     //
!//       th_{JMAX,n} = th0          for n*Dt >= tc,                                    //
!//                                                                                     //
!// where tc = 6 days, and the inititial condition                                      //
!//                                                                                     //
!//       th_{j,0} = 0.                                                                 //
!//                                                                                     //
!// We may also choose to use the boundary conditions                                   //
!//                                                                                     //
!//       th_{1,n} = 0,                                                                 // 
!//       th_{JMAX,n} = th0*tanh(gamma*(n*Dt - tc)/tc) for all n,                       //
!//                                                                                     //
!// where gamma = 1.5.                                                                  //
!//                                                                                     //
!/////////////////////////////////////////////////////////////////////////////////////////
!
PROGRAM diffu
  IMPLICIT NONE
!//-------------------------------------------------------------------------------------//
!// Specify dimensions and declare variables
!
!// Fixed parameters and constants
  INTEGER,PARAMETER		:: NTMAX=1201	! Maximum # of time levels NTMAX = 1201 
                                            	! => 30 days 
  INTEGER,PARAMETER		:: NPRINT=5	! Total number of printouts
  INTEGER,PARAMETER		:: NM=2		! Number of time levels saved at any time
  INTEGER,PARAMETER		:: olun=11	! Logical Unit Number (lun) for referencing
						! the file holding the results (outfile)
  INTEGER,PARAMETER		:: llun=12	! Logical Unit Number (lun) for referencing
						! the logrun file (logfile)
  INTEGER,PARAMETER		:: JMAX=27	! Maximum # of grid points along z-axis
  REAL				:: D=100.	! Depth (or height) of mixed layer
  REAL				:: th0=10.	! Maximum temperature at the boundary for 
						! forced case and maximum temperature 
						! anomaly for unforced case (unit oC)
  REAL				:: TDA=86400.	! Number of seconds in a day
  REAL				:: TC=518400.	! Number of seconds in 6 days. Used to de-
						! termine how fast the surface temperature  
                                      	    	! reaches its final value in the forced case
  REAL				:: ARG0=1.5	! Dimensionless constant used to determine 
						! how fast (in time) the hyperbolic tangent
						! approaches the maximum value.
  REAL				:: KOP=0.003	! Oceanic mixing coefficient (m^2/s)
  REAL				:: KAP=30.	! Atmospheric mixing coefficient (m^2/s)
!
!// Integer variables
  INTEGER			:: N1		! 1st time level for printout
  INTEGER			:: N2		! 2nd time level for printout
  INTEGER			:: N3		! 3rd time level for printout
  INTEGER			:: N4		! 4th time level for printout
  INTEGER			:: NTIME		! Time step counters
  INTEGER			:: J		! Counter along z-axis  
  INTEGER			:: JJ		! Index # for last wet grid point along  
						! z-axis (=JMAX-1)
  INTEGER 			:: ISAVE		! Auxiliary integer
  INTEGER			:: INEW		! New time step (n+1)
  INTEGER			:: ICUR		! Current time step (n)
!
!// Real variables
  REAL,DIMENSION(JMAX)        :: Z		! Depth (or height) coordinate (unit m)
  REAL,DIMENSION(JMAX,NM)     :: th		! Temperature variable (deg C) 
  REAL,DIMENSION(JMAX,NPRINT) :: thP		! Temperatures for printout (deg C)
!
!// Real constants
  REAL                        :: ARG        ! Auxiliary argument 
  REAL                        :: PI         ! The mathematical constant pi
  REAL                        :: ZZ         ! Dimensionless vertical coordinate along 
                                            ! z-axis (=z/D)
  REAL                        :: DZ         ! Increment along z-axis (m)
  REAL                        :: DT         ! Time step (s)
  REAL                        :: K          ! Dimensionless auxiliary constant
                                            ! (=0.45 or =0.55)
  REAL                        :: T1         ! 1st time for printout (unit s)
  REAL                        :: T2         ! 2nd time for printout (unit s)
  REAL                        :: T3         ! 3rd time for printout (unit s)
  REAL                        :: T4         ! 4th time for printout (unit s)
  REAL                        :: TIME       ! Time in seconds
  REAL                        :: THOURS     ! Time in hours
  REAL                        :: TDAYS      ! Time in days
!
!//Flags
  INTEGER                     :: IAO        ! Discriminates between the atmospheric
                                            ! application (= 1) or the oceanic 
                                            ! application (= 2)
  INTEGER                     :: ISTAB      ! Discriminates between the stable case (=1)  
                                            ! and unstable case (=2)
  INTEGER                     :: ITANH      ! Discriminates between a linear (=1)   
                                            ! and a tanh case (=2) in the ocean
                                            ! application
!
!// Character strings to hold the prompts for communicating with the user
  CHARACTER(LEN=80)           :: prompt1, prompt2, prompt3
!
!// A single character to hold the answers which is either A (atmosphere) or O (ocean)
!// for answer1, S (stable) or U (unstable) for answer2, and either N (no) or Y (yes)
!// for answer3 
  CHARACTER                   :: answer1,answer2,answer3
!
!// Declare character strings to hold the name of the file containing the results
  CHARACTER(LEN=80)           :: outfile     ! Holds the results
!
!// Declare character strings to hold the name of the file containing the logrun
  CHARACTER(LEN=80)           :: logfile     ! Holds the logrun info
!
!// A status variable to hold the result of file operations
  INTEGER :: res
!
!// Assign values to the prompts
  prompt1 = 'Specify application (A: atmosphere, O: ocean)'
  prompt2 = 'Specify case (S: Stable, K=0.45, U: Unstable, K=0.55)'
  prompt3 = 'Specify linear forcing or not (Y: Yes, N: No)'
!
!//-------------------------------------------------------------------------------------//
!// Assign values to constants, flags, etc.:
  JJ=JMAX-1
  PI=4.*ATAN(1.)
!
!// Ask whether this is an atmospheric application or an oceanic application, read 
!// the answer (answer1) from the keyboard input, and assign values to flag IAO
  PRINT *, TRIM(prompt1)
  READ(*,*) answer1
  IF(answer1.EQ.'A') THEN               ! Atmospheric application
    IAO = 1
  ELSEIF(answer1.EQ.'O') THEN           ! Ocean application
    IAO = 2
  ELSE                                  ! Default = Atmosphere
    PRINT *, 'Error in input! Using default = atmosphere application'
    IAO = 1
    answer1 = 'A'
  END IF
!
!// Next ask whether this is a unforced case or a forced case, read in the answer 
!// (answer2) from the keybord input, and assign values to flag ISTAB
  PRINT *, TRIM(prompt2)
  READ(*,*) answer2
  IF(answer2.EQ.'S') THEN               ! Stable case    
    ISTAB = 1
  ELSEIF(answer2.EQ.'U') THEN           ! Unstable case
    ISTAB = 2
  ELSE                                  ! Default = Atmosphere
    PRINT *, 'Error in input! Using default = stable case'
    ISTAB = 1
    answer2 = 'S'
  END IF
!
!// Next ask whether this is a linear forcing case or not, read in the answer 
!// (answer3) from the keybord input, and assign values to flag ITANH. Only 
  IF(IAO.EQ.2) THEN                          ! Ocean application                         
    PRINT *, TRIM(prompt3)
    READ(*,*) answer3
    IF(answer3.EQ.'Y') THEN                  ! Linear increase in boundary forcing
      ITANH = 1
    ELSEIF(answer3.EQ.'N') THEN              ! tanh increase in boundary forcing 
      ITANH = 2
    ELSE                                     ! Default = Linear
      PRINT *, 'Error in input! Using default = Linear forcing'
      ITANH = 1
      answer3 = 'Y'
    END IF
  ENDIF              
!
!//-------------------------------------------------------------------------------------//
!// Assign a meaningful filename to the outfile and logfile
  IF(IAO.EQ.1) THEN				! Atmospheric application
    IF(ISTAB.EQ.1) THEN				! Stable case
      outfile = "AtmosStable.txt"
      logfile = "AtmosStable_log.txt"
    ELSE					! Unstable case
      outfile = "AtmosUnstable.txt" 
      logfile = "AtmosUnstable_log.txt" 
    ENDIF
  ELSE						! Oceanic application
    IF(ISTAB.EQ.1) THEN				! Stable case
      IF(ITANH.EQ.1) THEN			! Linear increase in time
        outfile = "OceanStableLin.txt"
        logfile = "OceanStableLin_log.txt"
      ELSE					! tanh increase in time
        outfile = "OceanStableTanh.txt"
        logfile = "OceanStableTanh_log.txt"
      ENDIF
    ELSE					! Unstable case
      IF(ITANH.EQ.1) THEN			! Linear increase in time
        outfile = "OceanUnstableLin.txt"
        logfile = "OceanUnstableLin_log.txt"
      ELSE					! tanh increase in time
        outfile = "OceanUnstableTanh.txt"
        logfile = "OceanUnstableTanh_log.txt"
      ENDIF
    ENDIF
  ENDIF   
!
!// Open the outfile and logfile and test the success
  OPEN(UNIT=olun,FILE=outfile,FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
    PRINT *, "Error in opening output file, status: ", res
    STOP					! Stop the program
  END IF       
!
  OPEN(UNIT=llun,FILE=logfile,FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
    PRINT *, "Error in opening output file, status: ", res
    STOP					! Stop the program
  END IF       
!
!//-------------------------------------------------------------------------------------//
!// Compute space increment DZ and location of vertical levels
  DZ=D/JJ
  IF(IAO.EQ.1) THEN	! Atmosphere
    DO J=1,JMAX
      Z(J) = (J-1)*DZ
    ENDDO
  ELSE			! Ocean
    DO J=1,JMAX
      Z(J) = - D + (J-1)*DZ
    ENDDO
  ENDIF
!
!//-------------------------------------------------------------------------------------//
!// Specify K and compute time step DT 
  IF(ISTAB.EQ.1) THEN	! Stable
    K=0.45
  ELSE			! Unstable
    K=0.55
  ENDIF
  IF(IAO.EQ.1) THEN	! Atmospheric application
    DT = K*DZ*DZ/KAP
  ELSE			! Oceanic application
    DT = K*DZ*DZ/KOP
  ENDIF
!
!// Print application, mixing coefficient, the dimensionless constant K, forcing case
!// and time step to the terminal window
  PRINT *, answer1 
  PRINT *, answer2 
  IF(IAO.EQ.2) THEN	! Ocean application
    PRINT *, answer3
  ENDIF
  IF(IAO.EQ.1) THEN	! Atmosphere
    PRINT *, 'KAP: ', KAP,'m^2/s'
  ELSE			! Ocean
    PRINT *, 'KOP: ', KOP,'m^2/s'
  ENDIF
  PRINT *, 'K  = ', K
  PRINT *, 'DT = ', DT,'s'
!
!// Specify time levels for printout (< NTNAX-1)
  IF(ISTAB.EQ.1) THEN	! Stable case
    N1=100
    N2=200
    N3=600
    N4=1200
  ELSE			! Unstable
    N1=30
    N2=60
    N3=90
    N4=100
  ENDIF
!
!// Convert time levels for printout to time in minutes for the atmospheric
!// application and days for the oceanic application
  IF(IAO.EQ.1) THEN	! Atmospheric application
    T1=N1*DT/60.
    T2=N2*DT/60.
    T3=N3*DT/60.
    T4=N4*DT/60.
  ELSE			! Oceanic application
    T1=N1*DT/TDA
    T2=N2*DT/TDA
    T3=N3*DT/TDA
    T4=N4*DT/TDA
  ENDIF    
!
!// Write runtime variables to logfile and test whether the operation was successful
  WRITE(UNIT=llun,FMT='(1X,A49,1X,A1)',IOSTAT=res) &
    'Application (A: Atmospheric, O: Oceanic)        : ', answer1
  WRITE(UNIT=llun,FMT='(1X,A49,1X,A1)',IOSTAT=res) &
    'Stable or unstable case (S: Stable, U: Unstable): ', answer2
  IF(IAO.EQ.2) THEN               ! Ocean application
    WRITE(UNIT=llun,FMT='(1X,A49,1X,A1)',IOSTAT=res) &
    'Linear or tanh forcing (Y: Linear, N: tanh)?    : ', answer3
  ENDIF
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A5)',IOSTAT=res) &
    'Mixing coefficient atmosphere                   : ', KAP,'m^2/s'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A5)',IOSTAT=res) &
    'Mixing coefficient ocean                        : ', KOP,'m^2/s'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A3)',IOSTAT=res) &
    '1st time level for plot after                   : ', T1, 'min'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A3)',IOSTAT=res) &
    '2nd time level for plot after                   : ', T2, 'min'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A3)',IOSTAT=res) &
    '3rd time level for plot after                   : ', T3, 'min'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A3)',IOSTAT=res) &
    '4th time level for plot after                   : ', T4, 'min'
  WRITE(UNIT=llun,FMT='(1X,A49,1X,F9.4,1X,A1)',IOSTAT=res) &
    'Time step DT                                    : ', DT,'s'
  IF(res /= 0) THEN                 ! An error has occured
    PRINT *, "Error in writing file, status: ", res
    STOP                            ! Stop the program
  END IF 
!
!//-------------------------------------------------------------------------------------//
!// Initialize time counters and time
  ISAVE=0
  ICUR=1
  INEW=2
  NTIME=0
  TIME=0.
  THOURS=0.
  TDAYS=0.
!
!// Initial temperature distribution
  DO J=1,JMAX
    IF(IAO.EQ.1) THEN                         ! Atmospheric application
      th(J,ICUR) = th0*SIN(PI*Z(J)/D)         
    ELSE                                      ! Oceanic application
      th(J,ICUR)=0.0                          
    ENDIF
!
!// Ensure that there are zeroes in the new temperature initially
    th(J,INEW) = 0.0
!
!// Extract the initial values for plotting
    thP(J,1)=th(J,ICUR)/th0
  END DO
!//-------------------------------------------------------------------------------------//
!// Start integration
  DO NTIME=1,NTMAX                            ! Start time loop
    TIME=NTIME*DT
    THOURS=TIME/3600.
    TDAYS=THOURS/24.
!
!// Calculate the interior temperature at the next time level
    DO J=2,JJ
      th(J,INEW) = th(J,ICUR) + K*(th(J-1,ICUR) - 2.*th(J,ICUR) + th(J+1,ICUR))
    END DO
!
!// Then at boundaries:
    th(1,INEW) = 0.0
    IF(IAO.EQ.1) THEN                         ! Atmospheric application
      th(JMAX,INEW) = 0.0
    ELSE                                      ! Oceanic application
      IF(ITANH.EQ.1) THEN
        IF(TIME.LT.TC) THEN
          th(JMAX,INEW) = th0*TIME/TC
        ELSE
          th(JMAX,INEW) = th0
        ENDIF
      ELSE 
        ARG=ARG0*TIME/TC
        th(JMAX,INEW) = th0*TANH(ARG)
      ENDIF
    ENDIF  
!
!// Plot this time step?
    IF(NTIME.EQ.N1) THEN
      DO J=1,JMAX
        thP(J,2)=th(J,INEW)/th0
      END DO
    ELSEIF(NTIME.EQ.N2) THEN
      DO J=1,JMAX
        thP(J,3)=th(J,INEW)/th0
      END DO
    ELSEIF(NTIME.EQ.N3) THEN
      DO J=1,JMAX
        thP(J,4)=th(J,INEW)/th0
      END DO
    ELSEIF(NTIME.EQ.N4) THEN
      DO J=1,JMAX
        thP(J,5)=th(J,INEW)/th0
      END DO
    ENDIF
!
!// More time steps?
    IF(NTIME.GE.NTMAX) THEN                   ! Jump out of the loop
      EXIT
    ELSE                                      ! Swap old and new time steps
      ISAVE=INEW
      INEW=ICUR
      ICUR=ISAVE
    ENDIF
!
!// Do next integration step
  END DO                                      ! End time loop
!
!//-------------------------------------------------------------------------------------//
!// Write results to outfile and test whether the operation was successful
  DO J=1,JMAX
    ZZ = Z(J)/D
    WRITE(UNIT=olun,FMT='(1X,F10.4,2X,5F12.4)',IOSTAT=res) &
      ZZ,thP(J,1),thP(J,2),thP(J,3),thP(J,4),thP(J,5)
    IF(res /= 0) THEN                 ! An error has occured
      PRINT *, "Error in writing file, status: ", res
      EXIT                            ! Exit the loop
    END IF
  END DO 
!
!//-------------------------------------------------------------------------------------//
!// Close the output file
  CLOSE(UNIT=olun)
  CLOSE(UNIT=llun)
!
!// End program in an orderly fashion
  STOP  'after successful run'
END PROGRAM diffu



