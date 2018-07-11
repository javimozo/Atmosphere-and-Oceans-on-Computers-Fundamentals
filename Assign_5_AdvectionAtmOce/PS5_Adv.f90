

PROGRAM Adv

  implicit none

  REAL,PARAMETER                     ::  L=50.0                         !domain x=[-L,L]
  REAL                               ::  sigma                          !width of gaussian
  REAL                               ::  C                              !Courant number
  REAL,PARAMETER                     ::  Dx=1        !0.1*sigma         !space step
  INTEGER,PARAMETER                  ::  J=2*L/Dx                       !number of grid points
  REAL                               ::  num_Dt                         !number of time steps within a cicle
  REAL,DIMENSION(J+1)                ::  X                              !matrix with space domain with J points
  REAL,ALLOCATABLE,DIMENSION(:,:)    ::  PHI                            !matrix for the function
  INTEGER                            ::  num_Cy                         !number of cycles
  INTEGER                            ::  jj,n                           !counters
  INTEGER                            ::  ios                            !
  INTEGER                            ::  olun                           !Logical Unit Number (lun) for referencing the file holding the results (outfile)
  INTEGER                            ::  res                            !A status variable to hold the result of file operations
  CHARACTER(LEN=80)                  ::  answer1,answer2,answer3        !to questions about C, sigma and scheme
  CHARACTER(LEN=80)                  ::  prompt1,prompt2,prompt3        !questions about C, sigma and scheme
  CHARACTER(LEN=80)                  ::  outfile                        !file to save the results

  prompt1 = 'Specify C (A: C=1, B: C=0.5)'
  prompt2 = 'Specify width of gaussian (W: wide, N: narrow)'
  prompt3 = 'Specify scheme (LF: leapfrog, UP: upwind/upstream, LW: Lax-Wendroff)'

  PRINT *, TRIM(prompt1)               !select value of C
  READ(*,*) answer1
  IF(answer1.EQ.'A') THEN              
    C = 1
  ELSEIF(answer1.EQ.'B') THEN           
    C = 0.5
  ELSE                                  
    PRINT *, 'Error in input! Using default: C=1'
    C = 1
    answer1 = 'A'
  END IF

  PRINT *, TRIM(prompt2)               !select width of gaussian (value of sigma)
  READ(*,*) answer2
  IF(answer2.EQ.'W') THEN             
    sigma = 2*L/10
  ELSEIF(answer2.EQ.'N') THEN          
    sigma = 2*L/1000
  ELSE                                 
    PRINT *, 'Error in input! Using default: wide'
    sigma = 2*L/10
    answer2 = 'W'
  END IF

  num_Dt=2*L/(Dx*C)                        !o tambien = J/C
  num_Cy=nint(10*num_Dt)
  ALLOCATE(PHI(J+1,num_Cy+1))

  DO jj=1,J+1                              !fill up matrix for domain with distance for each point
     X(jj)=-L+(jj-1)*Dx
  END DO

  DO jj=1,J+1                              !initial conditions _ 1st column
     PHI(jj,1)=exp(-(X(jj)/sigma)**2)
  END DO

  outfile="Adv.txt"

  OPEN(UNIT=olun,FILE=outfile,FORM="FORMATTED",IOSTAT=res)
  IF(res /= 0) THEN				! An error has occured
     PRINT *, "Error in opening output file, status: ", res
     STOP					! Stop the program
  END IF

  PRINT *, TRIM(prompt3)                  !select scheme to use
  READ(*,*) answer3
  IF(answer3.EQ.'LF') THEN              
  
     !LEAPFROG  ----------------------------------------------------------------------------------
     !because this is a 3-level scheme we cannot compute directly Phi at t=2 (we lack info at t=0) 
     !so we approximate Phi at that time step only through the FTCS scheme, from t=1 to t=2
     
     DO jj=2,J
        PHI(jj,2)=PHI(jj,1)-0.5*C*(PHI(jj+1,1)-PHI(jj-1,1))   !1st cycle _ 2nd column     
     END DO

     PHI(1,2)=PHI(1,1)+0.5*C*(PHI(2,1)-PHI(J,1))                                         
     PHI(J+1,2)=PHI(1,2)                                      !Boundary Condition

     !for phi at t=3 onwards we procede with the CTCS scheme
  
     DO n=2,num_Cy                                            !2nd cycle _ 3rd column and onwards

        DO jj=2,J
           PHI(jj,n+1)=PHI(jj,n-1)-C*(PHI(jj+1,n)-PHI(jj-1,n))
        END DO

        PHI(1,n+1)=PHI(1,n-1)-C*(PHI(2,n)-PHI(J,n))
        PHI(J+1,n+1)=PHI(1,n+1)                               !Boundary Condition

     END DO

  ELSEIF(answer3.EQ.'UP') THEN         
  
     !UPSTREAM/UPWIND  ------------------------------------------------------------------------------
  
     DO n=1,num_Cy

        DO jj=2,J+1
           PHI(jj,n+1)=PHI(jj,n)-C*(PHI(jj,n)-PHI(jj-1,n))
        END DO

        PHI(1,n+1)=PHI(J+1,n+1)                               !Boundary Condition

     END DO
  
  ELSEIF(answer3.EQ.'LW') THEN                        

     !LAX-WENDROFF -----------------------------------------------------------------------------------
     
     DO n=1,num_Cy

        DO jj=2,J
           PHI(jj,n+1)=PHI(jj,n)-0.5*C*(PHI(jj+1,n)-PHI(jj-1,n))+0.5*(C**2)*(PHI(jj+1,n)-2*PHI(jj,n)+PHI(jj-1,n))
        END DO

        PHI(1,n+1)=PHI(1,n)-0.5*C*(PHI(2,n)-PHI(J,n))+0.5*(C**2)*(PHI(2,n)-2*PHI(1,n)+PHI(J,n))
        PHI(J+1,n+1)=PHI(1,n+1)                                                                           !Boundary Condition

     END DO

  ELSE                                  
    PRINT *, 'Error in input! Using default: upwind'
        
     DO n=1,num_Cy

         DO jj=2,J+1
           PHI(jj,n+1)=PHI(jj,n)-C*(PHI(jj,n)-PHI(jj-1,n))
         END DO

         PHI(1,n+1)=PHI(J+1,n+1)                               !Boundary Condition

     END DO

  END IF

!print results

  DO jj=1,J+1
     WRITE(UNIT=olun,FMT='(1X,F10.4,2X,6F12.4)',IOSTAT=res) &
          X(jj)/L,PHI(jj,1),PHI(jj,1+nint(num_Dt/2)),PHI(jj,1+nint(num_Dt)),PHI(jj,1+nint(2*num_Dt)), &
          PHI(jj,1+nint(5*num_Dt)),PHI(jj,1+nint(10*num_Dt))
     IF(res /= 0) THEN                 ! An error has occured
        PRINT *, "Error in writing file, status: ", res
        EXIT                            ! Exit the loop
     END IF
  END DO

  CLOSE(UNIT=olun)                    !Close the output file

  STOP  'after successful run'        !End program in an orderly fashion

END PROGRAM Adv
