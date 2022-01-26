!=======================================================================
! Generated by  : PSCAD v4.6.3.0
!
! Warning:  The content of this file is automatically generated.
!           Do not modify, as any changes made here will be lost!
!-----------------------------------------------------------------------
! Component     : Main
! Description   : 
!-----------------------------------------------------------------------


!=======================================================================

      SUBROUTINE MainDyn()

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 'emtstor.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's2.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'fnames.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'matlab.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments

! Electrical Node Indices

! Control Signals
      INTEGER  sw, v_init
      REAL     var_L, i_L, v_L

! Internal Variables
      INTEGER  IVD1_1, IVD1_2
      REAL     RVD1_1, RVD1_2, RVD1_3, RVD1_4

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER ISTOI, ISTOF, IT_0                  ! Storage Indices
      INTEGER SS, INODE, IBRCH                    ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------

! Dsdyn <-> Dsout transfer index storage

      NTXFR = NTXFR + 1

      TXFR(NTXFR,1) = NSTOL
      TXFR(NTXFR,2) = NSTOI
      TXFR(NTXFR,3) = NSTOF
      TXFR(NTXFR,4) = NSTOC

! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      ISTOI     = NSTOI
      NSTOI     = NSTOI + 2
      ISTOF     = NSTOF
      NSTOF     = NSTOF + 3
      NPGB      = NPGB + 3
      INODE     = NNODE + 2
      NNODE     = NNODE + 6
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 6
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Transfers from storage arrays
!---------------------------------------

      var_L    = STOF(ISTOF + 1)
      i_L      = STOF(ISTOF + 2)
      v_L      = STOF(ISTOF + 3)
      sw       = STOI(ISTOI + 1)
      v_init   = STOI(ISTOI + 2)


!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


!---------------------------------------
! Configuration of Models
!---------------------------------------

      IF ( TIMEZERO ) THEN
         FILENAME = 'Main.dta'
         CALL EMTDC_OPENFILE
         SECTION = 'DATADSD:'
         CALL EMTDC_GOTOSECTION
      ENDIF
!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 10:[tbreakn] Timed Breaker Logic 
! Timed breaker logic
      IF ( TIMEZERO ) THEN
         sw = 1
      ELSE
         sw = 1
         IF ( TIME .GE. 0.0 ) sw = (1-1)
      ENDIF

! 20:[breaker1] Single Phase Breaker 'sw'
      IVD1_2 = NSTORI
      NSTORI = NSTORI + 1
      CALL E1PBRKR1_EXE(SS, (IBRCH+5),1.0e-09,1000000000.0,1,NINT(1.0-RE&
     &AL(sw)))
      IVD1_1 = 2*E_BtoI(OPENBR( (IBRCH+5),SS))
      IF (FIRSTSTEP .OR. (STORI(IVD1_2) .NE. IVD1_1)) THEN
         CALL PSCAD_AGI2(ICALL_NO,1210198031,IVD1_1,"BOpen")
      ENDIF
      STORI(IVD1_2) = 2*E_BtoI(OPENBR( (IBRCH+5),SS))

! 40:[tbreakn] Timed Breaker Logic 
! Timed breaker logic
      IF ( TIMEZERO ) THEN
         v_init = 0
      ELSE
         v_init = 0
         IF ( TIME .GE. 0.0 ) v_init = (1-0)
      ENDIF

! 90:[varrlc] Variable R, L or C  
      CALL COMPONENT_ID(ICALL_NO,1113893919)
      CALL E_VARRLC1_EXE(1 ,SS ,  (IBRCH+1), 0, var_L, 0.0)

! 100:[breaker1] Single Phase Breaker 'v_init'
      IVD1_2 = NSTORI
      NSTORI = NSTORI + 1
      CALL E1PBRKR1_EXE(SS, (IBRCH+2),1.0e-09,1000000000.0,1,NINT(1.0-RE&
     &AL(v_init)))
      IVD1_1 = 2*E_BtoI(OPENBR( (IBRCH+2),SS))
      IF (FIRSTSTEP .OR. (STORI(IVD1_2) .NE. IVD1_1)) THEN
         CALL PSCAD_AGI2(ICALL_NO,1679398847,IVD1_1,"BOpen")
      ENDIF
      STORI(IVD1_2) = 2*E_BtoI(OPENBR( (IBRCH+2),SS))

! 1:[source_1] Single Phase Voltage Source Model 2 'Source1'
! DC source with specified terminal conditions: Type: Ideal
      RVD1_1 = RTCF(NRTCF)
      RVD1_2 = RTCF(NRTCF+1)
      RVD1_3 = RTCF(NRTCF+2)
      RVD1_4 = RTCF(NRTCF+3)
      NRTCF = NRTCF + 4
      CALL EMTDC_1PVSRC(SS, (IBRCH+4),RVD1_4,0,RVD1_1,RVD1_2,RVD1_3)

!---------------------------------------
! Feedbacks and transfers to storage
!---------------------------------------

      STOF(ISTOF + 1) = var_L
      STOF(ISTOF + 2) = i_L
      STOF(ISTOF + 3) = v_L
      STOI(ISTOI + 1) = sw
      STOI(ISTOI + 2) = v_init


!---------------------------------------
! Transfer to Exports
!---------------------------------------

!---------------------------------------
! Close Model Data read
!---------------------------------------

      IF ( TIMEZERO ) CALL EMTDC_CLOSEFILE
      RETURN
      END

!=======================================================================

      SUBROUTINE MainOut()

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 'emtstor.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's2.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'fnames.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'matlab.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------

      REAL    VBRANCH       ! 
      REAL    EMTDC_VVDC    ! 

!---------------------------------------
! Variable Declarations
!---------------------------------------


! Electrical Node Indices
      INTEGER  NT_3

! Control Signals
      REAL     var_L, i_L, v_L

! Internal Variables
      INTEGER  IVD1_1

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER ISTOL, ISTOI, ISTOF, ISTOC, IT_0    ! Storage Indices
      INTEGER IPGB                                ! Control/Monitoring
      INTEGER SS, INODE, IBRCH                    ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------

! Dsdyn <-> Dsout transfer index storage

      NTXFR = NTXFR + 1

      ISTOL = TXFR(NTXFR,1)
      ISTOI = TXFR(NTXFR,2)
      ISTOF = TXFR(NTXFR,3)
      ISTOC = TXFR(NTXFR,4)

! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      IPGB      = NPGB
      NPGB      = NPGB + 3
      INODE     = NNODE + 2
      NNODE     = NNODE + 6
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 6
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Transfers from storage arrays
!---------------------------------------

      var_L    = STOF(ISTOF + 1)
      i_L      = STOF(ISTOF + 2)
      v_L      = STOF(ISTOF + 3)


!---------------------------------------
! Electrical Node Lookup
!---------------------------------------

      NT_3  = NODE(INODE + 3)

!---------------------------------------
! Configuration of Models
!---------------------------------------

      IF ( TIMEZERO ) THEN
         FILENAME = 'Main.dta'
         CALL EMTDC_OPENFILE
         SECTION = 'DATADSO:'
         CALL EMTDC_GOTOSECTION
      ENDIF
!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 20:[breaker1] Single Phase Breaker 'sw'
! Single phase breaker current
!

! 30:[multimeter] Multimeter 
      IVD1_1 = NRTCF
      NRTCF  = NRTCF + 5
      i_L = ( CBR((IBRCH+3), SS))
      v_L = EMTDC_VVDC(SS, NT_3, 0)

! 50:[xy_transfer_function] X-Y transfer function 
      CALL COMPONENT_ID(ICALL_NO,1026196919)
      CALL XYFUNC1_EXE(7,0,0,0.0,0.0,1.0,1.0,i_L,var_L)

! 60:[pgb] Output Channel 'inductance'

      PGB(IPGB+1) = var_L

! 70:[pgb] Output Channel 'current'

      PGB(IPGB+2) = i_L

! 80:[pgb] Output Channel 'voltage'

      PGB(IPGB+3) = v_L

! 100:[breaker1] Single Phase Breaker 'v_init'
! Single phase breaker current
!

!---------------------------------------
! Feedbacks and transfers to storage
!---------------------------------------

      STOF(ISTOF + 1) = var_L
      STOF(ISTOF + 2) = i_L
      STOF(ISTOF + 3) = v_L


!---------------------------------------
! Close Model Data read
!---------------------------------------

      IF ( TIMEZERO ) CALL EMTDC_CLOSEFILE
      RETURN
      END

!=======================================================================

      SUBROUTINE MainDyn_Begin()

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments

! Electrical Node Indices

! Control Signals

! Internal Variables

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER IT_0                                ! Storage Indices
      INTEGER SS, INODE, IBRCH                    ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------


! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      INODE     = NNODE + 2
      NNODE     = NNODE + 6
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 6
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Electrical Node Lookup
!---------------------------------------


!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 20:[breaker1] Single Phase Breaker 'sw'
      CALL COMPONENT_ID(ICALL_NO,1210198031)
      CALL E1PBRKR1_CFG(1.0e-09,1000000000.0,0.0)

! 90:[varrlc] Variable R, L or C  
      CALL E_VARRLC1_CFG(1 ,SS ,  (IBRCH+1), 0)

! 100:[breaker1] Single Phase Breaker 'v_init'
      CALL COMPONENT_ID(ICALL_NO,1679398847)
      CALL E1PBRKR1_CFG(1.0e-09,1000000000.0,0.0)

! 1:[source_1] Single Phase Voltage Source Model 2 'Source1'
      CALL E_1PVSRC_CFG(0,1,6,10.0,60.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)

      RETURN
      END

!=======================================================================

      SUBROUTINE MainOut_Begin()

!---------------------------------------
! Standard includes
!---------------------------------------

      INCLUDE 'nd.h'
      INCLUDE 'emtconst.h'
      INCLUDE 's0.h'
      INCLUDE 's1.h'
      INCLUDE 's4.h'
      INCLUDE 'branches.h'
      INCLUDE 'pscadv3.h'
      INCLUDE 'radiolinks.h'
      INCLUDE 'rtconfig.h'

!---------------------------------------
! Function/Subroutine Declarations
!---------------------------------------


!---------------------------------------
! Variable Declarations
!---------------------------------------


! Subroutine Arguments

! Electrical Node Indices
      INTEGER  NT_3

! Control Signals

! Internal Variables
      INTEGER  IVD1_1

! Indexing variables
      INTEGER ICALL_NO                            ! Module call num
      INTEGER IT_0                                ! Storage Indices
      INTEGER SS, INODE, IBRCH                    ! SS/Node/Branch/Xfmr


!---------------------------------------
! Local Indices
!---------------------------------------


! Define electric network subsystem number

      SS        = NODE(NNODE+1)

! Increment and assign runtime configuration call indices

      ICALL_NO  = NCALL_NO
      NCALL_NO  = NCALL_NO + 1

! Increment global storage indices

      INODE     = NNODE + 2
      NNODE     = NNODE + 6
      IBRCH     = NBRCH(SS)
      NBRCH(SS) = NBRCH(SS) + 6
      NCSCS     = NCSCS + 0
      NCSCR     = NCSCR + 0

!---------------------------------------
! Electrical Node Lookup
!---------------------------------------

      NT_3  = NODE(INODE + 3)

!---------------------------------------
! Generated code from module definition
!---------------------------------------


! 30:[multimeter] Multimeter 
      IVD1_1 = NRTCF
      NRTCF  = NRTCF + 5

! 50:[xy_transfer_function] X-Y transfer function 
      RTCF(NRTCF)    = -2.0
      RTCF(NRTCF+1)  = 0.0004444
      RTCF(NRTCF+2)  = -0.1
      RTCF(NRTCF+3)  = 0.0004444
      RTCF(NRTCF+4)  = -0.1
      RTCF(NRTCF+5)  = 0.02
      RTCF(NRTCF+6)  = 0.0
      RTCF(NRTCF+7)  = 0.02
      RTCF(NRTCF+8)  = 0.1
      RTCF(NRTCF+9)  = 0.02
      RTCF(NRTCF+10) = 0.1
      RTCF(NRTCF+11) = 0.0004444
      RTCF(NRTCF+12) = 2.0
      RTCF(NRTCF+13) = 0.0004444
      NRTCF = NRTCF + 14

! 60:[pgb] Output Channel 'inductance'

! 70:[pgb] Output Channel 'current'

! 80:[pgb] Output Channel 'voltage'

      RETURN
      END