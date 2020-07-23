  !
  !******************************************************************************
  !       The Coupled Routing and Excelss Storage (CREST) Version 2.1(10-01-2014)
  !                Developed by Dr. Xianwu Xue and Dr. Yang Hong
  !
  !  		    Previous versions(e.g. v1.6c) were developed
  !	                by Dr. Jiahu Wang and Dr. Yang Hong
  !                  Copyright (c) 2014 University of Oklahoma
  !                           All rights reserved
  !       THE UNIVERSITY OF OKLAHOMA MAKE NO REPRESENTATION OR WARRANTIES
  !       WITH RESPECT TO THE CONTENTS HEROF AND SPECIFICALLY DISCLAIM ANY
  !       IMPLIED WARRANTIES OR MERCHANTABILITY OR FITNESS FOR ANY
  !       PARTICULAR PURPOSE.
  !
  !       Further, we reserve the right to revise this software and/or
  !       Documentation and to make changes from time to time in the
  !       content without obligation to notify any person of such revision
  !       or change.
  !
  !       Please cite the following reference if you use CREST model:
  !       Wang, J., H. Yang, L. Li, et al. (2011), The coupled routing
  !          and excess storage(CREST) distributed hydrological model,
  !          Hydrological Sciences Journal, 56(1), 84 - 98.
  !       Xue, X., Y. Hong, A. S. Limaye, et al. (2013), Statistical and hydrological evaluation
  !          of TRMM-based Multi-satellite Precipitation Analysis over the Wangchu Basin of Bhutan:
  !          Are the latest satellite precipitation products 3B42V7 ready for use in ungauged basins?,
  !          Journal of Hydrology, 499, 91-99.
  !
  !       If you have any problem with this model, please e-mail:
  !       Dr. Xianwu Xue (xuexianwu@ou.edu) or Dr. Yang Hong (yanghong@ou.edu)
  !
  !******************************************************************************

!============================================================
!
!============================================================
  MODULE CREST_Project
  IMPLICIT NONE
  INTEGER :: IsDebug=1
  CHARACTER(LEN=14), PARAMETER:: g_CREST_Version = "2.1 10/31/2014"
  CHARACTER (LEN=19) :: g_strD

  CHARACTER(LEN=200):: g_ProjectName,g_ProjectWS
  CHARACTER(LEN=200):: g_PrjNP
  INTEGER :: g_CREST_LogFileID

  CHARACTER(LEN=1) :: g_TimeMark !y(year);m(month);d(day);h(hour);u(minute);s(second)
  INTEGER :: g_TimeStep
  INTEGER :: g_StartDate(1:6)

  INTEGER :: g_NumWarmup,g_ITMax

  INTEGER :: g_WarmupDate(1:6)
  INTEGER :: g_EndDate(1:6)
  LOGICAL :: g_LoadState
  LOGICAL :: g_SaveState

  CHARACTER(LEN=20):: g_sRunStyle
  INTEGER :: g_RunStyle ! simu--1; cali--2; REAL--3; repe --4; fore ---6
  CHARACTER(LEN=10) :: g_BasicFormat
  CHARACTER(LEN=200):: g_BasicPath

  CHARACTER(LEN=10) :: g_ParamFormat
  CHARACTER(LEN=200):: g_ParamPath

  CHARACTER(LEN=10) :: g_StateFormat
  CHARACTER(LEN=200):: g_StatePath

  CHARACTER(LEN=10) :: g_ICSFormat
  CHARACTER(LEN=200):: g_ICSPath

  CHARACTER(LEN=10) :: g_OBSFormat !Leave it to USE in furture
  CHARACTER(LEN=200):: g_OBSPath

  CHARACTER(LEN=10) :: g_RainFormat
  CHARACTER(LEN=200):: g_RainPath

  CHARACTER(LEN=10) :: g_PETFormat
  CHARACTER(LEN=200):: g_PETPath

  CHARACTER(LEN=10) :: g_ResultFormat
  CHARACTER(LEN=200):: g_ResultPath

  CHARACTER(LEN=10) :: g_CalibFormat
  CHARACTER(LEN=200):: g_CalibPath


  INTEGER :: g_NOutPixs

  CHARACTER(LEN=20) :: g_RoutingType ="JLR" !JLR-Jumped Linear Reservoir; CLR-Continous Linear Reservoir

  TYPE:: CREST_OutPix
    CHARACTER(LEN=200) :: Name
    INTEGER :: Col
    INTEGER :: Row
    LOGICAL :: bIsOut
    INTEGER, ALLOCATABLE :: Mask(:,:)

    DOUBLE PRECISION, ALLOCATABLE :: Rain(:)
    DOUBLE PRECISION, ALLOCATABLE :: PET(:)
    DOUBLE PRECISION, ALLOCATABLE :: EPot(:)
    DOUBLE PRECISION, ALLOCATABLE :: EAct(:)
    DOUBLE PRECISION, ALLOCATABLE :: W(:)
    DOUBLE PRECISION, ALLOCATABLE :: SM(:) ! Soil Moisture
    DOUBLE PRECISION, ALLOCATABLE :: R(:)
    DOUBLE PRECISION, ALLOCATABLE :: ExcS(:),ExcI(:)
    DOUBLE PRECISION, ALLOCATABLE :: RS(:),RI(:)
    DOUBLE PRECISION, ALLOCATABLE :: RObs(:)
    LOGICAL :: HasOBS
  END TYPE CREST_OutPix

  TYPE(CREST_OutPix), ALLOCATABLE :: g_tOutPix(:)


  TYPE:: CREST_Outlet
    CHARACTER(LEN=200) :: Name
    INTEGER :: Col
    INTEGER :: Row
    LOGICAL :: bIsOut
    INTEGER, ALLOCATABLE :: Mask(:,:)

    DOUBLE PRECISION, ALLOCATABLE :: Rain(:)
    DOUBLE PRECISION, ALLOCATABLE :: PET(:)
    DOUBLE PRECISION, ALLOCATABLE :: EPot(:)
    DOUBLE PRECISION, ALLOCATABLE :: EAct(:)
    DOUBLE PRECISION, ALLOCATABLE :: W(:)
    DOUBLE PRECISION, ALLOCATABLE :: SM(:) ! Soil Moisture
    DOUBLE PRECISION, ALLOCATABLE :: R(:)
    DOUBLE PRECISION, ALLOCATABLE :: ExcS(:),ExcI(:)
    DOUBLE PRECISION, ALLOCATABLE :: RS(:),RI(:)

    DOUBLE PRECISION, ALLOCATABLE :: RObs(:)
    LOGICAL :: HasOutlet
    LOGICAL :: HasOBS
  END TYPE CREST_Outlet

  TYPE(CREST_Outlet) :: g_tOutlet

  INTEGER :: g_NOutDTs
  INTEGER, ALLOCATABLE :: g_OutDTIndex(:)

  CHARACTER(LEN=5)::g_sGOVarName(0:10) !Grid Output Variable Name
  DATA g_sGOVarName /"Rain","PET","EPot","EAct","W","SM", &
    "R","ExcS","ExcI","RS","RI"/
  LOGICAL, ALLOCATABLE :: g_bGOVar(:)

  !For Lake======================================
  INTEGER :: g_NLakes
  INTEGER, ALLOCATABLE :: g_LakeMask(:,:)


  ! User Defined Functions
  INTEGER, EXTERNAL :: XXWReadLineInt
  DOUBLE PRECISION, EXTERNAL :: XXWReadLineDbl
  CHARACTER(LEN=200), EXTERNAL :: XXWReadLineStr
  LOGICAL, EXTERNAL :: XXWReadLineBln
  LOGICAL, EXTERNAL :: InBasin
  LOGICAL, EXTERNAL :: InBasin2

  INTEGER, EXTERNAL :: myDdiff
  CHARACTER(LEN=14),EXTERNAL :: myDtoStr
  INTEGER, EXTERNAL :: myDCompare,myDCompare2
  LOGICAL, EXTERNAL :: myDEqual
  CHARACTER(LEN=6), EXTERNAL:: PosIntMsg
  DOUBLE PRECISION, EXTERNAL :: QUnitFactor

  END MODULE CREST_Project

!============================================================
!
!============================================================
  MODULE CREST_Param
  IMPLICIT NONE

  TYPE CREST_Params
    DOUBLE PRECISION, ALLOCATABLE :: RainFact(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: Ksat(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: WM(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: B(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: IM(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KE(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: coeM(:,:)

    DOUBLE PRECISION, ALLOCATABLE :: expM(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: coeR(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: coeS(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KS(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: KI(:,:)

  END TYPE CREST_Params

  TYPE(CREST_Params):: g_tParams_Cali
  TYPE(CREST_Params):: g_tParams
  TYPE(CREST_Params):: g_tParamsAdj_Cali

  END MODULE CREST_Param

!============================================================
!
!============================================================
  MODULE CREST_ICS
  IMPLICIT NONE

  DOUBLE PRECISION, ALLOCATABLE :: g_WU0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: g_SS0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: g_SI0(:,:)

  END MODULE CREST_ICS


!============================================================
!
!============================================================
  MODULE CREST_Calib_SCEUA
  IMPLICIT NONE
  INTEGER :: g_iSeed
  INTEGER :: g_maxn
  INTEGER :: g_kstop
  DOUBLE PRECISION :: g_pcento
  INTEGER :: g_ngs
  INTEGER :: g_nopt

  INTEGER,ALLOCATABLE::g_RegMask(:,:)
  INTEGER :: g_RegNum
  INTEGER :: g_NCalibStas

  TYPE:: CREST_CalibStations
    CHARACTER(LEN=200):: Name
    INTEGER :: Value
    INTEGER :: Col
    INTEGER :: Row
    LOGICAL :: bIsOut

    INTEGER :: Area !Area or Number of the Cells in FAC
    INTEGER, ALLOCATABLE :: Mask(:,:)

    DOUBLE PRECISION, ALLOCATABLE :: RObs(:),RSim(:)

    DOUBLE PRECISION, ALLOCATABLE :: a(:)
    DOUBLE PRECISION, ALLOCATABLE :: bl(:)
    DOUBLE PRECISION, ALLOCATABLE :: bu(:)
    DOUBLE PRECISION, ALLOCATABLE :: x(:)

  END TYPE CREST_CalibStations

  TYPE(CREST_CalibStations), ALLOCATABLE :: g_tCalibSta(:)

  CHARACTER(LEN=10)::g_sParamName(1:12)
  DATA g_sParamName /"RainFact","Ksat","WM","B","IM", &
    "KE","coeM","expM","coeR","coeS", &
    "KS","KI"/

  END MODULE CREST_Calib_SCEUA

!============================================================
!
!============================================================
  MODULE CREST_Basic
  IMPLICIT NONE

  INTEGER :: g_NCols
  INTEGER :: g_NRows
  DOUBLE PRECISION :: g_xllCorner
  DOUBLE PRECISION :: g_yllCorner
  DOUBLE PRECISION :: g_CellSize
  DOUBLE PRECISION :: g_NoData_Value

  DOUBLE PRECISION, ALLOCATABLE :: g_DEM(:,:)
  INTEGER, ALLOCATABLE :: g_Mask(:,:) !Optional
  INTEGER, ALLOCATABLE :: g_FDR(:,:)
  INTEGER, ALLOCATABLE :: g_FAC(:,:)
  INTEGER, ALLOCATABLE :: g_Stream(:,:)

  DOUBLE PRECISION, ALLOCATABLE :: g_GridArea(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: g_NextLen(:,:)
  INTEGER, ALLOCATABLE :: g_NextR(:,:)
  INTEGER, ALLOCATABLE :: g_NextC(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: g_Slope(:,:)

  END MODULE CREST_Basic

!+ygan
!============================================================
!
!============================================================
  MODULE CREST_Diamond
      IMPLICIT NONE
      TYPE::swan_radarInfo
        CHARACTER(LEN=16):: RadarStationName(20)
        REAL*4 :: RadarLongitude(20)
        REAL*4 :: RadarLatitude(20)
        REAL*4 :: RadarAltitude(20)
        INTEGER*1 :: MosaicFlag(20)
      END TYPE swan_radarInfo

      TYPE::diamond131_header
        CHARACTER(LEN=12) :: ZonName
        CHARACTER(LEN=38) :: DataName
        CHARACTER(LEN=8) :: Flag
        CHARACTER(LEN=8) :: Version
        INTEGER*2 :: year
        INTEGER*2 :: month
        INTEGER*2 :: day
        INTEGER*2 :: hour
        INTEGER*2 :: minute
        INTEGER*2 :: interval
        INTEGER*2 :: XNumGrids
        INTEGER*2 :: YNumGrids
        INTEGER*2 :: ZNumGrids
        INTEGER*4 :: RadarCount
        REAL*4 :: StartLon
        REAL*4 :: StartLat
        REAL*4 :: CenterLon
        REAL*4 :: CenterLat
        REAL*4 :: XReso
        REAL*4 :: YReso
        REAL*4 :: ZhighGrids(40)
        TYPE(swan_radarInfo) :: radarsInfo
        INTEGER*1 Reserved(172)
      END TYPE diamond131_header

      INTEGER(2),PARAMETER::diam131_max_qpe=203*10
      INTEGER(2),PARAMETER::diam131_no_value=0

  END MODULE CREST_Diamond
!.ygan

!============================================================
!
!============================================================
  SUBROUTINE ExportGridVar(strDate,Rain,PET,EPot,EAct, &
    W,Runoff,ExcS,ExcI,RS,RI)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  CHARACTER*(*) strDate
  LOGICAL :: bIsError
  DOUBLE PRECISION :: Rain(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: PET(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EPot(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: W(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EAct(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: Runoff(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: dblTemp(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: iCountOfLakeMask,i
  DOUBLE PRECISION :: SumValueInLake

  !Output Grid Variables
  IF(g_bGOVar(0))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=Rain/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_"  &
      //TRIM(g_sGOVarName(0)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(1))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=PET/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(1)) // "_" //TRIM(strDate),   &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(2))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=EPot/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(2)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(3))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=EAct/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(3)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(4))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=W
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        ! SumValueInLake
        ! & =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        ! SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=g_tParams%WM*(1+g_tParams%B)
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_"  &
      //TRIM(g_sGOVarName(4)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(5))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=W*100.0/(g_tParams%WM*(1+g_tParams%B)) !Modified by Xianwu Xue 2011.4.20
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        ! SumValueInLake
        ! & =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        ! SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)
        SumValueInLake=100.0
        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_"  &
      //TRIM(g_sGOVarName(5)) // "_" //TRIM(strDate),   &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(6))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=Runoff
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_"   &
      //TRIM(g_sGOVarName(6)) // "_" //TRIM(strDate),   &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(7))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=ExcS/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(7)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(8))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=ExcI/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_"  &
      //TRIM(g_sGOVarName(8)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(9))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=RS/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(9)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF

  IF(g_bGOVar(10))THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=RI/g_TimeStep
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE

    ! Added by Xianwu Xue 2011.4.21
    IF(g_NLakes>0)THEN

      i=1
      DO WHILE(i<=g_NLakes)
        iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

        IF(iCountOfLakeMask==0)THEN
          i=i+1
          CYCLE
        END IF

        SumValueInLake &
          =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
        SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

        WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
          dblTemp=SumValueInLake
        END WHERE

        i=i+1
      END DO
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GOVar_" &
      //TRIM(g_sGOVarName(10)) // "_" //TRIM(strDate),  &
      dblTemp, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_ResultFormat)
  END IF
  RETURN
  END SUBROUTINE ExportGridVar


!============================================================
!
!============================================================
  SUBROUTINE ExportOutDTData(strDate,Rain,PET,EPot,EAct, &
    W,Runoff,ExcS,ExcI,RS,RI)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  CHARACTER*(*) strDate
  LOGICAL :: bIsError
  DOUBLE PRECISION :: Rain(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: PET(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EPot(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: W(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EAct(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: Runoff(0:g_NCols-1,0:g_NRows-1)


  DOUBLE PRECISION :: dblTemp(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: iCountOfLakeMask,i
  DOUBLE PRECISION :: SumValueInLake

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=Rain/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(0)) // "_" //TRIM(strDate),    &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=PET/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE

  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF

  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(1)) // "_" //TRIM(strDate),   &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=EPot/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE

  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF

  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(2)) // "_" //TRIM(strDate),  &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=EAct/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(3)) // "_" //TRIM(strDate),   &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=W
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      ! SumValueInLake
      ! & =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      ! SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=g_tParams%WM*(1+g_tParams%B)
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(4)) // "_" //TRIM(strDate),  &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)


  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=W*100.0/(g_tParams%WM*(1+g_tParams%B)) !Modified by Xianwu Xue 2011.4.20
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      ! SumValueInLake
      ! & =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      ! SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)
      SumValueInLake=100
      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_" &
    //TRIM(g_sGOVarName(5)) // "_" //TRIM(strDate),    &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=Runoff
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(6)) // "_" //TRIM(strDate),   &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=ExcS/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(7)) // "_" //TRIM(strDate),   &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)


  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=ExcI/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_" &
    //TRIM(g_sGOVarName(8)) // "_" //TRIM(strDate),  &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)


  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=RS/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_"  &
    //TRIM(g_sGOVarName(9)) // "_" //TRIM(strDate),   &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  WHERE(g_Mask/=g_NoData_Value)
    dblTemp=RI/g_TimeStep
  ELSEWHERE
    dblTemp=g_NoData_Value
  END WHERE
  ! Added by Xianwu Xue 2011.4.21
  IF(g_NLakes>0)THEN

    i=1
    DO WHILE(i<=g_NLakes)
      iCountOfLakeMask=COUNT((g_LakeMask==i) .AND.(dblTemp>=0))

      IF(iCountOfLakeMask==0)THEN
        i=i+1
        CYCLE
      END IF

      SumValueInLake &
        =SUM(dblTemp,(g_LakeMask==i) .AND.(dblTemp>=0))
      SumValueInLake=SumValueInLake/DBLE(iCountOfLakeMask)

      WHERE((g_LakeMask==i) .AND.(dblTemp>=0))
        dblTemp=SumValueInLake
      END WHERE

      i=i+1
    END DO
  END IF
  CALL WriteMatrixFile(TRIM(g_ResultPath)//"OutDT_" &
    //TRIM(g_sGOVarName(10)) // "_" //TRIM(strDate),  &
    dblTemp, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
    bIsError,g_ResultFormat)

  RETURN
  END SUBROUTINE ExportOutDTData

!============================================================
!
!============================================================
  SUBROUTINE CalculateOutletData(k,Rain,PET,EPot,EAct, &
    W,Runoff,ExcS,ExcI,RS,RI)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  INTEGER :: N,k
  DOUBLE PRECISION :: Rain(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: PET(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EPot(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: W(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EAct(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: Runoff(0:g_NCols-1,0:g_NRows-1)
  !-----------------------------------------------
  !Output the state of outlet
  N=COUNT(g_tOutlet%Mask/= g_NoData_Value)
  IF(g_tOutlet%bIsOut .EQV. .FALSE.)THEN
    g_tOutlet%Rain(k)=SUM(Rain/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%PET(k)=SUM(PET/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%EPot(k)=SUM(EPot/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%EAct(k)=SUM(EAct/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%W(k)=SUM(W, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%SM(k)=SUM(W/g_tParams%WM, &
      g_tOutlet%Mask/= g_NoData_Value)/N

    g_tOutlet%R(k)=Runoff(g_tOutlet%Col, &
      g_tOutlet%Row)

    g_tOutlet%ExcS(k)=SUM(ExcS/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%ExcI(k)=SUM(ExcI/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%RS(k)=SUM(RS/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
    g_tOutlet%RI(k)=SUM(RI/g_TimeStep, &
      g_tOutlet%Mask/= g_NoData_Value)/N
  ELSE
    g_tOutlet%Rain(k)=g_NoData_Value
    g_tOutlet%PET(k)=g_NoData_Value
    g_tOutlet%EPot(k)=g_NoData_Value
    g_tOutlet%EAct(k)=g_NoData_Value
    g_tOutlet%W(k)=g_NoData_Value
    g_tOutlet%SM(k)=g_NoData_Value
    g_tOutlet%R(k)=g_NoData_Value
    g_tOutlet%ExcS(k)=g_NoData_Value
    g_tOutlet%ExcI(k)=g_NoData_Value
    g_tOutlet%RS(k)=g_NoData_Value
    g_tOutlet%RI(k)=g_NoData_Value
  END IF
  !print*,g_tOutlet%Col, g_tOutlet%Row, g_tOutlet%R(k)
  RETURN
  END SUBROUTINE CalculateOutletData

!============================================================
!
!============================================================
  SUBROUTINE ExportOutletData()
  USE CREST_Project
  IMPLICIT NONE
  INTEGER :: fileid
  INTEGER :: i,k
  INTEGER :: dtTemp(1:6)
  DOUBLE PRECISION :: NSCE, Bias, CC

  CALL XXWGetFreeFile(fileid)


  IF (g_tOutlet%bIsOut .EQV. .TRUE.) THEN !Outlet is not in Basin area
    RETURN
  END IF

  OPEN(fileid,FILE=TRIM(g_ResultPath)//"Outlet_" &
    //TRIM(g_tOutlet%Name) // "_Results.csv", &
    FORM="FORMATTED")
  IF(g_tOutlet%HasOBS .EQV. .TRUE.)THEN
    WRITE(fileid,*) "DateTime,", &
      "Rain,", &
      "PET,", &
      "EPot,", &
      "EAct,", &
      "W,", &
      "SM,", &
      "RS,", &
      "RI,", &
      "ExcS,", &
      "ExcI,", &
      "R,","RObs"
  ELSE
    WRITE(fileid,*) "DateTime,", &
      "Rain,", &
      "PET,", &
      "EPot,", &
      "EAct,", &
      "W,", &
      "SM,", &
      "RS,", &
      "RI,", &
      "ExcS,", &
      "ExcI,", &
      "R"
  END IF

  dtTemp=g_StartDate

  DO k=0, UBOUND(g_tOutlet%Rain,1)
    WRITE(fileid,"(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)", &
      ADVANCE='no') &
      dtTemp(1),"-", &
      dtTemp(2),"-", &
      dtTemp(3)," ", &
      dtTemp(4),":", &
      dtTemp(5),":", &
      dtTemp(6),","

    WRITE(fileid,"(10(F10.3,A1),F10.3)",ADVANCE='no') &
      g_tOutlet%Rain(k),",", &
      g_tOutlet%PET(k),",", &
      g_tOutlet%EPot(k),",", &
      g_tOutlet%EAct(k),",", &
      g_tOutlet%W(k),",", &
      g_tOutlet%SM(k),",", &
      g_tOutlet%RS(k),",", &
      g_tOutlet%RI(k),",", &
      g_tOutlet%ExcS(k),",", &
      g_tOutlet%ExcI(k),",", &
      g_tOutlet%R(k)
    IF(g_tOutlet%HasOBS .EQV. .TRUE.)THEN
      WRITE(fileid,"(A1,F10.3)") &
        ",",g_tOutlet%RObs(k)
    ELSE
      WRITE(fileid,*)
    END IF
    CALL myDAdd(g_TimeMark,g_TimeStep,dtTemp)
  END DO
  CLOSE(fileid)

  IF(g_tOutlet%HasOBS .EQV. .TRUE.)THEN
    CALL XXWGetFreeFile(fileid)
    OPEN(fileid,FILE=TRIM(g_ResultPath)//"Outlet_" &
      //TRIM(g_tOutlet%Name) // "_Results_Statistics.csv", &
      FORM="FORMATTED")

    CALL GetNSCE(g_tOutlet%RObs, &
      g_tOutlet%R,g_NumWarmup,g_ITMax, &
      NSCE)
    CALL GetBias(g_tOutlet%RObs, &
      g_tOutlet%R,g_NumWarmup,g_ITMax, &
      Bias)
    CALL GetCC(g_tOutlet%RObs, &
      g_tOutlet%R,g_NumWarmup,g_ITMax, &
      CC)
    WRITE(fileid,*)"NSCE,",NSCE
    WRITE(fileid,*)"Bias(%),",Bias
    WRITE(fileid,*)"CC,",CC

    CLOSE(fileid)
  END IF
  RETURN
  END SUBROUTINE ExportOutletData

!============================================================
!
!============================================================
  SUBROUTINE CalculateOutPixData(k,Rain,PET,EPot,EAct, &
    W,Runoff,ExcS,ExcI,RS,RI)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  INTEGER :: N,k,i
  DOUBLE PRECISION :: Rain(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: PET(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EPot(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: W(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: ExcI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RI(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: RS(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: EAct(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: Runoff(0:g_NCols-1,0:g_NRows-1)
  !-----------------------------------------------
  !Output the state of the Pixes

  DO i = 0, g_NOutPixs-1
    N=COUNT(g_tOutPix(i)%Mask/= g_NoData_Value)
    IF (g_tOutPix(i)%bIsOut .EQV. .FALSE.) THEN
      g_tOutPix(i)%Rain(k)=SUM(Rain/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%PET(k)=SUM(PET/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%EPot(k)=SUM(EPot/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%EAct(k)=SUM(EAct/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%W(k)=SUM(W, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%SM(k)=SUM(W/g_tParams%WM, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N

      g_tOutPix(i)%R(k)=Runoff(g_tOutPix(i)%Col, &
        g_tOutPix(i)%Row)

      g_tOutPix(i)%ExcS(k)=SUM(ExcS/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%ExcI(k)=SUM(ExcI/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%RS(k)=SUM(RS/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
      g_tOutPix(i)%RI(k)=SUM(RI/g_TimeStep, &
        g_tOutPix(i)%Mask/= g_NoData_Value)/N
    ELSE
      g_tOutPix(i)%Rain(k)=g_NoData_Value
      g_tOutPix(i)%PET(k)=g_NoData_Value
      g_tOutPix(i)%EPot(k)=g_NoData_Value
      g_tOutPix(i)%EAct(k)=g_NoData_Value
      g_tOutPix(i)%W(k)=g_NoData_Value
      g_tOutPix(i)%SM(k)=g_NoData_Value
      g_tOutPix(i)%R(k)=g_NoData_Value
      g_tOutPix(i)%ExcS(k)=g_NoData_Value
      g_tOutPix(i)%ExcI(k)=g_NoData_Value
      g_tOutPix(i)%RS(k)=g_NoData_Value
      g_tOutPix(i)%RI(k)=g_NoData_Value
    END IF
  END DO
  RETURN
  END SUBROUTINE CalculateOutPixData

!============================================================
!
!============================================================
  SUBROUTINE ExportPixData()
  USE CREST_Project
  IMPLICIT NONE
  INTEGER :: fileid
  INTEGER :: i,k
  INTEGER :: dtTemp(1:6)
  DOUBLE PRECISION :: NSCE, Bias, CC

  DO i = 0, g_NOutPixs-1
    IF (g_tOutPix(i)%bIsOut .EQV. .TRUE.) THEN !OutPix is not in Basin area
      CYCLE
    END IF

    CALL XXWGetFreeFile(fileid)

    OPEN(fileid,FILE=TRIM(g_ResultPath)//"OutPix_"  &
      //TRIM(g_tOutPix(i)%Name) // "_Results.csv",  &
      FORM="FORMATTED")
    IF(g_tOutPix(i)%HasOBS .EQV. .TRUE.)THEN
      WRITE(fileid,*) "DateTime,", &
        "Rain,", &
        "PET,", &
        "EPot,", &
        "EAct,", &
        "W,", &
        "SM,", &
        "RS,", &
        "RI,", &
        "ExcS,", &
        "ExcI,", &
        "R,","RObs"
    ELSE
      WRITE(fileid,*) "DateTime,", &
        "Rain,", &
        "PET,", &
        "EPot,", &
        "EAct,", &
        "W,", &
        "SM,", &
        "RS,", &
        "RI,", &
        "ExcS,", &
        "ExcI,", &
        "R"
    END IF

    dtTemp=g_StartDate

    DO k=0, UBOUND(g_tOutPix(i)%Rain,1)
      WRITE(fileid,"(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1)", &
        ADVANCE='no') &
        dtTemp(1),"-", &
        dtTemp(2),"-", &
        dtTemp(3)," ", &
        dtTemp(4),":", &
        dtTemp(5),":", &
        dtTemp(6),","

      WRITE(fileid,"(10(F10.3,A1),F10.3)",ADVANCE='no') &
        g_tOutPix(i)%Rain(k),",", &
        g_tOutPix(i)%PET(k),",", &
        g_tOutPix(i)%EPot(k),",", &
        g_tOutPix(i)%EAct(k),",", &
        g_tOutPix(i)%W(k),",", &
        g_tOutPix(i)%SM(k),",", &
        g_tOutPix(i)%RS(k),",", &
        g_tOutPix(i)%RI(k),",", &
        g_tOutPix(i)%ExcS(k),",", &
        g_tOutPix(i)%ExcI(k),",", &
        g_tOutPix(i)%R(k)

      IF(g_tOutPix(i)%HasOBS .EQV. .TRUE.)THEN
        WRITE(fileid,"(A1,F10.3)") &
          ",",g_tOutPix(i)%RObs(k)
      ELSE
        WRITE(fileid,*)
      END IF


      CALL myDAdd(g_TimeMark,g_TimeStep,dtTemp)
    END DO

    CLOSE(fileid)

    IF(g_tOutPix(i)%HasOBS .EQV. .TRUE.)THEN
      CALL XXWGetFreeFile(fileid)
      OPEN(fileid,FILE=TRIM(g_ResultPath)//"OutPix_"  &
        //TRIM(g_tOutPix(i)%Name) // "_Results_Statistics.csv", &
        FORM="FORMATTED")

      CALL GetNSCE(g_tOutPix(i)%RObs, &
        g_tOutPix(i)%R,g_NumWarmup,g_ITMax, &
        NSCE)
      CALL GetBias(g_tOutPix(i)%RObs, &
        g_tOutPix(i)%R,g_NumWarmup,g_ITMax, &
        Bias)
      CALL GetCC(g_tOutPix(i)%RObs, &
        g_tOutPix(i)%R,g_NumWarmup,g_ITMax, &
        CC)
      WRITE(fileid,*)"NSCE,",NSCE
      WRITE(fileid,*)"Bias(%),",Bias
      WRITE(fileid,*)"CC,",CC

      CLOSE(fileid)
    END IF
  END DO

  RETURN
  END SUBROUTINE ExportPixData

!============================================================
!
!============================================================
  FUNCTION PosIntMsg(I)
  IMPLICIT NONE
  InteGEr ::I
  Character(LEN=6):: PosIntMsg
  IF (I<10) THEN
    WRITE(PosIntMsg, "(I1)") I
  ELSEIF (I<100) THEN
    WRITE(PosIntMsg, "(I2)") I
  ELSEIF (I<1000) THEN
    WRITE(PosIntMsg, "(I3)") I
  ELSEIF (I<10000) THEN
    WRITE(PosIntMsg, "(I4)") I
  ELSEIF (I<100000) THEN
    WRITE(PosIntMsg, "(I5)") I
  ELSE
    WRITE(PosIntMsg, "(I6)") I
  END IF
  END FUNCTION PosIntMsg

!============================================================
!
!============================================================
  SUBROUTINE myStrtoD(Str, myD, Mark)
  INTEGER(4):: myD(1:6), I
  CHARACTER(14):: Str, msg
  CHARACTER(1)::Mark

  myD=0

  SELECT CASE (Mark)
  CASE ("y","Y")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    myD(2)=1
    myD(3)=1
    Str(5:14)="0101000000"
  CASE ("m","M")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    msg=Str(5:6)
    READ(msg,*) myD(2)
    myD(3)=1
    Str(7:14)="01000000"
  CASE ("d","D")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    msg=Str(5:6)
    READ(msg,*) myD(2)
    msg=Str(7:8)
    READ(msg,*) myD(3)
    Str(9:14)="000000"
  CASE ("h","H")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    msg=Str(5:6)
    READ(msg,*) myD(2)
    msg=Str(7:8)
    READ(msg,*) myD(3)
    msg=Str(9:10)
    READ(msg,*) myD(4)
    Str(11:14)="0000"
  CASE ("u","U")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    msg=Str(5:6)
    READ(msg,*) myD(2)
    msg=Str(7:8)
    READ(msg,*) myD(3)
    msg=Str(9:10)
    READ(msg,*) myD(4)
    msg=Str(11:12)
    READ(msg,*) myD(5)
    Str(13:14)="00"
  CASE ("s","S")
    msg=Str(1:4)
    READ(msg,*) myD(1)
    msg=Str(5:6)
    READ(msg,*) myD(2)
    msg=Str(7:8)
    READ(msg,*) myD(3)
    msg=Str(9:10)
    READ(msg,*) myD(4)
    msg=Str(11:12)
    READ(msg,*) myD(5)
    msg=Str(13:14)
    READ(msg,*) myD(6)
  END SELECT

  END SUBROUTINE myStrtoD

!============================================================
!
!============================================================
  INTEGER FUNCTION myDdiff(Mark, myDAout, myDBout)
  IMPLICIT NONE
  CHARACTER(1):: Mark
  INTEGER(4):: myDAout(1:6), myDBout(1:6)
  INTEGER(4):: myDA(1:6), myDB(1:6), DX(1:6), Num, Sign
  INTEGER(4), EXTERNAL:: myDcompare

  myDA=myDAout
  myDB=myDBout
  Sign=myDcompare(myDA, myDB,Mark)
  IF (Sign==-1) THEN
    myDA=myDBout
    myDB=myDAout
  END IF

  SELECT CASE (Mark)
  CASE ("h","H", "u", "U", "s", "S")
    DX=myDA
    Num=0
    DO
      IF (DX(1)==myDB(1) .AND. DX(2)==myDB(2) &
        .AND. DX(3)==myDB(3)) EXIT
      CALL myDadd("d", 1, DX)
      Num=Num+1
    END DO

    SELECT CASE (Mark)
    CASE ("h","H")
      myDdiff=Num*24+myDB(4)-myDA(4)
    CASE ("u","U")
      myDdiff=(Num*24+myDB(4)-myDA(4))*60+myDB(5)-myDA(5)
    CASE ("s","S")
      myDdiff=((Num*24+myDB(4)-myDA(4))*60 &
        +myDB(5)-myDA(5))*60+myDB(6)-myDA(6)
    END SELECT
    CASE DEFAULT
    DX=myDA
    Num=0
    DO
      SELECT CASE (Mark)
      CASE ("y", "Y")
        IF (DX(1)==myDB(1)) EXIT
      CASE ("m","M")
        IF (DX(1)==myDB(1) .AND. DX(2)==myDB(2)) EXIT
      CASE ("d", "D")
        IF (DX(1)==myDB(1) .AND. DX(2)==myDB(2) &
          .AND. DX(3)==myDB(3)) EXIT
      END SELECT
      CALL myDadd(Mark, 1, DX)
      Num=Num+1
    END DO
    myDdiff=Num
  END SELECT
  myDdiff=myDdiff*Sign
  END FUNCTION myDdiff

!============================================================
!
!============================================================
  SUBROUTINE myDAdd(Mark, Num, myD)
  IMPLICIT NONE
  CHARACTER(1):: Mark
  INTEGER(4):: myD(1:6)
  INTEGER(4):: Num, I, Days

  SELECT CASE (Mark)
  CASE ("y","Y"); I=1
  CASE ("m","M"); I=2
  CASE ("d","D"); I=3
  CASE ("h","H"); I=4
  CASE ("u","U"); I=5
  CASE ("s","S"); I=6
  END SELECT
  myD(I)=myD(I)+Num

  IF (I==6) THEN ! second
    IF (myD(6)<0) THEN
      DO WHILE (myD(6)<0)
        myD(6)=myD(6)+60
        myD(5)=myD(5)-1
      END DO
    ELSE
      DO WHILE (myD(6)>=60)
        myD(6)=myD(6)-60
        myD(5)=myD(5)+1
      END DO
    END IF
  END IF
  IF (I>=5) THEN ! minute
    IF (myD(5)<0) THEN
      DO WHILE (myD(5)<0)
        myD(5)=myD(5)+60
        myD(4)=myD(4)-1
      END DO
    ELSE
      DO WHILE (myD(5)>=60)
        myD(5)=myD(5)-60
        myD(4)=myD(4)+1
      END DO
    END IF
  END IF
  IF (I>=4) THEN ! hour
    IF (myD(4)<0) THEN
      DO WHILE (myD(4)<0)
        myD(4)=myD(4)+24
        myD(3)=myD(3)-1
      END DO
    ELSE
      DO WHILE (myD(4)>=24)
        myD(4)=myD(4)-24
        myD(3)=myD(3)+1
      END DO
    END IF
  END IF
  IF (I>=3) THEN ! day
    IF (myD(3)<0) THEN
      DO WHILE (myD(3)<0)
        SELECT CASE (myD(2))
        CASE (1,3,5,7,8,10,12)
          Days=31
        CASE (4,6,9,11)
          Days=30
        CASE (2)
          IF (Mod(myD(1), 4)==0) THEN
            Days=29
          ELSE
            Days=28
          END IF
        END SELECT
        myD(3)=myD(3)+days
        myD(2)=myD(2)-1
      END DO
    ELSE IF(myD(3)==0) THEN
      myD(2)=myD(2)-1
      IF (myD(2)==0) THEN
        myD(2)=12
        myD(1)=myD(1)-1
      END IF
      SELECT CASE (myD(2))
      CASE (1,3,5,7,8,10,12)
        Days=31
      CASE (4,6,9,11)
        Days=30
      CASE (2)
        IF (Mod(myD(1), 4)==0) THEN
          Days=29
        ELSE
          Days=28
        END IF
      END SELECT
      myD(3)=Days
    ELSE
      DO
        SELECT CASE (myD(2))
        CASE (1,3,5,7,8,10,12)
          Days=31
        CASE (4,6,9,11)
          Days=30
        CASE (2)
          IF (Mod(myD(1), 4)==0) THEN
            Days=29
          ELSE
            Days=28
          END IF
        END SELECT
        IF (myD(3)<=Days) THEN
          EXIT
        ELSE
          myD(3)=myD(3)-Days
          myD(2)=myD(2)+1
        END IF
      END DO
    END IF
  END IF
  IF (I>=2) THEN ! Month
    IF (myD(2)<=0) THEN
      DO WHILE (myD(2)<0)
        myD(2)=myD(2)+12
        myD(1)=myD(1)-1
      END DO
    ELSEIF (myD(2)==0) THEN
      myD(2)=12
      myD(1)=myD(1)-1
    ELSE
      DO WHILE (myD(2)>12)
        myD(1)=myD(1)+1
        myD(2)=myD(2)-12
      END DO
    END IF
  END IF
  END SUBROUTINE myDAdd


!============================================================
!
!============================================================
  FUNCTION myDcompare(DA, DB,Mark)
  IMPLICIT NONE
  INTEGER(4):: DA(1:6), DB(1:6), myDcompare,I
  CHARACTER(1):: Mark

  SELECT CASE (Mark)
  CASE ("y","Y"); I=1
  CASE ("m","M"); I=2
  CASE ("d","D"); I=3
  CASE ("h","H"); I=4
  CASE ("u","U"); I=5
  CASE ("s","S"); I=6
  END SELECT

  myDcompare=1
  IF (DA(1)>DB(1) .AND. I>=1) THEN
    myDcompare=-1
  ElseIF (DA(1)==DB(1) .AND. DA(2)>DB(2) .AND. I>=2) THEN
    myDcompare=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)>DB(3) .AND. I>=3) THEN
  myDcompare=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)>DB(4) .AND. I>=4) THEN
  myDcompare=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)==DB(4) &
    .AND. DA(5)>DB(5) .AND. I>=5) THEN
  myDcompare=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)==DB(4) &
    .AND. DA(5)==DB(5) &
    .AND. DA(6)>DB(6) .AND. I>=6) THEN
  myDcompare=-1
  END IF
  END FUNCTION myDcompare

!============================================================
!
!============================================================
  FUNCTION myDcompare2(DA, DB,Mark)
  IMPLICIT NONE
  INTEGER(4):: DA(1:6), DB(1:6), myDcompare2,I
  CHARACTER(1):: Mark

  SELECT CASE (Mark)
  CASE ("y","Y"); I=1
  CASE ("m","M"); I=2
  CASE ("d","D"); I=3
  CASE ("h","H"); I=4
  CASE ("u","U"); I=5
  CASE ("s","S"); I=6
  END SELECT

  myDcompare2=1
  IF (DA(1)>=DB(1) .AND. I>=1) THEN
    myDcompare2=-1
  ElseIF (DA(1)==DB(1) .AND. DA(2)>=DB(2) .AND. I>=2) THEN
    myDcompare2=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)>=DB(3) .AND. I>=3) THEN
  myDcompare2=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)>=DB(4) .AND. I>=4) THEN
  myDcompare2=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)==DB(4) &
    .AND. DA(5)>=DB(5) .AND. I>=5) THEN
  myDcompare2=-1
  ELSEIF (DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)==DB(4) &
    .AND. DA(5)==DB(5) &
    .AND. DA(6)>=DB(6) .AND. I>=6) THEN
  myDcompare2=-1
  END IF
  END FUNCTION myDcompare2

!============================================================
!
!============================================================
  FUNCTION myDtoStr(myD, Mark)
  IMPLICIT NONE
  INTEGER(4):: myD(1:6)
  CHARACTER(1):: Mark
  CHARACTER(14):: myDtoStr, msg

  myDtoStr=""
  WRITE(msg,"(I4)") myD(1)
  myDtoStr(1:4)=msg(1:4)

  SELECT CASE (Mark)
  CASE ("m","M")
    WRITE(msg,"(I2)") myD(2)
    myDtoStr(5:6)=msg(1:2)
    IF (myD(2)<10) myDtoStr(5:5)="0"
  CASE ("d","D")
    WRITE(msg,"(I2)") myD(2)
    myDtoStr(5:6)=msg(1:2)
    IF (myD(2)<10) myDtoStr(5:5)="0"
    WRITE(msg,"(I2)") myD(3)
    myDtoStr(7:8)=msg(1:2)
    IF (myD(3)<10) myDtoStr(7:7)="0"
  CASE ("h","H")
    WRITE(msg,"(I2)") myD(2)
    myDtoStr(5:6)=msg(1:2)
    IF (myD(2)<10) myDtoStr(5:5)="0"
    WRITE(msg,"(I2)") myD(3)
    myDtoStr(7:8)=msg(1:2)
    IF (myD(3)<10) myDtoStr(7:7)="0"
    WRITE(msg,"(I2)") myD(4); myDtoStr(9:10)=msg(1:2)
    IF (myD(4)<10) myDtoStr(9:9)="0"
  CASE ("u","U")
    WRITE(msg,"(I2)") myD(2)
    myDtoStr(5:6)=msg(1:2)
    IF (myD(2)<10) myDtoStr(5:5)="0"
    WRITE(msg,"(I2)") myD(3)
    myDtoStr(7:8)=msg(1:2)
    IF (myD(3)<10) myDtoStr(7:7)="0"
    WRITE(msg,"(I2)") myD(4); myDtoStr(9:10)=msg(1:2)
    IF (myD(4)<10) myDtoStr(9:9)="0"
    WRITE(msg,"(I2)") myD(5); myDtoStr(11:12)=msg(1:2)
    IF (myD(5)<10) myDtoStr(11:11)="0"
  CASE ("s","S")
    WRITE(msg,"(I2)") myD(2)
    myDtoStr(5:6)=msg(1:2)
    IF (myD(2)<10) myDtoStr(5:5)="0"
    WRITE(msg,"(I2)") myD(3)
    myDtoStr(7:8)=msg(1:2)
    IF (myD(3)<10) myDtoStr(7:7)="0"
    WRITE(msg,"(I2)") myD(4); myDtoStr(9:10)=msg(1:2)
    IF (myD(4)<10) myDtoStr(9:9)="0"
    WRITE(msg,"(I2)") myD(5); myDtoStr(11:12)=msg(1:2)
    IF (myD(5)<10) myDtoStr(11:11)="0"
    WRITE(msg,"(I2)") myD(6); myDtoStr(13:14)=msg(1:2)
    IF (myD(6)<10) myDtoStr(13:13)="0"
  END SELECT
  END FUNCTION myDtoStr

!============================================================
!
!============================================================
  SUBROUTINE myNow(DX)
  IMPLICIT NONE
  INTEGER(4):: DX(1:6)
  CHARACTER(20):: A, B
  CHARACTER(14):: strD

  CALL DATE_AND_TIME(A, B)
  strD(1:8)=A(1:8)
  strD(9:14)=B(1:6)
  CALL myStrtoD(strD, DX, "s")
  END SUBROUTINE myNow

!============================================================
!
!============================================================
  FUNCTION myDEqual(DA, DB)
  INTEGER(4):: DA(1:6), DB(1:6)
  LOGICAL:: myDEqual
  myDEqual=(DA(1)==DB(1) &
    .AND. DA(2)==DB(2) &
    .AND. DA(3)==DB(3) &
    .AND. DA(4)==DB(4) &
    .AND. DA(5)==DB(5) &
    .AND. DA(6)==DB(6))
  END FUNCTION myDEqual


!============================================================
!
!============================================================
  SUBROUTINE UtiGetAndDisplayBeginTime(IBDT)

  INTEGER, INTent(out):: IBDT(8)
  INTEGER :: i

  CALL DATE_AND_TIME(VALUES=IBDT)

  !IF(IBDT(1)>=2015)THEN
  ! WRITE(*,*)"Your license expired, please contract xuexianwu@ou.edu!"
  ! WRITE(*,*)
  ! STOP
  !END IF

  WRITE(*,*)
  WRITE(*,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
2 FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ', &
    I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
  WRITE(*,*)

  END SUBROUTINE UtiGetAndDisplayBeginTime

!============================================================
!
!============================================================
  SUBROUTINE UtiGetAndDisplayBeginTimeToFile(IBDT)
  USE CREST_Project

  INTEGER, INTent(in):: IBDT(8)
  INTEGER :: i


2 FORMAT(1X,'Run start date and time (yyyy/mm/dd hh:mm:ss): ', &
    I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)

  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,2) (IBDT(I),I=1,3),(IBDT(I),I=5,7)
  WRITE(g_CREST_LogFileID,*)

  END SUBROUTINE UtiGetAndDisplayBeginTimeToFile

!============================================================
!
!============================================================
  SUBROUTINE UtiGetAndDisplayEndTime(IOUT,IBDT,IPRTIM)
  USE CREST_Project

  INTEGER IBDT(8), IEDT(8), IDPM(12)
  DATA IDPM/31,28,31,30,31,30,31,31,30,31,30,31/
  DATA NSPD/86400/

  CALL DATE_AND_TIME(VALUES=IEDT)

  WRITE(*,*)
  WRITE(*,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)

  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)

1000 FORMAT(1X,'Run end date and time (yyyy/mm/dd hh:mm:ss): ', &
    I4,'/',I2.2,'/',I2.2,1X,I2,':',I2.2,':',I2.2)
  IF(IPRTIM.GT.0) THEN
    WRITE(IOUT,'(1X)')
    WRITE(IOUT,1000) (IEDT(I),I=1,3),(IEDT(I),I=5,7)
  END IF
  !
  ! Calculate elapsed time in days and seconds
  NDAYS=0
  LEAP=0
  IF (MOD(IEDT(1),4).EQ.0) LEAP = 1
  IBD = IBDT(3) ! BEGIN DAY
  IED = IEDT(3) ! END DAY
  ! FIND DAYS
  IF (IBDT(2).NE.IEDT(2)) THEN
    ! MONTHS DIFFER
    MB = IBDT(2) ! BEGIN MONTH
    ME = IEDT(2) ! END MONTH
    NM = ME-MB+1 ! NUMBER OF MONTHS TO LOOK AT
    IF (MB.GT.ME) NM = NM+12
    MC=MB-1
    DO 10 M=1,NM
      MC=MC+1 ! MC IS CURRENT MONTH
      IF (MC.EQ.13) MC = 1
      IF (MC.EQ.MB) THEN
        NDAYS = NDAYS+IDPM(MC)-IBD
        IF (MC.EQ.2) NDAYS = NDAYS + LEAP
      ELSEIF (MC.EQ.ME) THEN
        NDAYS = NDAYS+IED
      ELSE
        NDAYS = NDAYS+IDPM(MC)
        IF (MC.EQ.2) NDAYS = NDAYS + LEAP
      ENDIF
10 CONTINUE
  ELSEIF (IBD.LT.IED) THEN
    ! START AND END IN SAME MONTH, ONLY ACCOUNT FOR DAYS
    NDAYS = IED-IBD
  ENDIF
  ELSEC=NDAYS*NSPD
  !
  ! ADD OR SUBTRACT SECONDS
  ELSEC = ELSEC+(IEDT(5)-IBDT(5))*3600.0
  ELSEC = ELSEC+(IEDT(6)-IBDT(6))*60.0
  ELSEC = ELSEC+(IEDT(7)-IBDT(7))
  ELSEC = ELSEC+(IEDT(8)-IBDT(8))*0.001
  !
  ! CONVERT SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
  NDAYS = ELSEC/NSPD
  RSECS = MOD(ELSEC,86400.0)
  NHOURS = RSECS/3600.0
  RSECS = MOD(RSECS,3600.0)
  NMINS = RSECS/60.0
  RSECS = MOD(RSECS,60.0)
  NSECS = RSECS
  RSECS = MOD(RSECS,1.0)
  MSECS = NINT(RSECS*1000.0)
  NRSECS = NSECS
  IF (RSECS.GE.0.5) NRSECS=NRSECS+1
  !
  ! WRITE elapsed time to screen
  IF (NDAYS.GT.0) THEN
    WRITE(*,1010) NDAYS,NHOURS,NMINS,NRSECS

    WRITE(g_CREST_LogFileID,1010) NDAYS,NHOURS,NMINS,NRSECS

1010 FORMAT(1X,'Elapsed run time: ',I3,' Days, ',I2,' Hours, ',I2, &
      ' Minutes, ',I2,' Seconds',/)
  ELSEIF (NHOURS.GT.0) THEN
    WRITE(*,1020) NHOURS,NMINS,NRSECS

    WRITE(g_CREST_LogFileID,1020) NHOURS,NMINS,NRSECS

1020 FORMAT(1X,'Elapsed run time: ',I2,' Hours, ',I2, &
      ' Minutes, ',I2,' Seconds',/)
  ELSEIF (NMINS.GT.0) THEN
    WRITE(*,1030) NMINS,NSECS,MSECS

    WRITE(g_CREST_LogFileID,1030) NMINS,NSECS,MSECS

1030 FORMAT(1X,'Elapsed run time: ',I2,' Minutes, ', &
      I2,'.',I3.3,' Seconds',/)
  ELSE
    WRITE(*,1040) NSECS,MSECS

    WRITE(g_CREST_LogFileID,1040) NSECS,MSECS

1040 FORMAT(1X,'Elapsed run time: ',I2,'.',I3.3,' Seconds',/)
  ENDIF
  !
  RETURN
  END SUBROUTINE UtiGetAndDisplayEndTime



!============================================================
! Main Program of CREST
!============================================================
  PROGRAM CREST_Main

  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  USE CREST_ICS

  IMPLICIT NONE

  INTEGER :: IBDT(8),IEDT(8)
  INTEGER :: I
  LOGICAL :: bIsError

  CALL GetProjectFNameAndPath(g_ProjectName, &
    g_ProjectWS,g_PrjNP)

  g_PrjNP=TRIM(g_ProjectWS) // TRIM(g_ProjectName) // ".Project"

  ! Write to the Screen
  WRITE(*,*)
  WRITE(*,*)
  WRITE(*,"(42X,A)") "CREST"
  WRITE(*,"(21X,A)") "COUPLED ROUTING " &
    // "AND EXCESS STORAGE (OU & NASA)"
  WRITE(*,"(31X,'Version ',A)") g_CREST_Version

  CALL UtiGetAndDisplayBeginTime(IBDT)
  WRITE(g_strD,"(I4,'.',I2.2,'.',I2.2,'-',I2.2,'.',I2.2,'.',I2.2)") &
    (IBDT(I),I=1,3),(IBDT(I),I=5,7)

  !Write to the log file

  CALL XXWGetFreeFile(g_CREST_LogFileID)
  OPEN(g_CREST_LogFileID,FILE &
    =TRIM(g_ProjectWS)//TRIM(g_ProjectName)//"_" &
    //TRIM(g_strD)//"_CREST.log",FORM='FORMATTED')
  WRITE(g_CREST_LogFileID,"(42X,A)") "CREST"
  WRITE(g_CREST_LogFileID,"(21X,A)") "COUPLED ROUTING " &
    // "AND EXCESS STORAGE (OU & NASA)"
  WRITE(g_CREST_LogFileID,"(31X,'Version ',A)") g_CREST_Version

  !Get the begin time and then display on Screen
  CALL UtiGetAndDisplayBeginTimeToFile(IBDT)

  CALL CREST_Main_Pre(bIsError)
  IF(bIsError .EQV. .TRUE.)THEN
    STOP
  END IF

  WRITE(*,*)
  WRITE(g_CREST_LogFileID,*)

  SELECT CASE (TRIM(g_sRunStyle))
  CASE ("SIMU")
    WRITE(*,"(2X,A)") &
      "Running Style is Simulation!"
    WRITE(g_CREST_LogFileID,"(2X,A)") &
      "Running Style is Simulation!"

    WRITE(*,*)"  Model Start Time: "&
      //myDtoStr(g_StartDate, g_TimeMark)
    WRITE(*,*)"  Model End Time: " &
      //myDtoStr(g_EndDate, g_TimeMark)
    WRITE(g_CREST_LogFileID,*)"  Model Start Time: " &
      //myDtoStr(g_StartDate, g_TimeMark)
    WRITE(g_CREST_LogFileID,*)"  Model End Time: " &
      //myDtoStr(g_EndDate, g_TimeMark)
    WRITE(*,*)
    WRITE(g_CREST_LogFileID,*)

    CALL CREST_Simu()

  CASE ("CALI_SCEUA")
    WRITE(*,"(2X,A)") &
      "Running Style is Calibration using SCE-UA!"
    WRITE(g_CREST_LogFileID,"(2X,A)") &
      "Running Style is Calibration using SCE-UA!"

    CALL CREST_SCEUA()

  CASE ("REALTIME")
    WRITE(*,"(2X,A)") &
      "Running Style is Real Time!"
    WRITE(g_CREST_LogFileID,"(2X,A)") &
      "Running Style is Real Time!"

    CALL CREST_RealTime()

  CASE ("REPE")
    WRITE(*,"(2X,A)") &
      "Running Style is Return Period!"
    WRITE(g_CREST_LogFileID,"(2X,A)") &
      "Running Style is Return Period!"

    CALL CREST_ReturnPeriod()

  END SELECT

  WRITE(*,*)
  WRITE(g_CREST_LogFileID,*)

  !Get the end time and then display on Screen
  CALL UtiGetAndDisplayEndTime(0,IBDT,0)

  WRITE(*,*)"Project: "// TRIM(g_ProjectName) //" is finished!"
  ! WRITE(*,*)"Press any key to continue!"
  ! READ(*,*)

  CLOSE(g_CREST_LogFileID)
  STOP
  END PROGRAM CREST_Main

!============================================================
!
!============================================================
  SUBROUTINE GetNSCE(RObs,RSim,indexFrom,ITMax,NSCE)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT):: NSCE
  INTEGER :: i,indexFrom,indexL,ITMax,iCount
  DOUBLE PRECISION::RObs(0:ITMax-1),RSim(0:ITMax-1)
  DOUBLE PRECISION::RObs_Ave,RObs_Sum,RSim_Sum

  !print*, "RObs = ", RObs
  !print*, "RSim = ", RSim
  indexL=LBOUND(Robs,1)

  IF(indexFrom<indexL)THEN
    indexFrom=indexL
  END IF

  RObs_Ave=0.0
  iCount=0
  DO i=indexFrom,ITMax-1
    IF(RObs(i)<0 .OR. RSim(i)<0)THEN
      CYCLE
    END IF
    RObs_Ave=RObs_Ave+RObs(i)
    iCount=iCount+1
  END DO
  RObs_Ave=RObs_Ave/iCount

  RObs_Sum=0.0
  RSim_Sum=0.0
  DO i=indexFrom,ITMax-1
    IF(RObs(i)<0 .OR. RSim(i)<0)THEN
      CYCLE
    END IF
    RObs_Sum=RObs_Sum+(RObs(i)-RObs_Ave)**2.0
    RSim_Sum=RSim_Sum+(RObs(i)-RSim(i))**2.0
  END DO

  NSCE=1.0-RSim_Sum/RObs_Sum

  RETURN

  END SUBROUTINE GetNSCE

!============================================================
!
!============================================================
  SUBROUTINE GetBias(RObs,RSim,indexFrom,ITMax,Bias)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT):: Bias
  INTEGER :: i,indexFrom,indexL,ITMax,iCount
  DOUBLE PRECISION::RObs(0:ITMax-1),RSim(0:ITMax-1)
  DOUBLE PRECISION::RObs_Sum,RSim_Sum

  indexL=LBOUND(Robs,1)

  IF(indexFrom<indexL)THEN
    indexFrom=indexL
  END IF

  RObs_Sum=0.0
  RSim_Sum=0.0
  iCount=0
  DO i=indexFrom,ITMax-1
    IF(RObs(i)<0 .OR. RSim(i)<0)THEN
      CYCLE
    END IF
    RObs_Sum=RObs_Sum+RObs(i)
    RSim_Sum=RSim_Sum+RSim(i)
  END DO

  Bias=(RSim_Sum/RObs_Sum-1.0)*100

  RETURN

  END SUBROUTINE GetBias

!============================================================
!
!============================================================
  SUBROUTINE GetCC(RObs,RSim,indexFrom,ITMax,CC)
  ! orrelation coefficient
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT):: CC
  INTEGER :: i,indexFrom,indexL,ITMax,iCount
  DOUBLE PRECISION::RObs(0:ITMax-1),RSim(0:ITMax-1)
  DOUBLE PRECISION::RObs_Ave,RSim_Ave
  DOUBLE PRECISION::RAbove_Sum,RBelow_Sum
  DOUBLE PRECISION::RObs_Sum,RSim_Sum

  indexL=LBOUND(Robs,1)

  IF(indexFrom<indexL)THEN
    indexFrom=indexL
  END IF

  RObs_Ave=0.0
  RSim_Ave=0.0
  iCount=0
  DO i=indexFrom,ITMax-1
    IF(RObs(i)<0 .OR. RSim(i)<0)THEN
      CYCLE
    END IF
    RObs_Ave=RObs_Ave+RObs(i)
    RSim_Ave=RSim_Ave+RSim(i)
    iCount=iCount+1
  END DO
  RObs_Ave=RObs_Ave/iCount
  RSim_Ave=RSim_Ave/iCount

  RAbove_Sum=0.0
  RBelow_Sum=0.0
  RObs_Sum=0.0
  RSim_Sum=0.0
  DO i=indexFrom,ITMax-1
    IF(RObs(i)<0 .OR. RSim(i)<0)THEN
      CYCLE
    END IF
    RAbove_Sum=RAbove_Sum &
      +(RObs(i)-RObs_Ave)*(RSim(i)-RSim_Ave)
    RObs_Sum=RObs_Sum &
      +(RObs(i)-RObs_Ave)**2.0
    RSim_Sum=RSim_Sum &
      +(RSim(i)-RSim_Ave)**2.0
  END DO
  RBelow_Sum=sqrt(RObs_Sum*RSim_Sum)
  CC=RAbove_Sum/RBelow_Sum
  RETURN

  END SUBROUTINE GetCC


!============================================================
!
!============================================================
! Revised: Yanjun Gan, 03/12/2015

  SUBROUTINE ReadBasicFile(bIsError)
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE
  LOGICAL :: bIsError
  DOUBLE PRECISION,ALLOCATABLE :: AreaFact(:,:)

  ALLOCATE(g_DEM(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading DEM File!"
  CALL ReadMatrixFile(TRIM(g_BasicPath) // "DEM",g_DEM,  &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  ALLOCATE(g_FDR(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading FDR File!"
  CALL ReadMatrixFile_Int(TRIM(g_BasicPath) // "FDR", g_FDR, &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  ! Convert to Standard Flow Direction Map
  IF(COUNT(g_FDR==3)/=0 &
    .OR. COUNT(g_FDR==5)/=0 &
    .OR. COUNT(g_FDR==7)/=0)THEN
    WRITE(*,*)"  Converting DDM to FDR!"
    CALL ConvDDMToFDR()

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"FDR", &
      g_FDR, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)
  END IF

  ALLOCATE(g_FAC(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading FAC File!"
  CALL ReadMatrixFile_Int(TRIM(g_BasicPath) // "FAC", g_FAC, &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")

  ALLOCATE(g_GridArea(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(g_NextLen(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(g_NextC(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(g_NextR(0:g_NCols-1,0:g_NRows-1))

  WRITE(*,*)"  Reading GridArea File!"
  CALL ReadMatrixFile(TRIM(g_BasicPath) // "GridArea",g_GridArea,  &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  CALL AssignNextGroup(g_FDR,bIsError)

  IF(IsDebug==1)THEN
    !XianwuXue 20141122
    CALL WriteMatrixFile_Int(TRIM(g_ResultPath) // "NextC",g_NextC,  &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
    CALL WriteMatrixFile_Int(TRIM(g_ResultPath) // "NextR",g_NextR,  &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
  END IF

  WRITE(*,*)"  Reading AreaFact File!"
  ALLOCATE(AreaFact(0:g_NCols-1,0:g_NRows-1))
  CALL ReadMatrixFile(TRIM(g_BasicPath) // "AreaFact",AreaFact,  &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .FALSE.)THEN
    WHERE(g_Mask/=g_NoData_Value)
      g_GridArea=AreaFact*g_GridArea
    ENDWHERE
  END IF

!+ygan
  IF(IsDebug==1)THEN
    CALL WriteMatrixFile(TRIM(g_ResultPath)//"GridArea",g_GridArea,  &
      g_NCols,g_NRows,g_XLLCorner,g_YLLCorner,  &
      g_CellSize, g_NoData_Value,bIsError,g_BasicFormat)
  END IF
!.ygan

  !Read Stream Map
  ALLOCATE(g_Stream(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading Stream File!"
  CALL ReadMatrixFile_Int(TRIM(g_BasicPath) // "STREAM", g_Stream,  &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    CALL GetStream()

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"STREAM",  &
      g_Stream, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)
  END IF

  !Read Mask Map
  ALLOCATE(g_Mask(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading Mask File!"
  CALL ReadMatrixFile_Int(TRIM(g_BasicPath) // "MASK", g_Mask, &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    IF(g_tOutlet%HasOutlet)THEN
      !Get Mask Map by Outlet
      CALL GetMask(g_NCols,g_NRows,g_tOutlet%Col,g_tOutlet%Row, &
        g_NoData_Value,g_FDR,g_NextC,g_NextR, &
        g_Mask,InBasin)
    ELSE
      CALL GetMaskByNoData(g_NCols,g_NRows, &
        g_NoData_Value,g_DEM,g_Mask)
    END IF

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"MASK",  &
      g_Mask, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)
  END IF

  !Read Slope Map
  ALLOCATE(g_Slope(0:g_NCols-1,0:g_NRows-1))
  WRITE(*,*)"  Reading Slope File!"
  CALL ReadMatrixFile(TRIM(g_BasicPath) // "Slope",g_Slope,  &
    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN

    CALL CalSlope_Tiger(bIsError)

    IF(bIsError .EQV. .TRUE.)THEN
      CALL CalSlope_Tiger_XXW()
    END IF

    CALL WriteMatrixFile(TRIM(g_ResultPath)//"Slope",   &
      g_Slope, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)
  END IF
  WHERE(ABS(g_Slope-0.0)<0.000001)
    g_Slope=0.000001
  END WHERE


  RETURN
  END SUBROUTINE ReadBasicFile


!============================================================
!
!============================================================
  SUBROUTINE XXWReadRunoffObs(sFileNP,dblArray,bIsError)
  USE CREST_Project
  IMPLICIT NONE
  CHARACTER(*)::sFileNP
  LOGICAL :: bIsError
  LOGICAL :: fExist
  INTEGER :: fileid
  INTEGER :: error
  CHARACTER(LEN=200)::sTemp,sTemp2
  DOUBLE PRECISION :: dblTemp
  INTEGER :: dtTemp(1:6)
  CHARACTER(LEN=14)::strDate
  DOUBLE PRECISION :: dblArray(0:g_ITMax-1)
  INTEGER :: k

  dblArray=-9999

  bIsError=.FALSE.
  INQUIRE(FILE=TRIM(sFileNP), EXIST=fExist)
  IF(fExist .EQV. .FALSE.)THEN
    bIsError=.TRUE.
    RETURN
  END IF
  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0

  READ(fileid,*,IOSTAT=error,END=999)sTemp,sTemp2 !omit the header

  strDate=myDtoStr(g_StartDate, g_TimeMark)

  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error,END=999)sTemp,dblTemp
    IF(dblTemp<0)dblTemp=-9999
    CALL myStrtoD(sTemp,dtTemp,g_TimeMark)
    k=myDdiff(g_TimeMark, g_StartDate, dtTemp)/ g_TimeStep
    IF(k>=0 .AND. k<g_ITMax)THEN
      dblArray(k)=dblTemp
    END IF

  END DO
999 CONTINUE
  CLOSE(fileid)

  RETURN
  END SUBROUTINE XXWReadRunoffObs

!============================================================
!
!============================================================
  SUBROUTINE ReadICSFile()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_ICS

  IMPLICIT NONE
  INTEGER :: error
  LOGICAL :: bIsDIS,bIsError
  CHARACTER(LEN=200):: FileName,sParamName
  DOUBLE PRECISION :: dblTemp
  FileName=TRIM(g_ICSPath) // "InitialConditions.txt"

  ALLOCATE(g_WU0(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(g_SS0(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(g_SI0(0:g_NCols-1,0:g_NRows-1))
  g_WU0=g_NoData_Value
  g_SS0=g_NoData_Value
  g_SI0=g_NoData_Value

  sParamName="WU0"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ICSPath) // TRIM(sParamName), &
      g_WU0, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ICSFormat,"")
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_WU0=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="SS0"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ICSPath) // TRIM(sParamName),  &
      g_SS0, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ICSFormat,"")
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_SS0=dblTemp
    END WHERE
  END IF


  !-----------------------------------------------------------------------
  sParamName="SI0"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ICSPath) // TRIM(sParamName),  &
      g_SI0, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ICSFormat,"")
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_SI0=dblTemp
    END WHERE
  END IF

  RETURN
  END SUBROUTINE ReadICSFile

!============================================================
!
!============================================================
  SUBROUTINE ReadParamFile()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  INTEGER :: error
  LOGICAL :: bIsDIS,bIsError
  CHARACTER(LEN=200):: FileName,sParamName
  DOUBLE PRECISION :: dblTemp
  FileName=TRIM(g_ParamPath) // "Parameters.txt"

  CALL InitParamsType(g_tParams)

  ! Physical Parameters
  sParamName="Ksat"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%Ksat, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%Ksat=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%Ksat=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="RainFact"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",   &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),   &
      g_tParams%RainFact, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%RainFact=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%RainFact=dblTemp
    END WHERE
  END IF


  !-----------------------------------------------------------------------
  sParamName="WM"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%WM, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%WM=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%WM=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="B"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%B, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%B=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%B=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="IM"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%IM, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%IM=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%IM=dblTemp
    END WHERE
  END IF
  !-----------------------------------------------------------------------
  sParamName="KE"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%KE, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%KE=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%KE=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="coeM"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName), &
      g_tParams%coeM, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%coeM=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%coeM=dblTemp
    END WHERE
  END IF

  ! onceptual Parameters
  !-----------------------------------------------------------------------
  sParamName="expM"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%expM, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%expM=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%expM=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="coeR"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%coeR, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%coeR=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%coeR=dblTemp
    END WHERE
  END IF
  !-----------------------------------------------------------------------
  sParamName="coeS"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%coeS, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%coeS=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%coeS=dblTemp
    END WHERE
  END IF

  !-----------------------------------------------------------------------
  sParamName="KS"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type", &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%KS, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%KS=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%KS=dblTemp
    END WHERE
  END IF
  !-----------------------------------------------------------------------
  sParamName="KI"
  bIsDIS=XXWReadLineBln(FileName,TRIM(sParamName)//"Type",  &
    error,"Distributed")
  IF(bIsDIS)THEN
    CALL ReadMatrixFile(TRIM(g_ParamPath) // TRIM(sParamName),  &
      g_tParams%KI, &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_ParamFormat,"")
    WHERE(g_Mask==g_NoData_Value)
      g_tParams%KI=g_NoData_Value
    END WHERE
  ELSE
    dblTemp=XXWReadLineDbl(FileName,TRIM(sParamName), error)
    WHERE(g_Mask/=g_NoData_Value)
      g_tParams%KI=dblTemp
    END WHERE
  END IF

  RETURN
  END SUBROUTINE ReadParamFile



!============================================================
!
!============================================================
  SUBROUTINE CalSlope_Tiger(bIsError)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  INTEGER :: i,j,ii,jj
  DOUBLE PRECISION :: GM(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: GMValue,GMValueIn

  LOGICAL :: fExist,bIsError
  INTEGER :: fileid
  CHARACTER(LEN=200):: FileName
  DOUBLE PRECISION :: dblTemp
  FileName=TRIM(g_BasicPath) // "Slope.def"

  bIsError=.FALSE.
  INQUIRE(FILE=TRIM(FileName), EXIST=fExist)

  IF(fExist .EQV. .FALSE.)THEN
    bIsError=.TRUE.
    RETURN
  END IF
  CALL XXWGetFreeFile(fileid)
  OPEN(fileid,FILE=TRIM(FileName), FORM="FORMATTED")
  READ(fileid,*)GMValue
  CLOSE(fileid)

  ! GMValue=1.185728
  GM=GMValue

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      ii=g_NextR(j,i)
      jj=g_NextC(j,i)

      IF(g_Mask(j,i)/=g_NoData_Value)THEN
        IF(InBasin(jj,ii) .EQV. .FALSE.)THEN
          g_Slope(j,i)=GM(j,i)/g_NextLen(j,i)
        ELSE
          IF(g_DEM(j,i)>g_DEM(jj,ii))THEN
            g_Slope(j,i)=(g_DEM(j,i)-g_DEM(jj,ii)) &
              /g_NextLen(j,i)
          ELSE
            g_Slope(j,i) &
              =GM(j,i)/g_NextLen(j,i)
          END IF
        END IF
      ELSE
        g_Slope(j,i)=g_NoData_Value
      END IF
    END DO
  END DO

  RETURN
  END SUBROUTINE CalSlope_Tiger

!============================================================
!
!============================================================
  SUBROUTINE CalSlope_Tiger_XXW()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  INTEGER :: i,j,ii,jj
  DOUBLE PRECISION,EXTERNAL :: XXWCell_i_j_Slope_All

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1

      IF(g_Mask(j,i)/=g_NoData_Value)THEN
        ii=g_NextR(j,i)
        jj=g_NextC(j,i)
        IF(InBasin(jj,ii) .EQV. .FALSE.)THEN
          g_Slope(j,i)=XXWCell_i_j_Slope_All(i,j)

        ELSE
          IF(g_DEM(j,i)>g_DEM(jj,ii))THEN
            g_Slope(j,i)=(g_DEM(j,i)-g_DEM(jj,ii)) &
              /g_NextLen(j,i)
          ELSE ! If the direction map is correct, this condition will not happen
            g_Slope(j,i)=XXWCell_i_j_Slope_All(i,j)
          END IF
        END IF
        IF(g_Slope(j,i)<0)THEN
          g_Slope(j,i)=-g_Slope(j,i)
        END IF
        IF(ABS(g_Slope(j,i)-0.0)<0.000001)THEN
          g_Slope(j,i)=0.000001
        END IF
      ELSE
        g_Slope(j,i)=g_NoData_Value
      END IF

    END DO
  END DO

  RETURN
  END SUBROUTINE CalSlope_Tiger_XXW

!============================================================
!
!============================================================
  SUBROUTINE CalSlope_XXW()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  ! INTEGER :: i,j,ii,jj
  LOGICAL :: bGCS ,bIsNotExist ! Examine whether the GridArea.asc EXIST
  INTEGER :: j,i
  DOUBLE PRECISION :: LenSN, LenEW

  DOUBLE PRECISION, EXTERNAL :: XXWCell_i_j_Slope

  !IdentIFy whether it is a Geographic Coordinate System (GCS)
  ! or Projected Coordinate System (PCS)
  IF(g_XLLCorner>=-180.0 .AND. g_XLLCorner<=180 &
    .AND. g_YLLcorner>=-90 .AND. g_YLLCorner<=90)THEN
    bGCS=.TRUE.
  ELSE
    bGCS=.FALSE.
  END IF

  IF(bGCS .EQV. .TRUE.)THEN
    LenSN=g_CellSize*110574.0
  ELSE
    LenSN=g_CellSize
  END IF

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      IF(g_DEM(j,i)/=g_NoData_Value)THEN
        IF(bGCS .EQV. .TRUE.)THEN
          LenEW=g_YLLCorner+(g_NRows-i-0.5)*g_CellSize
          LenEW=LenSN*cos(LenEW*4.0*atan(1.0)/180.0)
        ELSE
          LenEW=g_CellSize
        END IF
        g_Slope(j,i)=XXWCell_i_j_Slope(g_NCols,g_NRows,i,j, &
          g_CellSize, g_NoData_Value,g_DEM,LenSN, LenEW)
        IF(g_Slope(j,i)<0)THEN
          g_Slope(j,i)=-g_Slope(j,i)
        END IF
        IF(g_Slope(j,i)==0)THEN
          g_Slope(j,i)=0.00000001
        END IF
      END IF
    END DO
  END DO

  RETURN
  END SUBROUTINE CalSlope_XXW

!============================================================
!
!============================================================
  DOUBLE PRECISION FUNCTION XXWCell_i_j_Slope_All(i,j)
  USE CREST_Basic
  IMPLICIT NONE

  INTEGER :: i,j
  DOUBLE PRECISION :: LenSN, LenEW
  DOUBLE PRECISION, EXTERNAL :: XXWCell_i_j_Slope
  LOGICAL :: bGCS

  IF(g_XLLCorner>=-180.0 .AND. g_XLLCorner<=180 &
    .AND. g_YLLcorner>=-90 .AND. g_YLLCorner<=90)THEN
    bGCS=.TRUE.
  ELSE
    bGCS=.FALSE.
  END IF

  IF(bGCS .EQV. .TRUE.)THEN
    LenSN=g_CellSize*110574.0
  ELSE
    LenSN=g_CellSize
  END IF

  IF(bGCS .EQV. .TRUE.)THEN
    LenEW=g_YLLCorner+(g_NRows-i-0.5)*g_CellSize
    LenEW=LenSN*cos(LenEW*4.0*atan(1.0)/180.0)
  ELSE
    LenEW=g_CellSize
  END IF

  XXWCell_i_j_Slope_All=XXWCell_i_j_Slope(g_NCols,g_NRows,i,j, &
    g_CellSize, g_NoData_Value,g_DEM,LenSN, LenEW)
  END FUNCTION XXWCell_i_j_Slope_All


!============================================================
!
!============================================================
  DOUBLE PRECISION FUNCTION XXWCell_i_j_Slope(NCols,NRows,i,j, &
    CellSize, NoData_Value,dblDEMMat,LenSN, LenEW)
  IMPLICIT NONE

  INTEGER :: NCols,NRows,i,j
  DOUBLE PRECISION :: CellSize, NoData_Value
  DOUBLE PRECISION :: dblDEMMat(0:NCols-1,0:NRows-1)
  DOUBLE PRECISION :: aa,bb,cc,dd,ee,ff,gg,hh,ii
  DOUBLE PRECISION :: LenSN, LenEW,dzdx,dzdy

  IF(dblDEMMat(j,i)==NoData_Value)THEN
    XXWCell_i_j_Slope=NoData_Value
    RETURN
  END IF

  ee=dblDEMMat(j,i)

  IF(i==0)THEN
    IF(j==0)THEN
      aa=ee
      bb=ee
      cc=ee
      dd=ee

      ff=dblDEMMat(j+1,i)
      gg=ee
      hh=dblDEMMat(j,i+1)
      ii=dblDEMMat(j+1,i+1)
    ELSEIF(j==NCols-1)THEN
      aa=ee
      bb=ee
      cc=ee
      dd=dblDEMMat(j-1,i)

      ff=ee
      gg=dblDEMMat(j-1,i+1)
      hh=dblDEMMat(j,i+1)
      ii=ee
    ELSE
      aa=ee
      bb=ee
      cc=ee
      dd=dblDEMMat(j-1,i)

      ff=dblDEMMat(j+1,i)
      gg=dblDEMMat(j-1,i+1)
      hh=dblDEMMat(j,i+1)
      ii=dblDEMMat(j+1,i+1)
    END IF

  ELSEIF(i==NRows-1)THEN
    IF(j==0)THEN
      aa=ee
      bb=dblDEMMat(j,i-1)
      cc=dblDEMMat(j+1,i-1)
      dd=ee

      ff=dblDEMMat(j+1,i)
      gg=ee
      hh=ee
      ii=ee
    ELSEIF(j==NCols-1)THEN
      aa=dblDEMMat(j-1,i-1)
      bb=dblDEMMat(j,i-1)
      cc=ee
      dd=dblDEMMat(j-1,i)

      ff=ee
      gg=ee
      hh=ee
      ii=ee
    ELSE
      aa=dblDEMMat(j-1,i-1)
      bb=dblDEMMat(j,i-1)
      cc=dblDEMMat(j+1,i-1)
      dd=dblDEMMat(j-1,i)

      ff=dblDEMMat(j+1,i)
      gg=ee
      hh=ee
      ii=ee
    END IF
  ELSE
    IF(j==0)THEN
      aa=ee
      bb=dblDEMMat(j,i-1)
      cc=dblDEMMat(j+1,i-1)
      dd=ee

      ff=dblDEMMat(j+1,i)
      gg=ee
      hh=dblDEMMat(j,i+1)
      ii=dblDEMMat(j+1,i+1)
    ELSEIF(j==NCols-1)THEN
      aa=dblDEMMat(j-1,i-1)
      bb=dblDEMMat(j,i-1)
      cc=ee
      dd=dblDEMMat(j-1,i)

      ff=ee
      gg=dblDEMMat(j-1,i+1)
      hh=dblDEMMat(j,i+1)
      ii=ee
    ELSE
      aa=dblDEMMat(j-1,i-1)
      bb=dblDEMMat(j,i-1)
      cc=dblDEMMat(j+1,i-1)
      dd=dblDEMMat(j-1,i)

      ff=dblDEMMat(j+1,i)
      gg=dblDEMMat(j-1,i+1)
      hh=dblDEMMat(j,i+1)
      ii=dblDEMMat(j+1,i+1)
    END IF
  END IF

  IF(aa==NoData_Value)THEN
    aa=ee
  END IF

  IF(bb==NoData_Value)THEN
    bb=ee
  END IF

  IF(cc==NoData_Value)THEN
    cc=ee
  END IF

  IF(dd==NoData_Value)THEN
    dd=ee
  END IF

  IF(ee==NoData_Value)THEN
    ee=ee
  END IF

  IF(ff==NoData_Value)THEN
    ff=ee
  END IF

  IF(gg==NoData_Value)THEN
    gg=ee
  END IF

  IF(hh==NoData_Value)THEN
    hh=ee
  END IF

  IF(ii==NoData_Value)THEN
    ii=ee
  END IF

  dzdx=(cc+2.0*ff+ii)-(aa+2.0*dd+gg)
  dzdx=dzdx/(8.0*LenEW)

  dzdy=(gg+2.0*hh+ii)-(aa+2.0*bb+cc)
  dzdy=dzdy/(8.0*LenSN)

  XXWCell_i_j_Slope=sqrt(dzdx**2+dzdy**2)
  ! XXWCell_i_j_Slope=XXWCell_i_j_Slope*4.0*atan(1.0)/180.0
  END FUNCTION XXWCell_i_j_Slope

!============================================================
!
!============================================================
  SUBROUTINE ConvDDMToFDR()
  USE CREST_Project
  USE CREST_Basic
  IMPLICIT NONE

  INTEGER :: i,j,ii,jj
  INTEGER :: numNoData,iRowNoData,iColNoData

  DO i=0,g_NRows-1
    DO j=0,g_NCols-1
      IF(g_FDR(j,i)==g_NoData_Value)THEN
        CYCLE
      END IF
      SELECT CASE(g_FDR(j,i))
      CASE(1)
        g_FDR(j,i)=64
      CASE(2)
        g_FDR(j,i)=128
      CASE(8)
        g_FDR(j,i)=32
      CASE(3)
        g_FDR(j,i)=1
      CASE(4)
        g_FDR(j,i)=2
      CASE(5)
        g_FDR(j,i)=4
      CASE(6)
        g_FDR(j,i)=8
      CASE(7)
        g_FDR(j,i)=16
      CASE(0) !Outlet location
        numNoData=0
        iRowNoData=g_NoData_Value
        iColNoData=g_NoData_Value
        DO ii=-1,1
          DO jj=-1,1
            IF((ii==0).AND.(jj==0))THEN
              CYCLE
            END IF
            IF(InBasin(j+jj,i+ii))THEN
              IF(g_FDR(j+jj,i+ii)/=g_NoData_Value)THEN
              ELSE
                numNoData=numNoData+1
                iRowNoData=ii
                iColNoData=jj
              END IF
            ELSE
              numNoData=numNoData+1
              iRowNoData=ii
              iColNoData=jj
            END IF
          END DO
        END DO
        IF(numNoData==0)THEN
          g_FDR(j,i)=256 ! This cell is a sink
        ELSEIF(numNoData==8)THEN
          g_FDR(j,i)=1 ! This cell is an islAND
        ELSE
          IF(iRowNoData==-1)THEN
            IF(iColNoData==-1)THEN
              g_FDR(j,i)=32
            ELSEIF(iColNoData==0)THEN
              g_FDR(j,i)=64
            ELSE
              g_FDR(j,i)=128
            END IF
          ELSEIF(iRowNoData==0)THEN
            IF(iColNoData==-1)THEN
              g_FDR(j,i)=16
            ELSEIF(iColNoData==0)THEN
              WRITE(*,*)"Wrong!!!!!!!!!!"
            ELSE
              g_FDR(j,i)=1
            END IF
          ELSE
            IF(iColNoData==-1)THEN
              g_FDR(j,i)=8
            ELSEIF(iColNoData==0)THEN
              g_FDR(j,i)=4
            ELSE
              g_FDR(j,i)=2
            END IF
          END IF
        END IF
      END SELECT
    END DO
  END DO
  RETURN
  END SUBROUTINE ConvDDMToFDR

!============================================================
!
!============================================================
  SUBROUTINE AssignNextGroup(FDR,bIsNotExist)
  USE CREST_Project
  USE CREST_Basic
  IMPLICIT NONE
  LOGICAL :: bGCS ,bIsNotExist ! Examine whether the GridArea.asc exist
  INTEGER :: FDR(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: j,i
  DOUBLE PRECISION :: LenSN, LenEW, LenCross

  ! Identify whether it is a Geographic Coordinate System (GCS)
  ! or Projected Coordinate System (PCS)
  IF(g_XLLCorner>=-180.0 .AND. g_XLLCorner<=180 &
    .AND. g_YLLcorner>=-90 .AND. g_YLLCorner<=90)THEN
    bGCS=.TRUE.
  ELSE
    bGCS=.FALSE.
  END IF

  IF(bGCS .EQV. .TRUE.)THEN
    LenSN=g_CellSize*110574.0
  ELSE
    LenSN=g_CellSize
  END IF

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      IF(FDR(j,i)/=g_NoData_Value)THEN
        IF(bGCS .EQV. .TRUE.)THEN
          LenEW=g_YLLCorner+(g_NRows-i-0.5)*g_CellSize
          LenEW=LenSN*cos(LenEW*4.0*atan(1.0)/180.0)
        ELSE
          LenEW=g_CellSize
        END IF
        LenCross=sqrt(LenEW**2+LenSN**2)

        IF(bIsNotExist .EQV. .TRUE.)THEN
          g_GridArea(j,i)=LenSN*LenEW/1.0e6 ! Convert to km^2
        END IF

        SELECT CASE (FDR(j,i))
        CASE (64)
          g_NextR(j,i)=i-1
          g_NextC(j,i)=j
          g_NextLen(j,i)=LenSN
        CASE (128)
          g_NextR(j,i)=i-1
          g_NextC(j,i)=j+1
          g_NextLen(j,i)=LenCross
        CASE (1)
          g_NextR(j,i)=i
          g_NextC(j,i)=j+1
          g_NextLen(j,i)=LenEW
        CASE (2)
          g_NextR(j,i)=i+1
          g_NextC(j,i)=j+1
          g_NextLen(j,i)=LenCross
        CASE (4)
          g_NextR(j,i)=i+1
          g_NextC(j,i)=j
          g_NextLen(j,i)=LenSN
        CASE (8)
          g_NextR(j,i)=i+1
          g_NextC(j,i)=j-1
          g_NextLen(j,i)=LenCross
        CASE (16)
          g_NextR(j,i)=i
          g_NextC(j,i)=j-1
          g_NextLen(j,i)=LenEW
        CASE (32)
          g_NextR(j,i)=i-1
          g_NextC(j,i)=j-1
          g_NextLen(j,i)=LenCross
        CASE (256) !Outlet
          g_NextR(j,i)=g_NRows
          g_NextC(j,i)=g_NCols
          g_NextLen(j,i)=LenCross
          CASE DEFAULT
          WRITE(*,*) "Something is wrong " &
            // "in your Flow Direction Map!"
        END SELECT

        IF(InBasin(g_NextC(j,i), &
          g_NextR(j,i)) .EQV. .FALSE.)THEN
         ! g_NextC(j,i)=g_NoData_Value
         ! g_NextR(j,i)=g_NoData_Value
          g_NextLen(j,i)=LenSN
        END IF
      ELSE
        g_NextC(j,i)=g_NoData_Value
        g_NextR(j,i)=g_NoData_Value
        g_NextLen(j,i)=g_NoData_Value
        IF(bIsNotExist .EQV. .TRUE.)THEN
          g_GridArea(j,i)=g_NoData_Value
        END IF
      END IF
    END DO
  END DO

  RETURN
  END SUBROUTINE AssignNextGroup



!============================================================
!
!============================================================
  SUBROUTINE URWORD(LINE,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)

  CHARACTER*(*) LINE
  CHARACTER*20 STRING
  CHARACTER*30 RW
  CHARACTER*1 TAB
  INTEGER :: N,IOUT,IN
  ! ------------------------------------------------------------------
  TAB=CHAR(9)
  !

  LINLEN=LEN(LINE)
  LINE(LINLEN:LINLEN)=' '
  ISTART=LINLEN
  ISTOP=LINLEN
  LINLEN=LINLEN-1
  IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
  !
  !2------Find start of word, which is indicated by first CHARACTER that
  !2------is not a blank, a comma, or a tab.
  DO 10 I=ICOL,LINLEN
    IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.',' &
      .AND. LINE(I:I).NE.TAB) GO TO 20
10 CONTINUE
  ICOL=LINLEN+1
  GO TO 100
  !
  !3------Found start of word. Look for end.
  !3A-----When word is quoted, only a quote can terminate it.
20 IF(LINE(I:I).EQ.'''') THEN
    I=I+1
    IF(I.LE.LINLEN) THEN
      DO 25 J=I,LINLEN
        IF(LINE(J:J).EQ.'''') GO TO 40
25      CONTINUE
    END IF
    !
    !3B-----When word is not quoted, space, comma, or tab will terminate.
  ELSE
    DO 30 J=I,LINLEN
      IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.',' &
        .OR. LINE(J:J).EQ.TAB) GO TO 40
30    CONTINUE
  END IF
  !
  !3C-----End of line without finding END of word; set END of word to
  !3C-----END of line.
  J=LINLEN+1
  !
  !4------Found end of word; set J to point to last character in word and
  !-------set icol to point to location for scanning for another word.
40 ICOL=J+1
  J=J-1
  IF(J.LT.I) GO TO 100
  ISTART=I
  ISTOP=J
  !
  !5------Convert word to upper case and return if NCODE is 1.
  IF(NCODE.EQ.1) THEN
    IDIFF=ICHAR('a')-ICHAR('A')
    DO 50 K=ISTART,ISTOP
      IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z') &
        LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50  CONTINUE
    RETURN
  END IF
  !
  !6------Convert word to a number IF rEQuested.
100 IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
    RW=' '
    L=30-ISTOP+ISTART
    IF(L.LT.1) GO TO 200
    RW(L:30)=LINE(ISTART:ISTOP)
    IF(NCODE.EQ.2) READ(RW,'(I30)',ERR=200) N
    IF(NCODE.EQ.3) READ(RW,'(F30.0)',ERR=200) R
  END IF
  RETURN
  !
  !7------Number conversion error.
200 IF(NCODE.EQ.3) THEN
    STRING= 'A REAL NUMBER'
    L=13
  ELSE
    STRING= 'AN INTEGER'
    L=10
  END IF
  !
  !7A-----If output unit is negative, set last character of string to 'E'.
  IF(IOUT.LT.0) THEN
    N=0
    R=0.
    LINE(LINLEN+1:LINLEN+1)='E'
    RETURN
    !
    !7B-----If output unit is positive; write a message to output unit.
  ELSE IF(IOUT.GT.0) THEN
    IF(IN.GT.0) THEN
      WRITE(IOUT,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
    ELSE
      WRITE(IOUT,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
    END IF
201 FORMAT(1X,/1X,'FILE UNIT ',I4,' : ERROR CONVERTING "',A, &
      '" TO ',A,' IN LINE:',/1X,A)
202 FORMAT(1X,/1X,'KEYBOARD INPUT : ERROR CONVERTING "',A, &
      '" TO ',A,' IN LINE:',/1X,A)
    !
    !7C-----If output unit is 0; write a message to default output.
  ELSE
    IF(IN.GT.0) THEN
      WRITE(*,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
    ELSE
      WRITE(*,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
    END IF
  END IF
  !
  END SUBROUTINE URWORD

!============================================================
!
!============================================================
  SUBROUTINE UPCASE(WORD)
  ! ******************************************************************
  ! CONVERT A CHARACTER STRING TO ALL UPPER CASE
  ! ******************************************************************
  ! SPECIFICATIONS:
  ! ------------------------------------------------------------------
  CHARACTER WORD*(*)
  !
  !1------Compute the difference between lowercase and uppercase.
  L = LEN(WORD)
  IDIFF=ICHAR('a')-ICHAR('A')
  !
  !2------Loop through the string and convert any lowercase characters.
  DO 10 K=1,L
    IF(WORD(K:K).GE.'a' .AND. WORD(K:K).LE.'z') &
      WORD(K:K)=CHAR(ICHAR(WORD(K:K))-IDIFF)
10 CONTINUE
  !
  !3------RETURN.
  RETURN
  END SUBROUTINE UPCASE

!============================================================
!
!============================================================
  SUBROUTINE GetMask(NCols,NRows,jCol,iRow,NoData_Value,FDR, &
    NextC,NextR,MaskOut,InBasin)
  IMPLICIT NONE
  INTEGER :: iRow,jCol !Objective Row/Col
  INTEGER :: NCols,NRows
  DOUBLE PRECISION :: NoData_Value
  INTEGER :: NextC(0:NCols-1,0:NRows-1)
  INTEGER :: NextR(0:NCols-1,0:NRows-1)
  INTEGER :: FDR(0:NCols-1,0:NRows-1)
  INTEGER, INTENT(OUT) :: MaskOut(0:NCols-1,0:NRows-1)
  INTEGER :: i,j,iiA,iiB,jjA,jjB
  LOGICAL,EXTERNAL :: InBasin

  MaskOut=NoData_Value
  MaskOut(jCol,iRow)=1
  DO i=0, NRows-1
    DO j=0, NCols-1
      IF(FDR(j,i)==NoData_Value)THEN
        CYCLE
      END IF
      iiA=NextR(j,i)
      jjA=NextC(j,i)
      DO WHILE(InBasin(jjA,iiA) .EQV. .TRUE.)
        IF(FDR(jjA,iiA)==NoData_Value)THEN
          EXIT
        END IF
        IF(iiA==iRow .AND. jjA==jCol)THEN
          MaskOut(j,i)=1
          EXIT
        END IF
        iiB=iiA
        jjB=jjA
        iiA=NextR(jjB,iiB)
        jjA=NextC(jjB,iiB)
      END DO
    END DO
  END DO
  RETURN
  END SUBROUTINE GetMask

!============================================================
!
!============================================================
  SUBROUTINE GetMaskByNoData(NCols,NRows,NoData_Value, &
    DEM,MaskOut)
  IMPLICIT NONE
  INTEGER :: NCols,NRows
  DOUBLE PRECISION :: NoData_Value
  DOUBLE PRECISION :: DEM(0:NCols-1,0:NRows-1)
  INTEGER, INTENT(OUT) :: MaskOut(0:NCols-1,0:NRows-1)
  INTEGER :: i,j

  DO i=0, NRows-1
    DO j=0, NCols-1
      IF(DEM(j,i)==NoData_Value)THEN
        MaskOut(j,i)=NoData_Value
      ELSE
        MaskOut(j,i)=1
      END IF
    END DO
  END DO
  RETURN
  END SUBROUTINE GetMaskByNoData

!============================================================
!
!============================================================
  SUBROUTINE GetStream()
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE
  DOUBLE PRECISION :: TH
  LOGICAL :: fExist

  INTEGER :: i,j
  INTEGER :: fileid
  CHARACTER(LEN=200):: FileName
  DOUBLE PRECISION :: dblTemp
  FileName=TRIM(g_BasicPath) // "Stream.def"

  INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
  IF(fExist .EQV. .FALSE.)THEN
    WRITE(*,*)"You should give a Stream Map or 'Stream.def' file!"
    WRITE(g_CREST_LogFileID,*) &
      "You should give a Stream Map or 'Stream.def' file!"
    STOP
  END IF
  CALL XXWGetFreeFile(fileid)
  OPEN(fileid,FILE=TRIM(FileName), FORM="FORMATTED")
  READ(fileid,*)TH
  CLOSE(fileid)

  DO i=0,g_NRows-1
    DO j=0,g_NCols-1
      IF(g_FAC(j,i)>TH)THEN
        g_Stream(j,i)=1
      ELSE
        g_Stream(j,i)=g_NoData_Value
      END IF
    END DO
  END DO

  RETURN
  END SUBROUTINE GetStream

!============================================================
!
!============================================================
  LOGICAL FUNCTION InBasin(jCol,iRow)
  USE CREST_Basic
  IMPLICIT NONE
  INTEGER :: iRow,jCol
  IF(jCol<0 &
    .OR. jCol>=g_NCols &
    .OR.iRow<0 &
    .OR. iRow>=g_NRows)THEN
    InBasin=.FALSE.
  ELSE
    ! IF(g_Mask(jCol,iRow)==g_NoData_Value)THEN
    ! InBasin=.FALSE.
    ! ELSE
    InBasin=.TRUE.
    ! END IF
  END IF
  RETURN
  END FUNCTION InBasin

!============================================================
!
!============================================================
  LOGICAL FUNCTION InBasin2(jCol,iRow)
  USE CREST_Basic
  IMPLICIT NONE
  INTEGER :: iRow,jCol
  IF(jCol<0 &
    .OR. jCol>=g_NCols &
    .OR.iRow<0 &
    .OR. iRow>=g_NRows)THEN
  InBasin2=.FALSE.
  ELSE
    ! IF(g_Mask(jCol,iRow)==g_NoData_Value)THEN
    ! InBasin=.FALSE.
    ! ELSE
    InBasin2=.TRUE.
    ! END IF
  END IF
  RETURN
  END FUNCTION InBasin2

!============================================================
!
!============================================================
! Revised: Yanjun Gan, 06/17/2015

  SUBROUTINE CREST_Simu()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_ICS
  USE CREST_Param
  USE CREST_Calib_SCEUA

  IMPLICIT NONE
  DOUBLE PRECISION :: NSCE, Bias, CC
  INTEGER :: ITMax
  INTEGER :: i,j,k,ii,jj,iii,jjj,iii2,jjj2
!+ygan
  INTEGER :: l, g_ProjectNameLength
  DOUBLE PRECISION, ALLOCATABLE :: Raina(:,:),PETa(:,:)
!.ygan
  INTEGER :: dtTemp(1:6)
  LOGICAL :: bIsError
  CHARACTER(14):: strDate
  CHARACTER(50):: sTemp
  INTEGER, ALLOCATABLE::OneRowA(:,:),OneRowB(:,:)
  INTEGER, ALLOCATABLE::OneColA(:,:),OneColB(:,:)
  INTEGER, ALLOCATABLE::TwoRowA(:,:),TwoRowB(:,:)
  INTEGER, ALLOCATABLE::TwoColA(:,:),TwoColB(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: OnePerA(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: OnePerB(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: TwoPerA(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: TwoPerB(:,:)

  DOUBLE PRECISION, ALLOCATABLE :: Rain(:,:),PET(:,:),EPot(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: W0(:,:),SS0(:,:),SI0(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: W(:,:),ExcS(:,:),ExcI(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: RI(:,:),RS(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: EAct(:,:),Runoff(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: Rain0(:,:),PET0(:,:)

  DOUBLE PRECISION :: Roff_s_tmp,Roff_i_tmp
  !-----------------------------------------------
  ALLOCATE(OneRowA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(OneRowB(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(OneColA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(OneColB(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(TwoRowA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(TwoRowB(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(TwoColA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(TwoColB(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(OnePerA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(OnePerB(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(TwoPerA(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(TwoPerB(0:g_NCols-1,0:g_NRows-1))

  !---------------------------------------------------
  ALLOCATE(Rain0(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(PET0(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(Rain(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(PET(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(EPot(0:g_NCols-1,0:g_NRows-1))

!+ygan
  ALLOCATE(Raina(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(PETa(0:g_NCols-1,0:g_NRows-1))
!.ygan

  ALLOCATE(W0(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(SS0(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(SI0(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(W(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(ExcS(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(ExcI(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(EAct(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(Runoff(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(RS(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(RI(0:g_NCols-1,0:g_NRows-1))

  SS0=g_SS0
  SI0=g_SI0

  Rain=g_NoData_Value
  PET=g_NoData_Value
  EPot=g_NoData_Value !Initial Potential Evapotranspiration
  EAct=g_NoData_Value
  W=g_NoData_Value
  Runoff=g_NoData_Value
  ExcS=g_NoData_Value
  ExcI=g_NoData_Value
  RS=g_NoData_Value
  RI=g_NoData_Value

  WHERE(g_Mask/=g_NoData_Value)
    W0=g_WU0*g_tParams%WM/100.0
    Runoff=0.0
  ELSEWHERE
    W0=g_NoData_Value
    Runoff=g_NoData_Value
  END WHERE

  !
  g_ITMax=myDdiff(g_TimeMark, g_StartDate, g_EndDate)/ g_TimeStep + 1
  ITMax=g_ITMax

  IF(g_RunStyle/=2)THEN
    IF(g_tOutlet%HasOutlet .EQV. .TRUE.)THEN
      ALLOCATE(g_tOutlet%Rain (0:ITMax-1))
      ALLOCATE(g_tOutlet%PET  (0:ITMax-1))
      ALLOCATE(g_tOutlet%EPot (0:ITMax-1))
      ALLOCATE(g_tOutlet%EAct (0:ITMax-1))
      ALLOCATE(g_tOutlet%W    (0:ITMax-1))
      ALLOCATE(g_tOutlet%SM   (0:ITMax-1))
      ALLOCATE(g_tOutlet%R    (0:ITMax-1))
      ALLOCATE(g_tOutlet%ExcS (0:ITMax-1))
      ALLOCATE(g_tOutlet%ExcI (0:ITMax-1))
      ALLOCATE(g_tOutlet%RS   (0:ITMax-1))
      ALLOCATE(g_tOutlet%RI   (0:ITMax-1))
    END IF

    DO k=0,g_NOutPixs-1
      ALLOCATE(g_tOutPix(k)%Rain  (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%PET   (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%EPot  (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%EAct  (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%W     (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%SM    (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%R     (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%ExcS  (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%ExcI  (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%RS    (0:ITMax-1))
      ALLOCATE(g_tOutPix(k)%RI    (0:ITMax-1))
    END DO
  END IF

  CALL CREST_RunoffRoutePre( &
    OneRowA,OneColA,OnePerA, &
    OneRowB,OneColB,OnePerB, &
    TwoRowA,TwoColA,TwoPerA, &
    TwoRowB,TwoColB,TwoPerB)


!--------------------------------------------------
! Load the State Data
!--------------------------------------------------
  IF(g_LoadState .EQV. .TRUE.)THEN !.AND. g_RunStyle/=2

    dtTemp=g_StartDate ! Modified by Xianwu Xue 2011.3.24
    CALL myDadd(g_TimeMark,-g_TimeStep,dtTemp)! Modified by Xianwu Xue 2011.3.24

    strDate=myDtoStr(dtTemp, g_TimeMark)! Modified by Xianwu Xue 2011.3.24

    WRITE(*,*) "Loading States at date: " // strDate
    WRITE(g_CREST_LogFileID,*)"Loading States at date: " // strDate

    CALL LoadStates(strDate,W0,SS0,SI0,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN ! File does not exist!
      WRITE(*,"(1X,A)") "*** Something wrong in your States file!"
      WRITE(g_CREST_LogFileID,"(1X,A)") "*** Something wrong in your States file!"

      RETURN
    END IF
  END IF

!--------------------------------------------------
! Begin the calculation (time->row->column)
!--------------------------------------------------
  dtTemp=g_StartDate
  DO k=0, ITMax-1   ! time loop
    strDate=myDtoStr(dtTemp, g_TimeMark)

    IF(g_RunStyle/=2)THEN  ! y(year);m(month);d(day);h(hour);u(minute);s(second)
      IF(TRIM(g_TimeMark)=="y")THEN
        sTemp=strDate(1:4)
      ELSEIF(TRIM(g_TimeMark)=="m")THEN
        sTemp=strDate(1:4)//"-"//strDate(5:6)
      ELSEIF(TRIM(g_TimeMark)=="d")THEN
        sTemp=strDate(1:4)//"-"//strDate(5:6)//"-"//strDate(7:8)
      ELSEIF(TRIM(g_TimeMark)=="h")THEN
        sTemp=strDate(1:4)//"-"//strDate(5:6)//"-"//strDate(7:8)&
          //" "//strDate(9:10)//":"//"00"//":"//"00"
      ELSEIF(TRIM(g_TimeMark)=="u")THEN
        sTemp=strDate(1:4)//"-"//strDate(5:6)//"-"//strDate(7:8) &
          //" "//strDate(9:10)//":"//strDate(11:12)//":"//"00"
      ELSEIF(TRIM(g_TimeMark)=="s")THEN
        sTemp=strDate(1:4)//"-"//strDate(5:6)//"-"//strDate(7:8) &
          //" "//strDate(9:10)//":"//strDate(11:12)//":"//strDate(13:14)
      ENDIF

      WRITE(*,*) k+1," ", TRIM(sTemp)
      WRITE(g_CREST_LogFileID,*) k+1," ", TRIM(sTemp)
    END IF


!+ygan
    g_ProjectNameLength = LEN(TRIM(g_ProjectName))
    IF (g_ProjectName(g_ProjectNameLength-2:).EQ. "acc") THEN
      dtTemp(3) = dtTemp(3) - 1   ! a day before

      DO l = 0, 23  ! daily accumulated
        CALL myDAdd(g_TimeMark,g_TimeStep,dtTemp)
        !print*,l+1,": ",dtTemp
        strDate=myDtoStr(dtTemp, g_TimeMark)
        WRITE(*,'(14X,A5,I2,A2,A10)'),"Hour ",l+1,": ",TRIM(strDate)

        !WRITE(*,*)"Reading Rain File!"
        CALL ReadMatrixFile(TRIM(g_RainPath), &
          Rain0, &
          g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
          g_CellSize,g_NoData_Value,bIsError,g_RainFormat, &
          TRIM(strDate))


        IF(bIsError .EQV. .TRUE.)THEN
          IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
          ELSE
            WRITE(*,"(15X,A)")"---Rain Data File Missed!"
            !WHERE(Rain>0) !Added by Dr. Xianwu Xue 2011.5.9
            !  Rain=Rain/g_TimeStep
            !END WHERE
          END IF
        ELSE
        END IF

        WHERE(Rain0>0)
          Rain=Rain0*g_TimeStep
        ELSEWHERE(Rain0<=0)
          Rain=0.0
        END WHERE

        !WRITE(*,*)"Reading PET File!"
        CALL ReadMatrixFile(TRIM(g_PETPath), &
          PET0, &
          g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
          g_CellSize,g_NoData_Value,bIsError,g_PETFormat, &
          TRIM(strDate))
        IF(bIsError .EQV. .TRUE.)THEN
          IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
          ELSE
            WRITE(*,"(15X,A)")"---PET Data File Missed!"
            !WHERE(PET>0) !Added by Dr. Xianwu Xue 2011.5.9
            ! PET=PET/g_TimeStep
            !END WHERE
          END IF
        ELSE
        END IF

        WHERE(PET0>0)
          PET=PET0*g_TimeStep
        ELSEWHERE(PET0<=0)
          PET=0.0
        END WHERE

        Raina = Raina + Rain
        PETa  = PETa  + PET
      END DO

      Rain = Raina
      PET  = PETa

    ELSE
      !WRITE(*,*)"Reading Rain File!"
      CALL ReadMatrixFile(TRIM(g_RainPath), &
        Rain0, &
        g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
        g_CellSize,g_NoData_Value,bIsError,g_RainFormat, &
        TRIM(strDate))


      IF(bIsError .EQV. .TRUE.)THEN
        IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
        ELSE
          WRITE(*,"(15X,A)")"---Rain Data File Missed!"
          !WHERE(Rain>0) !Added by Dr. Xianwu Xue 2011.5.9
          ! Rain=Rain/g_TimeStep
          !END WHERE
        END IF
      ELSE
      END IF

      WHERE(Rain0>0)
        Rain=Rain0*g_TimeStep
      ELSEWHERE(Rain0<=0)
        Rain=0.0
      END WHERE

      !WRITE(*,*)"Reading PET File!"
      CALL ReadMatrixFile(TRIM(g_PETPath), &
        PET0, &
        g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
        g_CellSize,g_NoData_Value,bIsError,g_PETFormat, &
        TRIM(strDate))
      IF(bIsError .EQV. .TRUE.)THEN
        IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
        ELSE
          WRITE(*,"(15X,A)")"---PET Data File Missed!"
          !WHERE(PET>0) !Added by Dr. Xianwu Xue 2011.5.9
          ! PET=PET/g_TimeStep
          !END WHERE
        END IF
      ELSE
      END IF

      WHERE(PET0>0)
        PET=PET0*g_TimeStep
      ELSEWHERE(PET0<=0)
        PET=0.0
      END WHERE

    END IF
!.ygan

!-ygan
  !  !WRITE(*,*)"Reading Rain File!"
  !  CALL ReadMatrixFile(TRIM(g_RainPath), &
  !    Rain0, &
  !    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
  !    g_CellSize,g_NoData_Value,bIsError,g_RainFormat, &
  !    TRIM(strDate))


  !  IF(bIsError .EQV. .TRUE.)THEN
  !    IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
  !    ELSE
  !      WRITE(*,"(15X,A)")"---Rain Data File Missed!"
  !      !WHERE(Rain>0) !Added by Dr. Xianwu Xue 2011.5.9
  !      ! Rain=Rain/g_TimeStep
  !      !END WHERE
  !    END IF
  !  ELSE
  !  END IF

  !  WHERE(Rain0>0)
  !    Rain=Rain0*g_TimeStep
  !  ELSEWHERE(Rain0<=0)
  !    Rain=0.0
  !  END WHERE

  !  !WRITE(*,*)"Reading PET File!"
  !  CALL ReadMatrixFile(TRIM(g_PETPath), &
  !    PET0, &
  !    g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
  !    g_CellSize,g_NoData_Value,bIsError,g_PETFormat, &
  !    TRIM(strDate))
  !  IF(bIsError .EQV. .TRUE.)THEN
  !    IF(g_RunStyle==2)THEN !Added by Dr. Xianwu Xue 2011.4.16
  !    ELSE
  !      WRITE(*,"(15X,A)")"---PET Data File Missed!"
  !      !WHERE(PET>0) !Added by Dr. Xianwu Xue 2011.5.9
  !      ! PET=PET/g_TimeStep
  !      !END WHERE
  !    END IF
  !  ELSE
  !  END IF

  !  WHERE(PET0>0)
  !    PET=PET0*g_TimeStep
  !  ELSEWHERE(PET0<=0)
  !    PET=0.0
  !  END WHERE
!.ygan

!--------------------------------------------------
! Runoff generation
!--------------------------------------------------
    DO i=0, g_NRows-1   ! row loop
      DO j=0, g_NCols-1   ! column loop
        !if(j==g_tOutlet%Col .and. i==g_tOutlet%Row) print*,g_tParams%RainFact(j,i), &
        !  g_tParams%Ksat(j,i),g_tParams%WM(j,i),  g_tParams%B(j,i),   g_tParams%IM(j,i), &
        !  g_tParams%KE(j,i),  g_tParams%coeM(j,i),g_tParams%expM(j,i),g_tParams%coeR(j,i), &
        !  g_tParams%coeS(j,i),g_tParams%KS(j,i),  g_tParams%KI(j,i)

        IF(g_Mask(j,i)==g_NoData_Value)THEN
          CYCLE
        END IF

        ! Convert unit of Rain and PET from mm/dt to mm
        !IF(Rain0(j,i)>0.0)THEN
        !  Rain(j,i)=Rain0(j,i)*g_TimeStep
        !ELSE
        !  Rain(j,i)=0.0
        !END IF
        !IF(PET0(j,i)>0.0)THEN
        !  PET(j,i)=PET0(j,i)*g_TimeStep
        !ELSE
        !  PET(j,i)=0.0
        !END IF

        CALL CREST_PrecipInt(Rain(j,i),g_tParams%RainFact(j,i), Rain(j,i))
        CALL CREST_EPotential(PET(j,i),g_tParams%KE(j,i), EPot(j,i))

        CALL CREST_RunoffGen(W0(j,i), &
          Rain(j,i),EPot(j,i), &
          g_tParams%WM(j,i),g_tParams%IM(j,i), &
          g_tParams%B(j,i),g_tParams%Ksat(j,i)*g_TimeStep, &
          W(j,i),ExcS(j,i),ExcI(j,i))

        CALL CREST_EAct(W0(j,i),Rain(j,i),EPot(j,i),W(j,i),EAct(j,i))

        SS0(j,i)=SS0(j,i)+ExcS(j,i)
        RS(j,i)=SS0(j,i)*g_tParams%KS(j,i)
        SS0(j,i)=SS0(j,i)*(1.0-g_tParams%KS(j,i))

        SI0(j,i)=SI0(j,i)+ExcI(j,i)
        RI(j,i)=SI0(j,i)*g_tParams%KI(j,i)
        SI0(j,i)=SI0(j,i)*(1.0-g_tParams%KI(j,i))

        W0(j,i)=W(j,i)

        !if(j==g_tOutlet%Col .and. i==g_tOutlet%Row) print*,g_tParams%RainFact(j,i), &
        !  g_tParams%Ksat(j,i),g_tParams%WM(j,i),  g_tParams%B(j,i),   g_tParams%IM(j,i), &
        !  g_tParams%KE(j,i),  g_tParams%coeM(j,i),g_tParams%expM(j,i),g_tParams%coeR(j,i), &
        !  g_tParams%coeS(j,i),g_tParams%KS(j,i),  g_tParams%KI(j,i)


        !WRITE(*,*)"Test",i,j
        IF(TRIM(g_RoutingType)=="CLR") THEN
          !Runoff(j,i)=((RS(j,i)+RI(j,i))/g_TimeStep)*g_GridArea(j,i)*QUnitFactor(g_TimeMark) !m3/s Original
          RS=1.0
          RI=1.0
          ! Revised by Ke Zhang, July, 2014
          !Runoff(j,i)=Runoff(j,i)+((RS(j,i)+RI(j,i))/g_TimeStep)*g_GridArea(j,i)*QUnitFactor(g_TimeMark)
          Runoff(j,i)=RS(j,i)+RI(j,i) !((RS(j,i)+RI(j,i)))*g_GridArea(j,i) ! Volume *QUnitFactor(g_TimeMark)/g_TimeStep
          !
          Roff_s_tmp=RS(j,i)*g_GridArea(j,i)! To Volume *QUnitFactor(g_TimeMark)/g_TimeStep !m3/s
          Roff_i_tmp=RI(j,i)*g_GridArea(j,i)! To Volume *QUnitFactor(g_TimeMark)/g_TimeStep !m3/s

          ! Revised by Xianwu Xue, December, 2014
          CALL ContinousLR(i,j, &
            oneRowA,oneColA,OnePerB, &
            g_NextR,g_NextC, &
            g_Mask, g_GridArea, &
            Runoff, g_NRows, &
            g_NCols,g_NoData_Value, &
            Roff_s_tmp)
          CALL ContinousLR(i,j, &
            twoRowA,twoColA,TwoPerB, &
            g_NextR,g_NextC, &
            g_Mask, g_GridArea, &
            Runoff, g_NRows, &
            g_NCols, g_NoData_Value, &
            Roff_i_tmp)

          ! Commented by Xianwu Xue, December, 2014
          !IF (oneRowA(j,i) /= g_NoData_Value .AND. oneColA(j,i) /= g_NoData_Value) THEN
          ! iii = g_NextR(j,i)
          ! jjj = g_NextC(j,i)
          ! DO WHILE ((InBasin(jjj,iii) .EQV. .TRUE.) &
          ! .AND. iii /= oneRowA(j,i) &
          ! .AND. jjj /= oneColA(j,i))
          ! Runoff(jjj,iii) = Runoff(jjj,iii)+Roff_s_tmp/g_GridArea(jjj,iii)
          ! iii2=g_NextR(jjj,iii)
          ! jjj2=g_NextC(jjj,iii)
          !
          ! iii=iii2
          ! jjj=jjj2
          !
          ! END DO
          !
          ! Runoff(oneColA(j,i), oneRowA(j,i))=Runoff(oneColA(j,i), oneRowA(j,i))+Roff_s_tmp*OnePerB(j,i)
          !
          !END IF
          !
          !IF (twoRowA(j,i) /= g_NoData_Value .AND. twoColA(j,i) /= g_NoData_Value) THEN
          ! iii = g_NextR(j,i)
          ! jjj = g_NextC(j,i)
          ! DO WHILE ((InBasin(jjj,iii).EQV. .TRUE.) &
          ! .AND. iii /= twoRowA(j,i) &
          ! .AND. jjj /= twoColA(j,i))
          ! Runoff(jjj,iii) = Runoff(jjj,iii)+Roff_i_tmp/g_GridArea(jjj,iii)
          ! iii2=g_NextR(jjj,iii)
          ! jjj2=g_NextC(jjj,iii)
          !
          ! iii=iii2
          ! jjj=jjj2
          !
          ! END DO
          !
          ! Runoff(twoColA(j,i), twoRowA(j,i))=Runoff(twoColA(j,i), twoRowA(j,i))+Roff_i_tmp*TwoPerB(j,i)
          !
          !END IF

        ELSE !JLR
          Runoff(j,i)=((RS(j,i)+RI(j,i))/g_TimeStep)*g_GridArea(j,i)*QUnitFactor(g_TimeMark) !m3/s Original
          !if(j==g_tOutlet%Col .and. i==g_tOutlet%Row) print*,Runoff(j,i),RS(j,i),RI(j,i)

        END If

      END DO    ! end column loop
    END DO    ! end row loop

    ! Revised by Xianwu Xue, December, 2014
    IF(TRIM(g_RoutingType)=="CLR") THEN
      !########################
      DO i=0, g_NRows-1
        DO j=0, g_NCols-1
          IF(g_Mask(j,i)/=1)THEN
            CYCLE
          END IF
          Runoff(j,i) = Runoff(j,i)*g_GridArea(j,i)*QUnitFactor(g_TimeMark)/g_TimeStep

        END DO
      END DO
    END IF

!--------------------------------------------------
! DownStream Routing
!--------------------------------------------------
    DO i=0, g_NRows-1
      DO j=0, g_NCols-1
        IF(g_Mask(j,i)==g_NoData_Value)THEN
          CYCLE
        END IF

        ii=OneRowA(j,i)
        jj=OneColA(j,i)
        IF(InBasin(jj,ii) .EQV. .TRUE.)THEN
          SS0(jj,ii)=SS0(jj,ii)+RS(j,i)*OnePerA(j,i) &
            *g_GridArea(j,i)/g_GridArea(jj,ii)
        END IF

        ii=OneRowB(j,i)
        jj=OneColB(j,i)
        IF(InBasin(jj,ii) .EQV. .TRUE.)THEN
          SS0(jj,ii)=SS0(jj,ii)+RS(j,i)*OnePerB(j,i) &
            *g_GridArea(j,i)/g_GridArea(jj,ii)
        END IF

        ii=TwoRowA(j,i)
        jj=TwoColA(j,i)
        IF(InBasin(jj,ii) .EQV. .TRUE.)THEN
          SI0(jj,ii)=SI0(jj,ii)+RI(j,i)*TwoPerA(j,i) &
            *g_GridArea(j,i)/g_GridArea(jj,ii)
        END IF

        ii=TwoRowB(j,i)
        jj=TwoColB(j,i)
        IF(InBasin(jj,ii) .EQV. .TRUE.)THEN
          SI0(jj,ii)=SI0(jj,ii)+RI(j,i)*TwoPerB(j,i) &
            *g_GridArea(j,i)/g_GridArea(jj,ii)
        END IF

      END DO
    END DO

!--------------------------------------------------
! Data output for grids and specific dates
!--------------------------------------------------
    IF(g_RunStyle/=2)THEN

      !Output the state of outlet
      IF(g_tOutlet%HasOutlet)THEN
        CALL CalculateOutletData(k,Rain,PET,EPot,EAct, &
          W,Runoff,ExcS,ExcI,RS,RI)
      END IF

      !Output the state of the Pixes
      CALL CalculateOutPixData(k,Rain,PET,EPot,EAct, &
        W,Runoff,ExcS,ExcI,RS,RI)

      !Grid Variable Outputs g_QPFBaseDate

      CALL ExportGridVar(strDate,Rain,PET,EPot,EAct, &
        W,Runoff,ExcS,ExcI,RS,RI)

      !Output specific Date
      Do i=LBOUND(g_OutDTIndex,1), UBOUND(g_OutDTIndex,1)
        IF(g_OutDTIndex(i)==k)THEN

          CALL ExportOutDTData(strDate,Rain,PET,EPot,EAct, &
            W,Runoff,ExcS,ExcI,RS,RI)

        END IF
      END DO

    END IF

    IF(g_RunStyle==2)THEN
      ii=g_tCalibSta(g_RegNum)%Row
      jj=g_tCalibSta(g_RegNum)%Col
      g_tCalibSta(g_RegNum)%RSim(k)=Runoff(jj,ii)
      !print*,g_tCalibSta(g_RegNum)%RSim(k)
      !print*,Runoff(g_tOutlet%Col, g_tOutlet%Row)
    END IF

    CALL myDAdd(g_TimeMark,g_TimeStep,dtTemp)

  END DO !End time loop


!--------------------------------------------------
! Data output for pixes and outlets
!--------------------------------------------------
  IF(g_RunStyle/=2)THEN

    CALL ExportPixData()

    ! Save the State Data
    CALL SaveStates(strDate,W0,SS0,SI0)

    IF(g_tOutlet%HasOutlet)THEN
      CALL ExportOutletData()

      IF(g_tOutlet%HasOBS .EQV. .TRUE.)THEN
        CALL GetNSCE(g_tOutlet%RObs, &
          g_tOutlet%R,g_NumWarmup,g_ITMax, &
          NSCE)
        CALL GetBias(g_tOutlet%RObs, &
          g_tOutlet%R,g_NumWarmup,g_ITMax, &
          Bias)
        CALL GetCC(g_tOutlet%RObs, &
          g_tOutlet%R,g_NumWarmup,g_ITMax, &
          CC)

        WRITE(*,*)
        WRITE(*,*)"  The results of the Outlet is: "
        WRITE(*,*)"    NSCE: ",NSCE
        WRITE(*,*)"    Bias(%): ",Bias
        WRITE(*,*)"    CC: ",CC

        WRITE(g_CREST_LogFileID,*)
        WRITE(g_CREST_LogFileID,*)"  The results of the Outlet is: "
        WRITE(g_CREST_LogFileID,*)"     NSCE: ",NSCE
        WRITE(g_CREST_LogFileID,*)"     Bias(%): ",Bias
        WRITE(g_CREST_LogFileID,*)"     CC: ",CC
      END IF
    END IF

  END IF

  ! DEALLOCATE
  DEALLOCATE(OneRowA)
  DEALLOCATE(OneRowB)
  DEALLOCATE(OneColA)
  DEALLOCATE(OneColB)

  DEALLOCATE(TwoRowA)
  DEALLOCATE(TwoRowB)
  DEALLOCATE(TwoColA)
  DEALLOCATE(TwoColB)

  DEALLOCATE(OnePerA)
  DEALLOCATE(OnePerB)
  DEALLOCATE(TwoPerA)
  DEALLOCATE(TwoPerB)

  RETURN

  CONTAINS

!============================================================
!
!============================================================
  SUBROUTINE ContinousLR(i,j, ToRowA,ToColA,ToPerB, NextR,NextC, Mask, GridArea, Runoff, NRows,NCols, NoData_Value,Roff_tmp)
  IMPLICIT NONE

  INTEGER, INTENT (IN), ALLOCATABLE :: ToRowA(:,:), ToColA(:,:)
  INTEGER, INTENT (IN), ALLOCATABLE :: NextR(:,:), NextC(:,:), Mask(:,:)
  REAL (KIND=8), INTENT (IN), ALLOCATABLE :: GridArea(:,:), ToPerB(:,:)
  REAL (KIND=8), INTENT (INOUT), ALLOCATABLE :: Runoff(:,:)


  INTEGER (KIND=4) :: NRows, NCols
  REAL (KIND=8) :: NoData_Value
  INTEGER (KIND=4),INTENT (IN) :: i,j
  INTEGER (KIND=4) :: k,ii,jj

  REAL (KIND=8),INTENT (IN) :: Roff_tmp

  INTEGER (KIND=4) :: iii,jjj,iii2,jjj2


  IF (ToRowA(j,i) /= INT(NoData_Value) .AND. ToColA(j,i) /= INT(NoData_Value)) THEN
    iii = NextR(j,i)
    jjj = NextC(j,i)
    DO WHILE (InBasin_4(iii,jjj,NRows, NCols) &
      .AND. iii /= ToRowA(j,i) &
      .AND. jjj /= ToColA(j,i))
    !WRITE(*,*)"test-1"
    Runoff(jjj,iii) = Runoff(jjj,iii)+Roff_tmp/GridArea(jjj,iii)
    iii2=NextR(jjj,iii)
    jjj2=NextC(jjj,iii)

    iii=iii2
    jjj=jjj2

    END DO
    !WRITE(*,*)"test-1"
    iii=ToRowA(j,i)
    jjj=ToColA(j,i)
    Runoff(jjj,iii) &
      =Runoff(jjj, iii)+(Roff_tmp/GridArea(jjj,iii))*ToPerB(j,i) !
  ELSE
    ! For the Downstream cells which will flow out of the basin
    ! It will
    !WRITE(*,*)"test-2",ToRowA(j,i)
    iii = NextR(j,i)
    jjj = NextC(j,i)
    DO WHILE (InBasin_4(iii,jjj,NRows, NCols) .AND. Mask(jjj,iii)==1)
      !WRITE(*,*)"test-3",pBasics%Mask(jjj,iii)
      Runoff(jjj,iii) = Runoff(jjj,iii)+Roff_tmp/GridArea(jjj,iii)
      iii2=NextR(jjj,iii)
      jjj2=NextC(jjj,iii)

      iii=iii2
      jjj=jjj2

    END DO

  END IF
  Return
  END SUBROUTINE ContinousLR

!============================================================
!
!============================================================
  FUNCTION InBasin_4(iRow, jCol,NRows, NCols) RESULT(IsInBasin)
  IMPLICIT NONE
  INTEGER (KIND=4) :: iRow,jCol,NRows, NCols
  LOGICAL :: IsInBasin
  IF(jCol < 0 .OR. jCol >= NCols .OR. iRow < 0 .OR. iRow >= NRows)THEN
    IsInBasin=.FALSE.
  ELSE
    ! IF(g_Mask(jCol,iRow)==g_NoData_Value)THEN
    ! InBasin=.FALSE.
    ! ELSE
    IsInBasin=.TRUE.
    ! END IF
  END IF
  RETURN
  END FUNCTION InBasin_4
  END SUBROUTINE CREST_Simu

!============================================================
!
!============================================================
  SUBROUTINE CREST_EAct(W0,P,EPot,W,EAct)

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(IN) :: W0,P,EPot,W
  DOUBLE PRECISION,INTENT(OUT) :: EAct


  IF(P>EPot)THEN
    EAct=EPot
  ELSE !P<=EPot
    EAct=W0-W
  END IF

  RETURN
  END SUBROUTINE CREST_EAct



!============================================================
!
!============================================================
  SUBROUTINE CREST_EPotential(E,EFact,EPot)

  IMPLICIT NONE

  DOUBLE PRECISION :: E,EFact
  DOUBLE PRECISION :: EPot

  EPot=E*EFact

  RETURN
  END SUBROUTINE CREST_EPotential

!============================================================
!
!============================================================
  SUBROUTINE CREST_PrecipInt(P,PFact,PInt)

  IMPLICIT NONE

  DOUBLE PRECISION :: P,PFact
  DOUBLE PRECISION :: PInt

  PInt=P*PFact

  RETURN
  END SUBROUTINE CREST_PrecipInt


!============================================================
!
!============================================================
  SUBROUTINE CREST_ReturnPeriod()
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE

  INTEGER :: i,j,k,l
  INTEGER :: lYmin, lYmax
  InteGEr :: D(1:6),DEnd(1:6)
  DOUBLE PRECISION,ALLOCATABLE:: Q(:,:,:), V(:,:)
  Character(80):: FileName, strDate
  Logical:: bIsError
  DOUBLE PRECISION :: temV

  If (g_StartDate(2)<=4) THEN
    lYmin=g_StartDate(1)
  ELSE
    lYmin=g_StartDate(1)+1
  END IF

  If (g_EndDate(2)>=9) THEN
    lYmax=g_EndDate(1)
  ELSE
    lYmax=g_EndDate(1)-1
    IF(g_EndDate(1)-1<g_StartDate(1))THEN
      lYmax=lYmin
    END IF
  END IF

  ALLOCATE(Q(lYmin:lYmax, 0:g_NCols-1, 0:g_NRows-1))
  ALLOCATE(V(0:g_NCols-1, 0:g_NRows-1))
  Q=g_NoData_Value

  DO k=lYmin, lYmax
    D(1)=k
    D(2)=1
    D(3)=1
    D(4)=0
    D(5)=0
    D(6)=0
    DEnd=D
    DEnd(1)=k+1

    DO WHILE (myDcompare2(D, DEnd,g_TimeMark)==1)
      strDate=myDtoStr(D, g_TimeMark)

      WRITE(*,*) "  Loading Runoff File: " // TRIM(strDate)
      WRITE(g_CREST_LogFileID,*) "  Loading Runoff File: " &
        // TRIM(strDate)

      FileName=TRIM(g_ResultPath)//"GOVar_"  &
        //TRIM("R") // "_" //TRIM(strDate)

      CALL ReadMatrixFile(FileName,V, &
        g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
        g_CellSize,g_NoData_Value,bIsError, &
        g_ResultFormat,"")
      IF(bIsError .EQV. .TRUE.)THEN
        CALL myDadd(g_TimeMark, g_TimeStep, D)
        CYCLE
      END IF

      DO i=0, g_NRows-1
        DO j=0, g_NCols-1
          IF (Q(k,j,i) < V(j,i))THEN
            Q(k,j,i) = V(j,i)
          END IF
        END DO
      END DO

      CALL myDadd(g_TimeMark, g_TimeStep, D)
    END DO
  END DO

  WRITE(*,*) " Sorting for Return Period!"
  WRITE(g_CREST_LogFileID,*) " Sorting for Return Period!"

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      DO k=lYmin, lYmax-1
        DO l=k+1, lYmax
          IF (Q(k, j, i)<Q(l, j, i)) THEN
            temV=Q(k, j, i)
            Q(k, j, i)=Q(l, j, i)
            Q(l, j, i)=temV
          END IF
        END DO
      END DO
    END DO
  END DO

  DO k=lYmin, lYmax
    l=lYmax-k+1
    WRITE (*,*) "Exporting Num", l
    WRITE (g_CREST_LogFileID,*) "Exporting Num", l

    DO i=0, g_NRows-1
      DO j=0, g_NCols-1
        V(j,i)=Q(k, j, i)
      END DO
    END DO
    FileName=TRIM(g_ResultPath) // "Level."   &
      // TRIM(PosIntMsg(l))

    CALL WriteMatrixFile(FileName, &
      V, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, &
      g_NoData_Value, &
      bIsError,g_ResultFormat)

  END DO

  RETURN
  END SUBROUTINE CREST_ReturnPeriod


!============================================================
!
!============================================================
  SUBROUTINE CREST_RunoffRoutePre( &
    OneRowA,OneColA,OnePerA, &
    OneRowB,OneColB,OnePerB, &
    TwoRowA,TwoColA,TwoPerA, &
    TwoRowB,TwoColB,TwoPerB)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  LOGICAL ::bIsNotExist
  DOUBLE PRECISION, ALLOCATABLE :: NextTime(:,:)   ! velocity (m/s)

  INTEGER :: OneRowA(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: OneRowB(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: OneColA(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: OneColB(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: TwoRowA(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: TwoRowB(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: TwoColA(0:g_NCols-1,0:g_NRows-1)
  INTEGER :: TwoColB(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: OnePerA(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: OnePerB(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: TwoPerA(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: TwoPerB(0:g_NCols-1,0:g_NRows-1)
  LOGICAL :: bIsError


  ALLOCATE(NextTime(0:g_NCols-1,0:g_NRows-1))

  CALL CalNextTime(NextTime)

  IF(IsDebug==1)THEN
    !XianwuXue 2014-12-22
    CALL WriteMatrixFile(TRIM(g_ResultPath) // "NextTime",NextTime,  &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
  END IF

  CALL CREST_RouteTreat(NextTime,OneRowA,OneColA,OnePerA, &
    OneRowB,OneColB,OnePerB)

  IF(IsDebug==1)THEN
    !XianwuXue 2014-12-22
    CALL WriteMatrixFile_Int(TRIM(g_ResultPath) // "ToRowA",OneRowA,  &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
    CALL WriteMatrixFile_Int(TRIM(g_ResultPath) // "ColA",OneColA,  &
      g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
  END IF

  WHERE(NextTime/=g_NoData_Value)
    NextTime=NextTime/g_tParams%coeS
  END WHERE

  CALL CREST_RouteTreat(NextTime,TwoRowA,TwoColA,TwoPerA, &
    TwoRowB,TwoColB,TwoPerB)

  !Deallocate the arrays
  DEALLOCATE(NextTime)

  RETURN
  END SUBROUTINE CREST_RunoffRoutePre


!============================================================
!
!============================================================
  SUBROUTINE CREST_RouteTreat(NextTime,toRowA,toColA,toPerA, &
    toRowB,toColB,toPerB)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE

  DOUBLE PRECISION :: NextTime(0:g_NCols-1,0:g_NRows-1)

  INTEGER, INTENT(OUT)::toRowA(0:g_NCols-1,0:g_NRows-1)
  INTEGER, INTENT(OUT)::toRowB(0:g_NCols-1,0:g_NRows-1)
  INTEGER, INTENT(OUT)::toColA(0:g_NCols-1,0:g_NRows-1)
  INTEGER, INTENT(OUT)::toColB(0:g_NCols-1,0:g_NRows-1)

  DOUBLE PRECISION, INTENT(OUT)::toPerA(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION, INTENT(OUT)::toPerB(0:g_NCols-1,0:g_NRows-1)

  INTEGER :: i,j


  toRowA=g_NoData_Value
  toRowB=g_NoData_Value
  toColA=g_NoData_Value
  toColB=g_NoData_Value

  toPerA=g_NoData_Value
  toPerB=g_NoData_Value

  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      IF(g_Mask(j,i)/=g_NoData_Value)THEN
        toRowB(j,i)=i
        toColB(j,i)=j
        toPerB(j,i)=0
        DO WHILE(toPerB(j,i)<g_TimeStep)
          toRowA(j,i)=toRowB(j,i)
          toColA(j,i)=toColB(j,i)
          toPerA(j,i)=toPerB(j,i)
          IF(toRowA(j,i)/=g_NoData_Value)THEN
            IF(InBasin(toColA(j,i),toRowA(j,i)))THEN
              toRowB(j,i) &
                =g_NextR(toColA(j,i),toRowA(j,i))
              toColB(j,i) &
                =g_NextC(toColA(j,i),toRowA(j,i))
              toPerB(j,i)=toPerB(j,i) &
                +NextTime(toColA(j,i),toRowA(j,i))
            ELSE
              toRowB(j,i)=g_NoData_Value
              toColB(j,i)=g_NoData_Value
              toPerB(j,i)=toPerB(j,i)+9999.0
            END IF

          ELSE
            toRowB(j,i)=g_NoData_Value
            toColB(j,i)=g_NoData_Value
            toPerB(j,i)=toPerB(j,i)+9999.0
          END IF
        END DO
        toPerB(j,i)=(g_TimeStep-toPerA(j,i)) &
          /(toPerB(j,i)-toPerA(j,i))
        toPerA(j,i)=1.0-toPerB(j,i)
      END IF
    END DO
  END DO


  RETURN
  END SUBROUTINE CREST_RouteTreat

!============================================================
!
!============================================================
  SUBROUTINE CalNextTime(NextTime)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param

  IMPLICIT NONE
  DOUBLE PRECISION :: SpeedVegLocal,SpeedVegNext

  INTEGER :: i,j,ii,jj
  DOUBLE PRECISION, ALLOCATABLE :: Speed(:,:)
  DOUBLE PRECISION :: NextTime(0:g_NCols-1,0:g_NRows-1)
  ALLOCATE(Speed(0:g_NCols-1,0:g_NRows-1))
  DO i=0, g_NRows-1
    DO j=0, g_NCols-1
      ii=g_NextR(j,i)
      jj=g_NextC(j,i)
      IF(g_Mask(j,i)/=g_NoData_Value)THEN
        SpeedVegLocal=0.5
        IF(InBasin(jj,ii) .EQV. .FALSE.)THEN
          SpeedVegNext=SpeedVegLocal
        ELSE
          SpeedVegNext=0.5
        END IF
        Speed(j,i)=g_tParams%coeM(j,i) &
          *((SpeedVegLocal+SpeedVegNext)/2.0) &
          *g_Slope(j,i)**g_tParams%expM(j,i)
        IF(g_Stream(j,i)==1)THEN !g_FAC(j,i)>g_tParams%TH(j,i))THEN
          Speed(j,i)=Speed(j,i)*g_tParams%coeR(j,i)
        END IF
        NextTime(j,i)=g_NextLen(j,i)/Speed(j,i)*QUnitFactor(g_TimeMark)/1000.0 !Unit=meter/second
      ELSE
        Speed(j,i)=g_NoData_Value
        NextTime(j,i)=g_NoData_Value
      END IF
    END DO
  END DO

  DEALLOCATE(Speed)

  RETURN
  END SUBROUTINE CalNextTime

!============================================================
!
!============================================================
  SUBROUTINE CREST_RunoffGen(W0,P,EPot,WM,IM,B,Ksat,W,ExcS,ExcI)

  IMPLICIT NONE

  DOUBLE PRECISION,INTENT(IN) :: W0,P,EPot,WM,IM,B,Ksat
  DOUBLE PRECISION,INTENT(OUT) :: W,ExcS,ExcI

  DOUBLE PRECISION :: PSoil,R,temX,WMM,A

  IF(P>EPot)THEN
    ! Calculate part of precip that goes into soil
    PSoil=(P-EPot)*(1.0-IM)

    IF(W0<WM)THEN
      WMM=WM*(1.0+B)
      A=WMM*(1.0-(1.0-W0/WM)**(1.0/(1.0+B)))
      IF(PSoil+A>=WMM)THEN
        R=PSoil-(WM-W0)
        W=WM
      ELSE
!-ygan
        !R=PSoil-WM*((1.0-A/WMM)**(1.0+B) &
        !  -(1.0-(A+PSoil)/WMM)**(1.0+B))
!.ygan
!+ygan  (same with the above equation)
        R=PSoil-(WM-W0)+WM*(1.0-(A+PSoil)/WMM)**(1.0+B)
!.ygan
        IF(R<0.0)THEN
          R=0.0
        END IF
        W=W0+PSoil-R
      END IF
    ELSE ! W0>WM the soil is full
      R=PSoil
      W=W0
    END IF
    ! Calculate how much water can infiltrate
    temX=((W0+W)/2.0)*KSat/WM
    IF(R<=temX)THEN
      ExcI=R
    ELSE
      ExcI=temX
    END IF
    ExcS=R-ExcI+(P-EPot)*IM
    IF(ExcS<0)THEN !Added by Dr.Xianwu Xue 2011.5.13
      ExcS=0.0
    END IF
  ELSE !P<=EPot
    ExcS=0.0
    ExcI=0.0
    temX=(EPot-P)*W0/WM
    IF(temX<W0)THEN
      W=W0-temX
    ELSE
      W=0.0
    END IF
  END IF

  RETURN
  END SUBROUTINE CREST_RunoffGen

!============================================================
!
!============================================================
  SUBROUTINE CREST_RealTime()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_ICS
  USE CREST_Param
  USE CREST_Calib_SCEUA

  IMPLICIT NONE

  LOGICAL :: bIsError
  CHARACTER(14):: strDate

  INTEGER :: dtTemp(1:6),dtNow(1:6),dtStart(1:6),dtEnd(1:6)
  INTEGER :: k,i

  CALL myNow(dtNow)
  dtTemp=dtNow

  SELECT CASE (g_TimeMark)
  CASE ("y","Y")
    i=1
  CASE ("m","M")
    i=2
  CASE ("d","D")
    i=3
  CASE ("h","H")
    i=4
  CASE ("u","U")
    i=5
  CASE ("s","S")
    i=6
  END SELECT
  DO k=i,UBOUND(dtTemp,1)
    dtTemp(k)=0
  END DO
  DO WHILE(dtTemp(i)<=dtNow(i))
    dtTemp(i)=dtTemp(i)+g_TimeStep
  END DO
  dtTemp(i)=dtTemp(i)-g_TimeStep

  dtEnd=g_EndDate
  IF(myDCompare(dtTemp,g_EndDate,g_TimeMark)==-1)THEN !indicate dtTemp>g_EndDate
    g_EndDate=dtTemp
  END IF

  !Get the end date of Rain file
  bIsError=.TRUE.
  DO WHILE(bIsError)

    strDate=myDtoStr(g_EndDate, g_TimeMark)
    CALL InquireMatrixFile(TRIM(g_RainPath),bIsError, &
      g_RainFormat,strDate)
    IF(myDEqual(g_StartDate,g_EndDate))THEN
      bIsError=.FALSE. !Indicate the file exist!
    ELSE
      IF(bIsError)THEN
        CALL myDadd(g_TimeMark,-g_TimeStep,g_EndDate)
      END IF
    END IF
  END DO

  !Get date of last status
  dtStart=g_StartDate
  g_StartDate=g_EndDate
  CALL myDadd(g_TimeMark,g_TimeStep,g_StartDate)
  bIsError=.TRUE.
  DO WHILE(bIsError)
    strDate=myDtoStr(g_StartDate, g_TimeMark)
    CALL InquireMatrixFile(TRIM(g_StatePath) &
      // "State_"//TRIM(strDate)// "_W0",  &
      bIsError,g_StateFormat,"")

    IF(myDEqual(g_StartDate,dtStart))THEN
      bIsError=.FALSE. !Indicate the file exist!
    ELSE
      IF(bIsError)THEN
        CALL myDadd(g_TimeMark,-g_TimeStep,g_StartDate)
      END IF
    END IF
  END DO

  CALL myDadd(g_TimeMark,g_TimeStep,g_StartDate) ! Modified by Xianwu Xue 2011.3.24

  IF(myDCompare(g_StartDate,g_EndDate,g_TimeMark)==1)THEN !When g_StartDate<g_EndDate

    WRITE(*,*)"  Model Start Time: " &
      //myDtoStr(g_StartDate, g_TimeMark)
    WRITE(*,*)"  Model End Time: " &
      //myDtoStr(g_EndDate, g_TimeMark)
    WRITE(g_CREST_LogFileID,*)"  Model Start Time: " &
      //myDtoStr(g_StartDate, g_TimeMark)
    WRITE(g_CREST_LogFileID,*)"  Model End Time: " &
      //myDtoStr(g_EndDate, g_TimeMark)
    WRITE(*,*)
    WRITE(g_CREST_LogFileID,*)

    g_LoadState=.TRUE.
    g_SaveState=.TRUE.

    g_WarmupDate=g_StartDate
    g_NumWarmup=0
    g_ITMax=myDdiff(g_TimeMark, g_StartDate, g_EndDate) &
      / g_TimeStep + 1
    IF(g_tOutlet%HasOBS .EQV. .TRUE.)THEN
      DEALLOCATE(g_tOutlet%RObs)
      ALLOCATE(g_tOutlet%RObs(0:g_ITMax-1))
      CALL XXWReadRunoffObs(TRIM(g_OBSPath) &
        //TRIM(g_tOutlet%Name)//"_Obs.csv",  &
        g_tOutlet%RObs,bIsError)
    END IF
    ! IF(bIsError .EQV. .TRUE.)THEN
    ! g_tOutlet%HasOBS=.FALSE.
    ! ELSE
    ! g_tOutlet%HasOBS=.TRUE.
    ! END IF

    Do i=LBOUND(g_tOutPix,1), UBOUND(g_tOutPix,1)
      DEALLOCATE(g_tOutPix(i)%RObs)
      ALLOCATE(g_tOutPix(i)%RObs(0:g_ITMax-1))
      CALL XXWReadRunoffObs(TRIM(g_OBSPath) &
        //TRIM(g_tOutPix(i)%Name)//"_Obs.csv",  &
        g_tOutPix(i)%RObs,bIsError)
      IF(bIsError .EQV. .TRUE.)THEN
        g_tOutPix(i)%HasOBS=.FALSE.
        bIsError=.FALSE.
      ELSE
        g_tOutPix(i)%HasOBS=.TRUE.
        bIsError=.FALSE.
      END IF
    END DO

    CALL CREST_Simu()

  END IF

  RETURN
  END SUBROUTINE CREST_RealTime

!============================================================
!
!============================================================
! Revised: Yanjun Gan, 06/01/2015

  SUBROUTINE CREST_Simu_For_Calib(x)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_ICS
  USE CREST_Param
  USE CREST_Calib_SCEUA
  IMPLICIT NONE
  DOUBLE PRECISION :: x(1:g_nopt)
  INTEGER :: j,k

  k=0
  DO j=1,UBOUND(g_sParamName,1)

    IF(g_tCalibSta(g_RegNum)%bu(j) &
      -g_tCalibSta(g_RegNum)%bl(j)>0)THEN
      k=k+1
      g_tCalibSta(g_RegNum)%x(j)=x(k)
      !print*,g_tCalibSta(g_RegNum)%x(j)
    END IF
  END DO

!-ygan
!  WHERE(g_RegMask==g_tCalibSta(g_RegNum)%Value)
!
!    g_tParamsAdj_Cali%RainFact = g_tCalibSta(g_RegNum)%x(1)
!    g_tParamsAdj_Cali%Ksat     = g_tCalibSta(g_RegNum)%x(2)
!    g_tParamsAdj_Cali%WM       = g_tCalibSta(g_RegNum)%x(3)
!    g_tParamsAdj_Cali%B        = g_tCalibSta(g_RegNum)%x(4)
!    g_tParamsAdj_Cali%IM       = g_tCalibSta(g_RegNum)%x(5)
!    g_tParamsAdj_Cali%KE       = g_tCalibSta(g_RegNum)%x(6)
!    g_tParamsAdj_Cali%coeM     = g_tCalibSta(g_RegNum)%x(7)
!
!    g_tParamsAdj_Cali%expM     = g_tCalibSta(g_RegNum)%x(8)
!    g_tParamsAdj_Cali%coeR     = g_tCalibSta(g_RegNum)%x(9)
!    g_tParamsAdj_Cali%coeS     = g_tCalibSta(g_RegNum)%x(10)
!    g_tParamsAdj_Cali%KS       = g_tCalibSta(g_RegNum)%x(11)
!    g_tParamsAdj_Cali%KI       = g_tCalibSta(g_RegNum)%x(12)
!
!  END WHERE
!
!  WHERE(g_Mask/=g_NoData_Value)
!    g_tParams%RainFact &
!      = g_tParams_Cali%RainFact * g_tParamsAdj_Cali%RainFact
!    g_tParams%Ksat &
!      = g_tParams_Cali%Ksat * g_tParamsAdj_Cali%Ksat
!    g_tParams%WM &
!      = g_tParams_Cali%WM * g_tParamsAdj_Cali%WM
!    g_tParams%B &
!      = g_tParams_Cali%B * g_tParamsAdj_Cali%B
!    g_tParams%IM &
!      = g_tParams_Cali%IM * g_tParamsAdj_Cali%IM
!    g_tParams%KE &
!      = g_tParams_Cali%KE * g_tParamsAdj_Cali%KE
!    g_tParams%coeM &
!      = g_tParams_Cali%coeM * g_tParamsAdj_Cali%coeM
!
!    g_tParams%expM &
!      = g_tParams_Cali%expM * g_tParamsAdj_Cali%expM
!    g_tParams%coeR &
!      = g_tParams_Cali%coeR * g_tParamsAdj_Cali%coeR
!    g_tParams%coeS &
!      = g_tParams_Cali%coeS * g_tParamsAdj_Cali%coeS
!    g_tParams%KS &
!      = g_tParams_Cali%KS * g_tParamsAdj_Cali%KS
!    g_tParams%KI &
!      = g_tParams_Cali%KI * g_tParamsAdj_Cali%KI
!
!  END WHERE
!.ygan

!+ygan
  WHERE(g_RegMask==g_tCalibSta(g_RegNum)%Value)

    g_tParams%RainFact = g_tCalibSta(g_RegNum)%x(1)
    g_tParams%Ksat     = g_tCalibSta(g_RegNum)%x(2)
    g_tParams%WM       = g_tCalibSta(g_RegNum)%x(3)
    g_tParams%B        = g_tCalibSta(g_RegNum)%x(4)
    g_tParams%IM       = g_tCalibSta(g_RegNum)%x(5)
    g_tParams%KE       = g_tCalibSta(g_RegNum)%x(6)
    g_tParams%coeM     = g_tCalibSta(g_RegNum)%x(7)

    g_tParams%expM     = g_tCalibSta(g_RegNum)%x(8)
    g_tParams%coeR     = g_tCalibSta(g_RegNum)%x(9)
    g_tParams%coeS     = g_tCalibSta(g_RegNum)%x(10)
    g_tParams%KS       = g_tCalibSta(g_RegNum)%x(11)
    g_tParams%KI       = g_tCalibSta(g_RegNum)%x(12)

  END WHERE
!.ygan

  !k=g_tCalibSta(g_RegNum)%Row
  !j=g_tCalibSta(g_RegNum)%Col
  !print*,g_tParams%RainFact(j,k), g_tParams%Ksat(j,k), g_tParams%WM(j,k)  , g_tParams%B(j,k)   , &
  !       g_tParams%IM(j,k)      , g_tParams%KE(j,k)  , g_tParams%coeM(j,k), g_tParams%expM(j,k), &
  !       g_tParams%coeR(j,k)    , g_tParams%coeS(j,k), g_tParams%KS(j,k)  , g_tParams%KI(j,k)

  CALL CREST_Simu()

  RETURN
  END SUBROUTINE CREST_Simu_For_Calib

!============================================================
!
!============================================================
  SUBROUTINE CREST_OutDISParamsToFile(sWS,sPrefix)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  LOGICAL :: bIsError
  CHARACTER(*) sWS,sPrefix

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"RainFact",  &
    g_tParams%RainFact, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"Ksat", &
    g_tParams%Ksat, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"WM",   &
    g_tParams%WM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"B",   &
    g_tParams%B, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"IM",   &
    g_tParams%IM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KE",  &
    g_tParams%KE, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeM",  &
    g_tParams%coeM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"expM",   &
    g_tParams%expM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeR",  &
    g_tParams%coeR, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeS",  &
    g_tParams%coeS, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KS",  &
    g_tParams%KS, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  CALL WriteMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KI",  &
    g_tParams%KI, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat)

  END SUBROUTINE CREST_OutDISParamsToFile

!============================================================
!
!============================================================
  SUBROUTINE CREST_ReadDISParamsFile(sWS,sPrefix,tParams,bIsError)
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  LOGICAL :: bIsError
  CHARACTER(*) sWS,sPrefix
  TYPE(CREST_Params):: tParams

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"RainFact",   &
    tParams%RainFact, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"Ksat",  &
    tParams%Ksat, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"WM",  &
    tParams%WM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"B",   &
    tParams%B, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"IM",  &
    tParams%IM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KE",  &
    tParams%KE, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeM",  &
    tParams%coeM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"expM",  &
    tParams%expM, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeR",  &
    tParams%coeR, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"coeS",  &
    tParams%coeS, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KS",   &
    tParams%KS, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(sWS)//TRIM(sPrefix)//"KI",  &
    tParams%KI, g_NCols,g_NRows, &
    g_XLLCorner,g_YLLCorner, &
    g_CellSize, g_NoData_Value,bIsError,g_ParamFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN
    RETURN
  END IF

  END SUBROUTINE CREST_ReadDISParamsFile


!============================================================
!
!============================================================
  SUBROUTINE CREST_SCEUA()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_ICS
  USE CREST_Param
  USE CREST_Calib_SCEUA
  IMPLICIT NONE
  INTEGER :: i,k,j,iOutputID

  DOUBLE PRECISION :: dblValue
  LOGICAL :: bIsError
  DOUBLE PRECISION, ALLOCATABLE :: a(:)
  DOUBLE PRECISION, ALLOCATABLE :: bl(:)
  DOUBLE PRECISION, ALLOCATABLE :: bu(:)
!-ygan
!  CHARACTER(LEN=80) :: strNum
!_ygan

  dblValue=1.0

  CALL InitParamsType(g_tParams_Cali)
  CALL InitParamsType(g_tParamsAdj_Cali)
  CALL AssignValueToParamsType(g_tParamsAdj_Cali,dblValue)
  CALL EqualParamsType(g_tParams,g_tParams_Cali)

  CALL XXWGetFreeFile(iOutputID)
  OPEN(UNIT=iOutputID,FILE=TRIM(g_ResultPath) &
    //TRIM("SCEUAOut_")//TRIM(g_strD)// '.dat',STATUS='UNKNOWN')

! Optimize parameters station by station
  DO i=0,g_NCalibStas-1

    g_RegNum=i

!-ygan
!    WRITE(strNum,*)g_tCalibSta(i)%Value
!.ygan

    g_nopt=COUNT(g_tCalibSta(i)%bu-g_tCalibSta(i)%bl>0)
    ALLOCATE(a(1:g_nopt))
    ALLOCATE(bl(1:g_nopt))
    ALLOCATE(bu(1:g_nopt))

    k=0
    DO j=1,UBOUND(g_sParamName,1)

      IF(g_tCalibSta(i)%bu(j)-g_tCalibSta(i)%bl(j)>0)THEN
        k=k+1
        a(k)=g_tCalibSta(i)%a(j)
        bl(k)=g_tCalibSta(i)%bl(j)
        bu(k)=g_tCalibSta(i)%bu(j)
        !print*,a(k),bl(k),bu(k)
      END IF
    END DO

    CALL GeneralSCEUA2(a,bl,bu, &
      g_nopt,g_maxn, &
      g_kstop,g_pcento,g_iseed,g_ngs,55, TRIM(g_ResultPath),iOutputID)

    DEALLOCATE(a)
    DEALLOCATE(bl)
    DEALLOCATE(bu)
  END DO
  CLOSE(iOutputID)

! Use the optimized parameters to run the model
  g_RunStyle=1
  CALL CREST_Simu()

! Write the optimizaed parameter values
  CALL CREST_OutDISParamsToFile(TRIM(g_ResultPath),"")

  RETURN

  END SUBROUTINE CREST_SCEUA

!============================================================
!
!============================================================
  DOUBLE PRECISION FUNCTION XXWCRESTFunctn(nopt,x)
  USE CREST_Project
  USE CREST_Calib_SCEUA
  IMPLICIT NONE
  INTEGER nopt
  DOUBLE PRECISION :: x(1:nopt)
  CHARACTER(LEN=20):: sPrefix

  CALL CREST_Simu_For_Calib(x)

  CALL GetNSCE(g_tCalibSta(g_RegNum)%RObs, &
    g_tCalibSta(g_RegNum)%RSim,g_NumWarmup,g_ITMax, &
    XXWCRESTFunctn)

  XXWCRESTFunctn=1.0-XXWCRESTFunctn

  !WRITE(*,"(4X,F20.15,3X,A15,I3)") 1.0-XXWCRESTFunctn, &
  ! "Region Number: ",g_tCalibSta(g_RegNum)%Value
  WRITE(*,*) " NSCE = ",1.0-XXWCRESTFunctn

  RETURN
  END FUNCTION XXWCRESTFunctn

!============================================================
!
!============================================================
  SUBROUTINE EqualParamsType(tParamsIn,tParamsOut)
  USE CREST_Param
  IMPLICIT NONE
  TYPE(CREST_Params) tParamsIn,tParamsOut

  tParamsOut%RainFact = tParamsIn%RainFact
  tParamsOut%Ksat     = tParamsIn%Ksat
  tParamsOut%WM       = tParamsIn%WM
  tParamsOut%B        = tParamsIn%B
  tParamsOut%IM       = tParamsIn%IM
  tParamsOut%KE       = tParamsIn%KE
  tParamsOut%coeM     = tParamsIn%coeM

  tParamsOut%expM     = tParamsIn%expM
  tParamsOut%coeR     = tParamsIn%coeR
  tParamsOut%coeS     = tParamsIn%coeS
  tParamsOut%KS       = tParamsIn%KS
  tParamsOut%KI       = tParamsIn%KI

  RETURN
  END SUBROUTINE EqualParamsType

!============================================================
!
!============================================================
  SUBROUTINE InitParamsType(tParams)
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  TYPE(CREST_Params) tParams

  ALLOCATE(tParams%RainFact(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%Ksat(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%WM(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%B(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%IM(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%KE(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%coeM(0:g_NCols-1,0:g_NRows-1))

  ALLOCATE(tParams%expM(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%coeR(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%coeS(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%KS(0:g_NCols-1,0:g_NRows-1))
  ALLOCATE(tParams%KI(0:g_NCols-1,0:g_NRows-1))

  tParams%RainFact=g_NoData_Value
  tParams%Ksat=g_NoData_Value
  tParams%WM=g_NoData_Value
  tParams%B=g_NoData_Value
  tParams%IM=g_NoData_Value
  tParams%KE=g_NoData_Value
  tParams%coeM=g_NoData_Value

  tParams%expM=g_NoData_Value
  tParams%coeR=g_NoData_Value
  tParams%coeS=g_NoData_Value
  tParams%KS=g_NoData_Value
  tParams%KI=g_NoData_Value

  RETURN
  END SUBROUTINE InitParamsType

!============================================================
!
!============================================================
  SUBROUTINE AssignValueToParamsType(tParams,dblValue)
  USE CREST_Basic
  USE CREST_Param
  IMPLICIT NONE
  TYPE(CREST_Params) tParams
  DOUBLE PRECISION:: dblValue

  WHERE(g_Mask==1)
    tParams%RainFact=dblValue
    tParams%Ksat=dblValue
    tParams%WM=dblValue
    tParams%B=dblValue
    tParams%IM=dblValue
    tParams%KE=dblValue
    tParams%coeM=dblValue

    tParams%expM=dblValue
    tParams%coeR=dblValue
    tParams%coeS=dblValue
    tParams%KS=dblValue
    tParams%KI=dblValue

  ELSEWHERE
    tParams%RainFact=g_NoData_Value
    tParams%Ksat=g_NoData_Value
    tParams%WM=g_NoData_Value
    tParams%B=g_NoData_Value
    tParams%IM=g_NoData_Value
    tParams%KE=g_NoData_Value
    tParams%coeM=g_NoData_Value

    tParams%expM=g_NoData_Value
    tParams%coeR=g_NoData_Value
    tParams%coeS=g_NoData_Value
    tParams%KS=g_NoData_Value
    tParams%KI=g_NoData_Value

  END WHERE
  RETURN
  END SUBROUTINE AssignValueToParamsType

!============================================================
!
!============================================================
! Revised: Yanjun Gan, 06/17/2015

  SUBROUTINE ReadCalib_SCEUAFile()
  USE CREST_Project
  USE CREST_Basic
  USE CREST_Calib_SCEUA
  IMPLICIT NONE
  CHARACTER(LEN=200)::sFileName,sTemp

  INTEGER :: error,i,j,k

  LOGICAL :: bIsColRow,bIsError
  DOUBLE PRECISION:: dblTemp

  sFileName=TRIM(g_CalibPath)//TRIM("Calibrations.txt")
  g_iseed=XXWReadLineInt(TRIM(sFileName),"iseed", error)
  g_maxn=XXWReadLineInt(TRIM(sFileName),"maxn", error)
  g_kstop=XXWReadLineInt(TRIM(sFileName),"kstop", error)
  g_pcento=XXWReadLineDbl(TRIM(sFileName),"pcento", error)
  g_ngs=XXWReadLineInt(TRIM(sFileName),"ngs", error)

  g_NCalibStas=XXWReadLineInt(TRIM(sFileName), &
    "NCalibStations", error)

  bIsColRow=XXWReadLineBln(TRIM(sFileName),"IsColRow", &
    error,"yes")

  ALLOCATE(g_tCalibSta(0:g_NCalibStas-1))


  DO i=0,g_NCalibStas-1
    WRITE(sTemp,*)i+1
    sTemp=ADJUSTL(sTemp)
    g_tCalibSta(i)%Name=XXWReadLineStr(TRIM(sFileName), &
      "Name_"//TRIM(sTemp), error)
    g_tCalibSta(i)%Value=XXWReadLineInt(TRIM(sFileName), &
      "Value_"//TRIM(sTemp), error)

    IF(bIsColRow .EQV. .TRUE.)THEN
      g_tCalibSta(i)%Col &
        =XXWReadLineInt(TRIM(sFileName), &
        "Col_"//TRIM(sTemp), error)

      g_tCalibSta(i)%Row &
        =XXWReadLineInt(TRIM(sFileName), &
        "Row_"//TRIM(sTemp), error)
    ELSE
      dblTemp &
        =XXWReadLineDbl(TRIM(sFileName), &
        "Long_"//TRIM(sTemp), error)
      g_tCalibSta(i)%Col &
        =INT((dblTemp-g_XLLCorner)/g_CellSize)

      dblTemp &
        =XXWReadLineDbl(TRIM(sFileName), &
        "Lati_"//TRIM(sTemp), error)
      g_tCalibSta(i)%Row &
        =INT((g_YLLCorner+g_NRows*g_CellSize-dblTemp) &
        /g_CellSize)

    END IF

    IF(InBasin2(g_tCalibSta(i)%Col,g_tCalibSta(i)%Row) &
      .EQV. .FALSE.)THEN
      g_tCalibSta(i)%bIsOut=.TRUE.
    ELSE
      g_tCalibSta(i)%bIsOut=.FALSE.
    END IF


    IF(g_tCalibSta(i)%bIsOut .EQV. .TRUE.)THEN
      CYCLE
    END IF

    ALLOCATE(g_tCalibSta(i)%a(1:UBOUND(g_sParamName,1)))
    ALLOCATE(g_tCalibSta(i)%bl(1:UBOUND(g_sParamName,1)))
    ALLOCATE(g_tCalibSta(i)%bu(1:UBOUND(g_sParamName,1)))
    ALLOCATE(g_tCalibSta(i)%x(1:UBOUND(g_sParamName,1)))

    g_tCalibSta(i)%a=1.0
    g_tCalibSta(i)%bl=1.0
    g_tCalibSta(i)%bu=1.0

    ! Input the Extent of Parameters
    DO j=1,UBOUND(g_sParamName,1)
      CALL XXWReadLineDbl3(TRIM(sFileName), &
        TRIM(g_sParamName(j))//"_"//TRIM(sTemp),error, &
        g_tCalibSta(i)%bl(j), &
        g_tCalibSta(i)%a(j), &
        g_tCalibSta(i)%bu(j))
        !print*,g_tCalibSta(i)%a(j)
    END DO

    g_tCalibSta(i)%x=g_tCalibSta(i)%a

    g_tCalibSta(i)%Area &
      =g_FAC(g_tCalibSta(i)%Col,g_tCalibSta(i)%Row)

    ALLOCATE(g_tCalibSta(i)%Mask(0:g_NCols-1,0:g_NRows-1))

    g_tCalibSta(i)%Mask=g_NoData_Value

    ! Get Station's Mask Matrix
    CALL GetMask(g_NCols,g_NRows,g_tCalibSta(i)%Col, &
      g_tCalibSta(i)%Row,g_NoData_Value,g_FDR, &
      g_NextC,g_NextR,g_tCalibSta(i)%Mask,InBasin)

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"Sta_" &
      //TRIM(g_tCalibSta(i)%Name)//"_Mask",  &
      g_tCalibSta(i)%Mask, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)

    ! Read Station's Observated Runoff
    ALLOCATE(g_tCalibSta(i)%RObs(0:g_ITMax-1))
    ALLOCATE(g_tCalibSta(i)%RSim(0:g_ITMax-1))

    CALL XXWReadRunoffObs(TRIM(g_OBSPath) &
      //TRIM(g_tCalibSta(i)%Name)//"_Obs.csv",  &
      g_tCalibSta(i)%RObs,bIsError)
    !print*,g_tCalibSta(i)%RObs

  END DO

  !Get g_RegMask Map
  ALLOCATE(g_RegMask(0:g_NCols-1,0:g_NRows-1))
  CALL ReadMatrixFile_Int(TRIM(g_CalibPath) // "CalibMask",   &
    g_RegMask,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_CalibFormat,"")

  IF(bIsError .EQV. .TRUE.)THEN
!+ygan
    g_RegMask=g_NoData_Value
!.ygan
    DO i=g_NCalibStas-1,0,-1
!-ygan
!      g_RegMask=g_NoData_Value
!.ygan
      WHERE(g_tCalibSta(i)%Mask==1)
        g_RegMask=g_tCalibSta(i)%Value
      END WHERE
    END DO
    CALL WriteMatrixFile_Int(TRIM(g_ResultPath) // "CalibMask",   &
      g_RegMask,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat)
  END IF

  RETURN
  END SUBROUTINE ReadCalib_SCEUAFile


!============================================================
!
!============================================================
  FUNCTION QUnitFactor(TimeMark)
  IMPLICIT NONE
  CHARACTER(1):: TimeMark
  DOUBLE PRECISION :: QUnitFactor

  QUnitFactor=1.0
  SELECT CASE (TimeMark)
  CASE ("y","Y"); !!Year
    QUnitFactor=1000.0/(3600.0*24.0*365.0)
  CASE ("m","M"); !!Month
    QUnitFactor=1000.0/(3600.0*24.0*30.0)
  CASE ("d","D"); !!Day
    QUnitFactor=1000.0/(3600.0*24.0)
  CASE ("h","H"); !!Hour
    QUnitFactor=1000.0/3600.0
  CASE ("u","U"); !!Minute
    QUnitFactor=1000.0/60.0
  CASE ("s","S"); !!Second
    QUnitFactor=1000.0
  END SELECT

  END FUNCTION QUnitFactor

!============================================================
!
!============================================================
  SUBROUTINE GetProjectFNameAndPath(sPFN,sPWS,sFNP)
  IMPLICIT NONE

  CHARACTER*(*) :: sPFN,sPWS,sFNP
  CHARACTER*200 :: COMLIN,sTemp
  LOGICAL :: EXISTS
  INTEGER :: ICOL,ISTART,ISTOP,N,iPoint,iBackSlash

  sFNP=' '
  COMLIN=' '
  sTemp=' '
  CALL GETARG(1,COMLIN)
  ICOL = 1
  IF(COMLIN.NE.' ') THEN
    sTemp=COMLIN
  ELSE
15  WRITE(*,*) ' Enter the name of the Project file: '
    READ(*,'(A)') COMLIN
    CALL URWORD(COMLIN,ICOL,ISTART,ISTOP,0,N,0,0)
    sTemp=COMLIN(ISTART:ISTOP)
    IF (sTemp.EQ.' ') GOTO 15
  ENDIF
  sFNP=TRIM(sTemp)
  INQUIRE (FILE=sFNP,EXIST=EXISTS)
  IF(.NOT.EXISTS) THEN
    sFNP=TRIM(sTemp)//".Project"
    INQUIRE (FILE=sFNP,EXIST=EXISTS)
    IF(.NOT.EXISTS) THEN
      WRITE (*,480) TRIM(sTemp),TRIM(sFNP)
480   FORMAT(1X,'Can not find Project file ',A,' or ',A)

      WRITE(*,*)"Press any key to continue!"
      READ(*,*)
      STOP
    ENDIF
  ENDIF
  iPoint=INDEX(sFNP,'.',.TRUE.)
  iBackSlash=INDEX(sFNP,'\\',.TRUE.)
  IF(iBackSlash==0)THEN
    iBackSlash=INDEX(sFNP,'/',.TRUE.)
  END IF
  sPFN=sFNP(iBackSlash+1:iPoint-1)
  sPWS=sFNP(1:iBackSlash)

  RETURN
  END SUBROUTINE GetProjectFNameAndPath

!============================================================
!
!============================================================
! Revised: Yanjun Gan, 06/01/2015

  SUBROUTINE CREST_Main_Pre(bIsError)

  USE CREST_Project
  USE CREST_Basic
  USE CREST_Param
  USE CREST_ICS

  IMPLICIT NONE

  INTEGER :: IBDT(8),IEDT(8)
  INTEGER :: k,i,j,ii,jj
  LOGICAL :: bIsError

  bIsError=.FALSE.
  WRITE(*,"(2X,A)")"Reading Project Data!"
  CALL ReadProjectFile(bIsError)

  IF(bIsError .EQV. .TRUE.)THEN
    WRITE(*,"(1X,A)") &
      "*** Something wrong in your project file!"
    WRITE(g_CREST_LogFileID,"(1X,A)") &
      "*** Something wrong in your project file!"
    RETURN
  END IF

  WRITE(*,"(2X,A)")"Reading Basic Data!"
  CALL ReadBasicFile(bIsError)
  IF(bIsError .EQV. .TRUE.)THEN
    WRITE(*,"(1X,A)") &
      "*** Something wrong in your basics file!"
    WRITE(g_CREST_LogFileID,"(1X,A)") &
      "*** Something wrong in your basics file!"
    RETURN
  END IF

  IF(g_tOutlet%HasOutlet)THEN
    g_tOutlet%Mask=g_NoData_Value
    ii=g_tOutlet%Row
    jj=g_tOutlet%Col

    WRITE(*,"(2X,A)")"Getting Mask Map of Outlet!"

    CALL GetMask(g_NCols,g_NRows,jj,ii,g_NoData_Value, &
      g_FDR,g_NextC,g_NextR,g_tOutlet%Mask,InBasin)

    WRITE(*,"(2X,A)")"Writing Mask Map of Outlet to File!"

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"Outlet_" &
      //TRIM(g_tOutlet%Name)//"_Mask",   &
      g_tOutlet%Mask, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)
!+ygan
    PRINT*,"  Outlet ",TRIM(g_tOutlet%Name)," located at Column ",jj,", Row ",ii,". FAC value is ",g_FAC(jj,ii)
!.ygan
  END IF

  DO k=0,g_NOutPixs-1
    g_tOutPix(k)%Mask=g_NoData_Value
    ii=g_tOutPix(k)%Row
    jj=g_tOutPix(k)%Col

    WRITE(*,"(2X,A,I2,A)") &
      "Getting Mask Map of OutPix Num ",k+1," !"

    CALL GetMask(g_NCols,g_NRows,jj,ii,g_NoData_Value, &
      g_FDR,g_NextC,g_NextR,g_tOutPix(k)%Mask,InBasin)

    WRITE(*,"(2X,A,I2,A)") &
      "Writing Mask Map of OutPix Num ",k+1," to File!"

    CALL WriteMatrixFile_Int(TRIM(g_ResultPath)//"OutPix_"   &
      //TRIM(g_tOutPix(k)%Name)//"_Mask",   &
      g_tOutPix(k)%Mask, g_NCols,g_NRows, &
      g_XLLCorner,g_YLLCorner, g_CellSize, g_NoData_Value, &
      bIsError,g_BasicFormat)

!+ygan
    PRINT*,"  OutPix ",TRIM(g_tOutPix(k)%Name)," located at Column ",jj,", Row ",ii,". FAC value is ",g_FAC(jj,ii)
!.ygan
  END DO

  WRITE(*,"(2X,A)")"Reading Prarmeters' Data!"
  CALL ReadParamFile()
  WRITE(*,"(2X,A)")"Reading Initial Conditions Data!"
  CALL ReadICSFile()

  SELECT CASE (TRIM(g_sRunStyle))
  CASE ("CALI_SCEUA")
    WRITE(*,"(2X,A)")"Reading calibration's Data!"
    CALL ReadCalib_SCEUAFile()
  END SELECT

  WRITE(*,*)
  RETURN
  END SUBROUTINE CREST_Main_Pre


!============================================================
!
!============================================================
  SUBROUTINE ReadProjectFile(bIsError)
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE

  INTEGER :: error,i
  CHARACTER(LEN=200):: sTemp
  CHARACTER(LEN=200):: sTemp_ForNum
  LOGICAL :: bOutPixColRow
  DOUBLE PRECISION:: dblTemp
  LOGICAL :: bIsError
  INTEGER :: dtTemp(1:6)

  bIsError=.FALSE.
  g_NCols=XXWReadLineInt(g_PrjNP,"NCols", error)
  g_NRows=XXWReadLineInt(g_PrjNP,"NRows", error)

  g_XLLCorner=XXWReadLineDbl(g_PrjNP,"XLLCorner", error)
  g_YLLCorner=XXWReadLineDbl(g_PrjNP,"YLLCorner", error)

  g_CellSize=XXWReadLineDbl(g_PrjNP,"CellSize", error)
  g_NoData_Value=XXWReadLineDbl(g_PrjNP,"NoData_Value", error)

  !WRITE(*,*)"---",g_NCols,g_NRows,g_XLLCorner,g_YLLCorner,g_CellSize,g_NoData_Value !--------------
  !Model Directory
  g_sRunStyle=XXWReadLineStr(g_PrjNP,"RunStyle", error)
  CALL UPCASE(g_sRunStyle)

  SELECT CASE(g_sRunStyle(1:4))
  CASE("SIMU")
    g_RunStyle=1
  CASE("CALI")
    g_RunStyle=2
  CASE("REAL")
    g_RunStyle=3
  CASE("REPE")
    g_RunStyle=4
  CASE("MCAT")
    g_RunStyle=5
  CASE("FORE")
    g_RunStyle=6
    CASE DEFAULT
    g_RunStyle=1
  END SELECT

  WRITE(*,*)"  "//"Reading Routing Type!"
  !===
  g_RoutingType=XXWReadLineStr(g_PrjNP,"RoutingType", error)
  CALL UPCASE(g_RoutingType)
  !WRITE(*,*)"1-"//g_RoutingType
  IF(error/=0)THEN
    g_RoutingType="JLR" !JLR CJR
  END IF
  !WRITE(*,*)g_RoutingType
  !READ(*,*)
  g_BasicFormat=XXWReadLineStr(g_PrjNP,"BasicFormat", error)
  CALL UPCASE(g_BasicFormat)
  g_BasicPath=XXWReadLineStr(g_PrjNP,"BasicPath", error)

  g_ParamFormat=XXWReadLineStr(g_PrjNP,"ParamFormat", error)
  CALL UPCASE(g_ParamFormat)
  g_ParamPath=XXWReadLineStr(g_PrjNP,"ParamPath", error)

  g_StateFormat=XXWReadLineStr(g_PrjNP,"StateFormat", error)
  CALL UPCASE(g_StateFormat)
  g_StatePath=XXWReadLineStr(g_PrjNP,"StatePath", error)

  g_ICSFormat=XXWReadLineStr(g_PrjNP,"ICSFormat", error)
  CALL UPCASE(g_ICSFormat)
  g_ICSPath=XXWReadLineStr(g_PrjNP,"ICSPath", error)

  g_OBSFormat=XXWReadLineStr(g_PrjNP,"OBSFormat", error)
  CALL UPCASE(g_OBSFormat)
  g_OBSPath=XXWReadLineStr(g_PrjNP,"OBSPath", error)

  g_RainFormat=XXWReadLineStr(g_PrjNP,"RainFormat", error)
  CALL UPCASE(g_RainFormat)
  g_RainPath=XXWReadLineStr(g_PrjNP,"RainPath", error)

  g_PETFormat=XXWReadLineStr(g_PrjNP,"PETFormat", error)
  CALL UPCASE(g_PETFormat)
  g_PETPath=XXWReadLineStr(g_PrjNP,"PETPath", error)

  g_ResultFormat=XXWReadLineStr(g_PrjNP,"ResultFormat", error)
  CALL UPCASE(g_PETFormat)
  g_ResultPath=XXWReadLineStr(g_PrjNP,"ResultPath", error)

  IF(g_sRunStyle(1:4)=="CALI")THEN
    g_CalibFormat=XXWReadLineStr(g_PrjNP,"CalibFormat",error)
    CALL UPCASE(g_CalibFormat)
    g_CalibPath=XXWReadLineStr(g_PrjNP,"CalibPath",error)
  END IF


  ! Model Run time Information
  g_TimeMark=XXWReadLineStr(g_PrjNP,"TimeMark", error)
  g_TimeStep=XXWReadLineDbl(g_PrjNP,"TimeStep", error)

  sTemp=XXWReadLineStr(g_PrjNP,"StartDate", error)
  CALL myStrtoD(sTemp,g_StartDate,g_TimeMark)

  sTemp=XXWReadLineStr(g_PrjNP,"WarmupDate", error)
  CALL myStrtoD(sTemp,g_WarmupDate,g_TimeMark)

  sTemp=XXWReadLineStr(g_PrjNP,"EndDate", error)
  CALL myStrtoD(sTemp,g_EndDate,g_TimeMark)

  g_NumWarmup=myDdiff(g_TimeMark, g_StartDate, g_WarmupDate) &
    / g_TimeStep
  g_ITMax=myDdiff(g_TimeMark, g_StartDate, g_EndDate) &
    / g_TimeStep + 1
  g_LoadState=XXWReadLineBln(g_PrjNP,"LoadState", error,"yes")

  g_SaveState=XXWReadLineBln(g_PrjNP,"SaveState", error,"yes")

  !WRITE(*,*)"---",g_NumWarmup,g_ITMax !--------------

  ! Outlet information
  g_tOutlet%HasOutlet=XXWReadLineBln(g_PrjNP,"HasOutlet", &
    error,"yes")
  IF(g_tOutlet%HasOutlet .NEQV. .TRUE.)THEN
    WRITE(*,*)"You did not assign a outlet!"
    GOTO 5000
  END IF
  bOutPixColRow=XXWReadLineBln(g_PrjNP,"OutletColRow", &
    error,"yes")
  g_tOutlet%Name &
    =XXWReadLineStr(g_PrjNP,"OutletName",error)

  IF(bOutPixColRow .EQV. .TRUE.)THEN
    g_tOutlet%Col &
      =XXWReadLineInt(g_PrjNP,"OutletCol", error)
    g_tOutlet%Row &
      =XXWReadLineInt(g_PrjNP,"OutletRow", error)
  ELSE
    dblTemp &
      =XXWReadLineDbl(g_PrjNP,"OutletLong", error)
    g_tOutlet%Col=INT((dblTemp-g_XLLCorner)/g_CellSize)

    dblTemp &
      =XXWReadLineDbl(g_PrjNP,"OutletLati", error)
    g_tOutlet%Row &
      =INT((g_YLLCorner+g_NRows*g_CellSize-dblTemp) &
      /g_CellSize)
  END IF

  IF(g_tOutlet%Col<0 &
    .OR. g_tOutlet%Col>=g_NCols &
    .OR. g_tOutlet%Row<0 &
    .OR. g_tOutlet%Row>=g_NRows)THEN
    g_tOutlet%bIsOut=.TRUE.
    WRITE(*,*) "Your outlet is out of the basin! " &
      //"Please check it!"
    STOP
  ELSE
    g_tOutlet%bIsOut=.FALSE.
  END IF

  ALLOCATE(g_tOutlet%RObs(0:g_ITMax-1))
  CALL XXWReadRunoffObs(TRIM(g_OBSPath) &
    //TRIM(g_tOutlet%Name)//"_Obs.csv", &
    g_tOutlet%RObs,bIsError)
  IF(bIsError .EQV. .TRUE.)THEN
    g_tOutlet%HasOBS=.FALSE.
    bIsError=.FALSE.
  ELSE
    g_tOutlet%HasOBS=.TRUE.
    bIsError=.FALSE.
  END IF
  ALLOCATE(g_tOutlet%Mask(0:g_NCols-1,0:g_NRows-1))

5000 CONTINUE

  !Output Pixels information
  g_NOutPixs=XXWReadLineInt(g_PrjNP,"NOutPixs", error)
  bOutPixColRow=XXWReadLineBln(g_PrjNP,"OutPixColRow", &
    error,"yes")

  ALLOCATE(g_tOutPix(0:g_NOutPixs-1))

  Do i=LBOUND(g_tOutPix,1), UBOUND(g_tOutPix,1)
    WRITE(sTemp_ForNum,"(I2)")i+1
    sTemp_ForNum=ADJUSTL(sTemp_ForNum)

    g_tOutPix(i)%Name &
      =XXWReadLineStr(g_PrjNP,"OutPixName" &
      //TRIM(sTemp_ForNum),error)

    IF(bOutPixColRow .EQV. .TRUE.)THEN
      g_tOutPix(i)%Col &
        =XXWReadLineInt(g_PrjNP,"OutPixCol" &
        //sTemp_ForNum, error)
      g_tOutPix(i)%Row &
        =XXWReadLineInt(g_PrjNP,"OutPixRow" &
        //sTemp_ForNum, error)
    ELSE
      dblTemp &
        =XXWReadLineDbl(g_PrjNP, &
        "OutPixLong" // sTemp_ForNum, error)
      g_tOutPix(i)%Col=INT((dblTemp-g_XLLCorner)/g_CellSize)

      dblTemp &
        =XXWReadLineDbl(g_PrjNP, &
        "OutPixLati"//sTemp_ForNum, error)
      g_tOutPix(i)%Row &
        =INT((g_YLLCorner+g_NRows*g_CellSize-dblTemp) &
        /g_CellSize)
    END IF

    IF(g_tOutPix(i)%Col<0 &
      .OR. g_tOutPix(i)%Col>=g_NCols &
      .OR. g_tOutPix(i)%Row<0 &
      .OR. g_tOutPix(i)%Row>=g_NRows)THEN
      g_tOutPix(i)%bIsOut=.TRUE.
    ELSE
      g_tOutPix(i)%bIsOut=.FALSE.
    END IF

    ALLOCATE(g_tOutPix(i)%RObs(0:g_ITMax-1))
    CALL XXWReadRunoffObs(TRIM(g_OBSPath) &
      //TRIM(g_tOutPix(i)%Name)//"_Obs.csv", &
      g_tOutPix(i)%RObs,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      g_tOutPix(i)%HasOBS=.FALSE.
      bIsError=.FALSE.
    ELSE
      g_tOutPix(i)%HasOBS=.TRUE.
      bIsError=.FALSE.
    END IF

    ALLOCATE(g_tOutPix(i)%Mask(0:g_NCols-1,0:g_NRows-1))
  END DO


  !Grid Outputs information
  ALLOCATE(g_bGOVar(LBOUND(g_sGOVarName,1) &
    :UBOUND(g_sGOVarName,1)))
  DO i=LBOUND(g_sGOVarName,1),UBOUND(g_sGOVarName,1)
    g_bGOVar(i)=XXWReadLineBln(g_PrjNP,"GOVar_" &
      //TRIM(g_sGOVarName(i)), &
      error,"yes")
  END DO

  ! Output the specified date
  g_NOutDTs=XXWReadLineInt(g_PrjNP,"NumOfOutputDates", error)

  ALLOCATE(g_OutDTIndex(0:g_NOutDTs-1))

  Do i=LBOUND(g_OutDTIndex,1), UBOUND(g_OutDTIndex,1)
    WRITE(sTemp_ForNum,"(I2)")i+1
    sTemp_ForNum=ADJUSTL(sTemp_ForNum)

    sTemp=XXWReadLineStr(g_PrjNP,"OutputDate_" &
      // sTemp_ForNum, error)
    CALL myStrtoD(sTemp,dtTemp,g_TimeMark)
    g_OutDTIndex(i)=myDdiff(g_TimeMark, g_StartDate, dtTemp) &
      / g_TimeStep
  END DO

  bIsError=.FALSE.

  !Read Lake Mask Data Added by Xianwu Xue 2011.4.21
  g_NLakes=XXWReadLineInt(g_PrjNP,"NumOfLakes", error)

  IF(g_NLakes>0)THEN
    ALLOCATE(g_LakeMask(0:g_NCols-1,0:g_NRows-1))

    CALL ReadMatrixFile_Int(TRIM(g_BasicPath) // "LakeMask",   &
      g_LakeMask,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_BasicFormat,"")
    IF(bIsError .EQV. .TRUE.)THEN
      WRITE(*,"(1X,A)") &
        "*** Something wrong in your Lake Mask file!"
      WRITE(g_CREST_LogFileID,"(1X,A)") &
        "*** Something wrong in your Lake Mask file!"
      RETURN
    END IF
  END IF


  RETURN
  END SUBROUTINE ReadProjectFile

!============================================================
!
!============================================================
  INTEGER FUNCTION XXWReadLineInt(sFileNP,sKey,error)
  IMPLICIT NONE
  CHARACTER(LEN=200)::sKey2
  CHARACTER*(*)::sKey,sFileNP
  INTEGER :: fileid
  INTEGER, INTent(out) :: error
  CHARACTER(LEN=200)::sTemp,sTemp2
  INTEGER :: ICOL=1,ISTART,ISTOP,NCODE=1,N,IOUT,IN
  INTEGER :: ISTART_in,ISTOP_in

  XXWReadLineInt=0 ! Added by Xianwu Xue 2011.4.21

  IOUT=-1
  IN=-1
  sKey2=sKey
  ICOL=1
  CALL URWORD(sKey2,ICOL,ISTART_in,ISTOP_in,NCODE,N,IOUT,IN)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0
  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error)sTemp
    ! WRITE(*,*)sTemp(1:1)
    IF (sTemp(1:1)=="#")THEN
      CYCLE
    END IF
    ICOL=1
    CALL URWORD(sTemp,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)
    IF(sTemp(ISTART:ISTOP)==sKey2(ISTART_in:ISTOP_in)) THEN
      BACKSPACE(fileid)
      READ(fileid,*,IOSTAT=error)sTemp,sTemp2,XXWReadLineInt
      CLOSE(fileid)
      RETURN
    END IF
  END DO

  CLOSE(fileid)

  RETURN
  END FUNCTION XXWReadLineInt

!============================================================
!
!============================================================
  DOUBLE PRECISION FUNCTION XXWReadLineDbl(sFileNP,sKey,error)
  IMPLICIT NONE
  CHARACTER(LEN=200)::sKey2
  CHARACTER*(*)::sKey,sFileNP
  INTEGER :: fileid
  INTEGER, INTent(out) :: error
  CHARACTER(LEN=200)::sTemp,sTemp2
  INTEGER :: ICOL=1,ISTART,ISTOP,NCODE=1,N,IOUT,IN
  INTEGER :: ISTART_in,ISTOP_in

  XXWReadLineDbl=0 ! Added by Xianwu Xue 2011.4.21

  IOUT=-1
  IN=-1
  sKey2=sKey
  ICOL=1
  CALL URWORD(sKey2,ICOL,ISTART_in,ISTOP_in,NCODE,N,IOUT,IN)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0
  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error)sTemp
    ! WRITE(*,*)sTemp(1:1)
    IF (sTemp(1:1)=="#")THEN
      CYCLE
    END IF
    ICOL=1
    CALL URWORD(sTemp,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)
    IF(sTemp(ISTART:ISTOP)==sKey2(ISTART_in:ISTOP_in)) THEN
      BACKSPACE(fileid)
      READ(fileid,*,IOSTAT=error)sTemp,sTemp2,XXWReadLineDbl
      CLOSE(fileid)
      RETURN
    END IF
  END DO

  CLOSE(fileid)

  RETURN
  END FUNCTION XXWReadLineDbl

!============================================================
!
!============================================================
  SUBROUTINE XXWReadLineDbl3(sFileNP,sKey,error,bl,a,bu)
  IMPLICIT NONE
  CHARACTER(LEN=200)::sKey2
  CHARACTER*(*)::sKey,sFileNP
  INTEGER :: fileid
  INTEGER, INTent(out) :: error
  CHARACTER(LEN=200)::sTemp,sTemp2
  INTEGER :: ICOL=1,ISTART,ISTOP,NCODE=1,N,IOUT,IN
  INTEGER :: ISTART_in,ISTOP_in
  DOUBLE PRECISION :: bl,a,bu

  IOUT=-1
  IN=-1
  sKey2=sKey
  ICOL=1
  CALL URWORD(sKey2,ICOL,ISTART_in,ISTOP_in,NCODE,N,IOUT,IN)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0
  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error)sTemp
    ! WRITE(*,*)sTemp(1:1)
    IF (sTemp(1:1)=="#")THEN
      CYCLE
    END IF
    ICOL=1
    CALL URWORD(sTemp,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)
    IF(sTemp(ISTART:ISTOP)==sKey2(ISTART_in:ISTOP_in)) THEN
      BACKSPACE(fileid)
      READ(fileid,*,IOSTAT=error)sTemp, &
        sTemp2,bl,a,bu
      CLOSE(fileid)
      RETURN
    END IF
  END DO

  CLOSE(fileid)
  !If it can't find the sKey's value
  bl=1.0
  a=1.0
  bu=1.0

  RETURN
  END SUBROUTINE XXWReadLineDbl3

!============================================================
!
!============================================================
  CHARACTER(LEN=200) FUNCTION XXWReadLineStr(sFileNP,sKey,error)
  IMPLICIT NONE
  CHARACTER(LEN=200)::sKey2
  CHARACTER*(*)::sKey,sFileNP
  INTEGER :: fileid
  INTEGER, INTent(out) :: error
  CHARACTER(LEN=200)::sTemp,sTemp2
  INTEGER :: ICOL=1,ISTART,ISTOP,NCODE=1,N,IOUT,IN
  INTEGER :: ISTART_in,ISTOP_in

  XXWReadLineStr=""

  IOUT=-1
  IN=-1
  sKey2=sKey
  ICOL=1
  CALL URWORD(sKey2,ICOL,ISTART_in,ISTOP_in,NCODE,N,IOUT,IN)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0
  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error)sTemp
    ! WRITE(*,*)sTemp(1:1)
    IF (sTemp(1:1)=="#")THEN
      CYCLE
    END IF
    ICOL=1
    CALL URWORD(sTemp,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)
    IF(sTemp(ISTART:ISTOP)==sKey2(ISTART_in:ISTOP_in)) THEN
      BACKSPACE(fileid)
      READ(fileid,*,IOSTAT=error)sTemp,sTemp2,XXWReadLineStr
      CLOSE(fileid)
      RETURN
    END IF
  END DO

  CLOSE(fileid)

  RETURN
  END FUNCTION XXWReadLineStr

!============================================================
!
!============================================================
  LOGICAL FUNCTION XXWReadLineBln(sFileNP,sKey,error,sKeyValue)
  IMPLICIT NONE
  CHARACTER(LEN=200)::sKey2,sKeyValue2
  CHARACTER*(*)::sKey,sKeyValue,sFileNP
  INTEGER :: fileid
  INTEGER, INTent(out) :: error
  CHARACTER(LEN=200)::sTemp,sTemp2,sTemp3
  INTEGER :: ICOL=1,ISTART,ISTOP,NCODE=1,N,IOUT,IN
  INTEGER :: ISTART_in,ISTOP_in

  sKeyValue2=sKeyValue
  CALL UPCASE(sKeyValue2)

  IOUT=-1
  IN=-1
  sKey2=sKey
  ICOL=1
  CALL URWORD(sKey2,ICOL,ISTART_in,ISTOP_in,NCODE,N,IOUT,IN)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(sFileNP), FORM="FORMATTED")

  error=0
  ! rewind(fileid,IOSTAT=error)
  DO WHILE (error==0)
    READ(fileid,*,IOSTAT=error)sTemp
    ! WRITE(*,*)sTemp(1:1)
    IF (sTemp(1:1)=="#")THEN
      CYCLE
    END IF
    ICOL=1
    CALL URWORD(sTemp,ICOL,ISTART,ISTOP,NCODE,N,IOUT,IN)
    IF(sTemp(ISTART:ISTOP)==sKey2(ISTART_in:ISTOP_in)) THEN
      BACKSPACE(fileid)
      READ(fileid,*,IOSTAT=error)sTemp,sTemp2,sTemp3
      CALL UPCASE(sTemp3)

      IF(TRIM(sTemp3)==TRIM(sKeyValue2))THEN
        XXWReadLineBln=.TRUE.
      ELSE
        XXWReadLineBln=.FALSE.
      END IF
      CLOSE(fileid)
      RETURN
    END IF
  END DO

  CLOSE(fileid)

  RETURN
  END FUNCTION XXWReadLineBln


!============================================================
!
!============================================================
! Revised: Yanjun Gan, 03/09/2015

  SUBROUTINE ReadMatrixFile(FileName_In, dblMat, &
    NCols_In, NRows_In,XLLCorner_In,YLLCorner_In, &
    CellSize_In,NoData_Value_In,bIsError,sFileFormat_In, &
    strDate_In)

  USE CREST_Project
  IMPLICIT NONE

  CHARACTER(*) :: FileName_In,sFileFormat_In
  CHARACTER(*) :: strDate_In
  CHARACTER(LEN=200):: FileName
  CHARACTER(14):: strDate,strTemp
  INTEGER :: NCols_In,NRows_In,intTemp
  DOUBLE PRECISION :: XLLCorner_In,YLLCorner_In
  DOUBLE PRECISION :: CellSize_In,NoData_Value_In
  LOGICAL :: bIsError,fExist

  INTEGER :: NCols, NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value

  DOUBLE PRECISION :: dblMat(0:NCols_In-1,0:NRows_In-1)
  DOUBLE PRECISION,ALLOCATABLE :: dblMatBig(:,:)
  CHARACTER(LEN=10) :: sFileFormat

  INTEGER :: dtTemp(1:6),iNum
  CHARACTER(20):: sTemp

  sFileFormat=TRIM(ADJUSTL(sFileFormat_In))
  strDate=TRIM(strDate_In)

  CALL UPCASE(sFileFormat)

  bIsError=.FALSE.

  !WRITE(*,*)"Format: ",TRIM(sFileFormat_In), TRIM(sFileFormat)
  SELECT CASE (TRIM(sFileFormat))

    ! #########################################################################
  CASE ("ASBIMO") !For PET--------------------------------------------
    FileName=TRIM(FileName_In)//"PET025." &
      //TRIM(strDate(5:6))// ".asc"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF
    CALL ReadASCFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    !WRITE(*,*)TRIM(FileName),NCols,NRows

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadASCFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF
    ! #########################################################################
  CASE ("BIBIMO") !For PET--------------------------------------------
    FileName=TRIM(FileName_In)//"PET025." &
      //TRIM(strDate(5:6))// ".bif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadBIFFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadBIFFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

!+ygan
    ! #########################################################################
  CASE ("DPET") ! FEWS NET Daily PET-------------------------------------------
    FileName=TRIM(FileName_In)//"et"//TRIM(strDate(3:8))// ".txt"
    !print*,FileName
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF
    CALL ReadASCFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    !WRITE(*,*)TRIM(FileName),NCols,NRows

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadASCFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF
!.ygan

    ! #########################################################################
  CASE ("ASC","TXT") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)// ".asc"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In)//TRIM(strDate)// ".txt"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF

    CALL ReadASCFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    !WRITE(*,*)bIsError
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadASCFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("DBIF") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)// ".dbif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadDBIFFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadDBIFFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("BIFFIT") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)//".bif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadBIFFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadBIFFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("BBIF") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)//".bbif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadBBIFFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadBBIFFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("TRMMRT") !--------------------------------------------
    FileName=TRIM(FileName_In) // "3B42RT." &
      //TRIM(strDate(1:10))//".7R2.bin"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In) // "3B42RT." &
        //TRIM(strDate(1:10))//".7.bin"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        FileName=TRIM(FileName_In) // "3B42RT."  &
          //TRIM(strDate(1:10))// ".6A.bin"
        INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
        IF(fExist .EQV. .FALSE.)THEN
          FileName=TRIM(FileName_In) // "3B42RT."  &
            //TRIM(strDate(1:10))// ".6.bin"
          INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
          IF(fExist .EQV. .FALSE.)THEN
            FileName=TRIM(FileName_In) // "3B42RT." &
              //TRIM(strDate(1:10))// ".bin"
            INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
            IF(fExist .EQV. .FALSE.)THEN
              bIsError=.TRUE.
              RETURN
            END IF
          END IF
        END IF
      END IF
    END IF

    CALL ReadTRMMRTFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadTRMMRTFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("TRMMV6") !--------------------------------------------
    strTemp=TRIM(strDate(9:10))
    READ(strTemp,*) intTemp
    IF(intTemp<10)THEN
      WRITE(strTemp,"(I1)")intTemp
    ELSE
      WRITE(strTemp,"(I2)")intTemp
    END IF

    FileName=TRIM(FileName_In) // "3B42." &
      //TRIM(strDate(3:8))//"."//TRIM(strTemp)//".6.HDF"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In) // "3B42."  &
        //TRIM(strDate(3:8))//"."//TRIM(strTemp)// ".6A.HDF"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF

    CALL ReadTRMMV6FileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadTRMMV6File(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("TRMMV7") !--------------------------------------------
    FileName=TRIM(FileName_In) // "3B42." &
      //TRIM(strDate(1:8))//"."//TRIM(strDate(9:10))//".7.HDF"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadTRMMV6FileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadTRMMV6File(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

    ! #########################################################################
  CASE ("NMQBIN") !------XXW: Not Validate--------------------------------------
    strTemp=TRIM(strDate(9:10))
    FileName=TRIM(FileName_In) // "1HRAD.HSR."  &
      //TRIM(strDate(1:8))//"."//TRIM(strTemp)//"00"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadNMQBINFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadNMQBINFile(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF

!+ygan
    ! #########################################################################
  CASE ("DIAMOND4") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)//".000"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In)//TRIM(strDate)//".txt"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF

    CALL ReadDiamond4FileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadDiamond4File(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
       DEALLOCATE(dblMatBig)
       RETURN
    END IF

    ! #########################################################################
  CASE ("DIAMOND131") !-------------------------------------------- SWAN QPE QPF
    FileName=TRIM(FileName_In)//TRIM(strDate)//"0000.bin"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    CALL ReadDiamond131FileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(dblMatBig(0:NCols-1,0:NRows-1))

    CALL ReadDiamond131File(TRIM(FileName), dblMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      DEALLOCATE(dblMatBig)
      RETURN
    END IF
!.ygan

  CASE DEFAULT

    WRITE(*,*) "ERROR!!! This version does not use this format!"
    bIsError=.TRUE.
    !DEALLOCATE(dblMatBig)
    RETURN
  END SELECT

! Interpolation
  IF(NCols_In /=NCols .OR. NRows_In /=NRows &
    .OR. CellSize_In /=CellSize)THEN

    CALL GetBfromA(NCols, NRows, XLLCorner, YLLCorner, &
      CellSize, NoData_Value, dblMatBig, &
      NCols_In, NRows_In, &
      XLLCorner_In, YLLCorner_In, &
      CellSize_In, NoData_Value_In, dblMat)

  ELSE
    dblMat=dblMatBig
  END IF

  IF(NoData_Value_In/=NoData_Value)THEN
    WHERE(dblMat==NoData_Value)
      dblMat=NoData_Value_In
    END WHERE
  END IF

  DEALLOCATE(dblMatBig)

  RETURN
  END SUBROUTINE ReadMatrixFile


  ! Change resolution and span, including interpolation function based on resampling technique
!============================================================
!
!============================================================
  SUBROUTINE GetBfromA(ANCols, ANRows, AXLLCorner, AYLLCorner, &
    ACellSize, ANoData_Value, AdblMat, &
    BNCols, BNRows, BXLLCorner, BYLLCorner, &
    BCellSize, BNoData_Value, BdblMat)
  IMPLICIT NONE
  INTEGER*4:: ANCols, ANRows, BNCols, BNRows, i, j,souR, souC
  DOUBLE PRECISION :: AXLLCorner, AYLLCorner
  DOUBLE PRECISION :: ACellSize, ANoData_Value
  DOUBLE PRECISION :: AdblMat(0:AnCols-1, 0:AnRows-1)

  DOUBLE PRECISION :: BXLLCorner, BYLLCorner
  DOUBLE PRECISION :: BCellSize, BNoData_Value
  DOUBLE PRECISION :: BdblMat(0:BnCols-1, 0:BnRows-1)
  DOUBLE PRECISION :: SR, SC


  BdblMat=BNoData_Value
  DO i=0, BNRows-1
    SR = BYLLCorner + (BNRows-i-0.5) * BCellSize
    souR = (AYLLCorner +ANRows * ACellSize - SR) / ACellSize
    IF (souR>=0 .AND. souR < ANRows) THEN
      DO j=0, BNCols-1
        SC = BXLLCorner+(j + 0.5) * BCellSize
        souC = (SC - AXLLCorner) / ACellSize
        IF (souC>=0 .AND. souC < ANCols) THEN
          IF (AdblMat(souC, souR) /= ANoData_Value) THEN
            BdblMat(j, i) = AdblMat(souC, souR)
          END IF
        END IF
      ENDDO
    END IF
  ENDDO
  RETURN
  END SUBROUTINE GetBfromA



!============================================================
!
!============================================================
  SUBROUTINE GetBfromA_Int(ANCols, ANRows, AXLLCorner, AYLLCorner, &
    ACellSize, ANoData_Value, AdblMat, &
    BNCols, BNRows, BXLLCorner, BYLLCorner, &
    BCellSize, BNoData_Value, BdblMat)
  IMPLICIT NONE
  INTEGER*4:: ANCols, ANRows, BNCols, BNRows, i, j,souR, souC
  DOUBLE PRECISION :: AXLLCorner, AYLLCorner
  DOUBLE PRECISION :: ACellSize, ANoData_Value
  INTEGER :: AdblMat(0:AnCols-1, 0:AnRows-1)

  DOUBLE PRECISION :: BXLLCorner, BYLLCorner
  DOUBLE PRECISION :: BCellSize, BNoData_Value
  INTEGER :: BdblMat(0:BnCols-1, 0:BnRows-1)
  DOUBLE PRECISION :: SR, SC


  BdblMat=BNoData_Value
  DO i=0, BNRows-1
    SR = BYLLCorner + (BNRows-i-0.5) * BCellSize
    souR = (AYLLCorner +ANRows * ACellSize - SR) / ACellSize
    IF (souR>=0 .AND. souR < ANRows) THEN
      DO j=0, BNCols-1
        SC = BXLLCorner+(j + 0.5) * BCellSize
        souC = (SC - AXLLCorner) / ACellSize
        IF (souC>=0 .AND. souC < ANCols) THEN
          IF (AdblMat(souC, souR) /= ANoData_Value) THEN
            BdblMat(j, i) = AdblMat(souC, souR)
          END IF
        END IF
      ENDDO
    END IF
  ENDDO
  RETURN
  END SUBROUTINE GetBfromA_Int

!============================================================
!
!============================================================
  SUBROUTINE ReadMatrixFile_Int(FileName_In, IntMat, &
    NCols_In, NRows_In,XLLCorner_In,YLLCorner_In, &
    CellSize_In,NoData_Value_In,bIsError,sFileFormat_In, &
    strDate_In)
  IMPLICIT NONE


  CHARACTER(*):: FileName_In,sFileFormat_In,strDate_In
  CHARACTER(LEN=200):: FileName
  CHARACTER(14):: strDate,strTemp
  INTEGER :: NCols_In,NRows_In
  DOUBLE PRECISION :: XLLCorner_In,YLLCorner_In
  DOUBLE PRECISION :: CellSize_In,NoData_Value_In
  LOGICAL :: bIsError,fExist

  INTEGER :: NCols, NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value

  INTEGER :: IntMat(0:NCols_In-1,0:NRows_In-1)
  INTEGER,ALLOCATABLE :: IntMatBig(:,:)
  CHARACTER(LEN=10) :: sFileFormat

  sFileFormat=sFileFormat_In
  strDate=TRIM(strDate_In)

  CALL UPCASE(sFileFormat)

  bIsError=.FALSE.

  SELECT CASE (TRIM(sFileFormat))
  CASE ("ASC","TXT") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)// ".asc"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In)//TRIM(strDate)// ".txt"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF

    CALL ReadASCFileHeader(TRIM(FileName), &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      RETURN
    END IF

    ALLOCATE(IntMatBig(0:NCols-1,0:NRows-1))

    CALL ReadASCFile_Int(TRIM(FileName),IntMatBig, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
    IF(bIsError .EQV. .TRUE.)THEN
      WRITE(*,*) "Some errors in your file: " // TRIM(FileName)
      DEALLOCATE(IntMatBig)
      RETURN
    END IF
    CASE DEFAULT
  END SELECT

  IF(NCols_In /=NCols .OR. NRows_In /=NRows)THEN
    CALL GetBfromA_Int(NCols, NRows, XLLCorner, YLLCorner, &
      CellSize, NoData_Value, IntMatBig, &
      NCols_In, NRows_In, &
      XLLCorner_In, YLLCorner_In, &
      CellSize_In, NoData_Value_In, IntMat)
  ELSE
    IntMat=IntMatBig
  END IF

  IF(NoData_Value_In/=NoData_Value)THEN
    WHERE(IntMat==NoData_Value)
      IntMat=NoData_Value_In
    END WHERE
  END IF

  RETURN
  END SUBROUTINE ReadMatrixFile_Int


!============================================================
!
!============================================================
  SUBROUTINE ReadASCFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)

  IMPLICIT NONE
  CHARACTER(*):: FileName
  LOGICAL :: bIsError
  INTEGER :: fileid

  INTEGER :: NCols,NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner
  DOUBLE PRECISION :: CellSize,NoData_Value

  CHARACTER(20):: sTemp

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),FORM='FORMATTED')
  READ(fileid,*) sTemp, NCols
  !WRITE(*,*) NCols

  READ(fileid,*) sTemp, NRows
  !WRITE(*,*) NRows

  READ(fileid,*) sTemp, XLLCorner
  !WRITE(*,*) XLLCorner

  READ(fileid,*) sTemp, YLLCorner
  !WRITE(*,*) YLLCorner

  READ(fileid,*) sTemp, CellSize
  !WRITE(*,*) CellSize

  READ(fileid,*) sTemp, NoData_Value
  !WRITE(*,*) NoData_Value

  CLOSE(fileid)

  RETURN

  END SUBROUTINE ReadASCFileHeader


!============================================================
!
!============================================================
  SUBROUTINE ReadASCFile(FileName,dblMat, &
    NCols,NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE
  CHARACTER(*):: FileName
  LOGICAL :: bIsError
  INTEGER :: fileid
  INTEGER :: i,j

  INTEGER :: NCols,NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner
  DOUBLE PRECISION :: CellSize,NoData_Value

  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)

  CHARACTER(20):: sTemp

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),FORM='FORMATTED')
  READ(fileid,*) sTemp, NCols
  READ(fileid,*) sTemp, NRows
  READ(fileid,*) sTemp, XLLCorner
  READ(fileid,*) sTemp, YLLCorner
  READ(fileid,*) sTemp, CellSize
  READ(fileid,*) sTemp, NoData_Value

  READ(fileid,*)((dblMat(j,i),j=0,NCols-1),i=0,NRows-1)

  CLOSE(fileid)
  RETURN
  END SUBROUTINE ReadASCFile

!============================================================
!
!============================================================
  SUBROUTINE ReadASCFile_Int(FileName,IntMat, &
    NCols,NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)

  IMPLICIT NONE
  CHARACTER(*):: FileName
  LOGICAL :: bIsError
  INTEGER :: fileid
  INTEGER :: i,j

  INTEGER :: NCols,NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner
  DOUBLE PRECISION :: CellSize,NoData_Value

  INTEGER :: IntMat(0:NCols-1,0:NRows-1)

  CHARACTER(20):: sTemp

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),FORM='FORMATTED')
  READ(fileid,*) sTemp, NCols
  READ(fileid,*) sTemp, NRows
  READ(fileid,*) sTemp, XLLCorner
  READ(fileid,*) sTemp, YLLCorner
  READ(fileid,*) sTemp, CellSize
  READ(fileid,*) sTemp, NoData_Value

  READ(fileid,*)((IntMat(j,i),j=0,NCols-1),i=0,NRows-1)

  CLOSE(fileid)
  RETURN
  END SUBROUTINE ReadASCFile_Int


!============================================================
!
!============================================================
! Revised: Yanjun Gan, 05/04/2015

  SUBROUTINE WriteASCFile(FileName, dblMat, nCols, nRows, XLLCorner, &
    YLLCorner, CellSize, NoData_Value,bIsError)
  IMPLICIT NONE
  CHARACTER(*):: FileName
  INTEGER :: lC, lR, nCols, nRows
  DOUBLE PRECISION :: dblMat(0:nCols-1, 0:nRows-1)
  DOUBLE PRECISION :: XLLCorner, YLLCorner
  DOUBLE PRECISION :: CellSize, NoData_Value
  INTEGER :: FileID
  LOGICAL :: bIsError

  bIsError=.FALSE.
  CALL XXWGetFreeFile(FileID)

  OPEN(FileID,FILE=TRIM(FileName),FORM='FORMATTED')
  WRITE(FileID,"('nCols 		')",ADVANCE='no')
  WRITE(FileID,"(i8)")nCols
  WRITE(FileID,"('nRows 		')",ADVANCE='no')
  WRITE(FileID,"(i8)")nRows
  WRITE(FileID,"('XLLCorner 	')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")XLLCorner
  WRITE(FileID,"('YLLCorner 	')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")YLLCorner
  WRITE(FileID,"('cellSize 	 ')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")CellSize
  WRITE(FileID,"('NODATA_value ')",ADVANCE='no')
  WRITE(FileID,"(f11.0)")NoData_Value
  DO lR=0, nRows-1
    DO lC=0, nCols-1
!-ygan
      !WRITE(FileID,'(F21.15,A)',ADVANCE='no') &
!.ygan
!+ygan
      WRITE(FileID,'(F12.5,A)',ADVANCE='no') &
!.ygan
        dblMat(lC,lR),' ';
    END DO
    WRITE(FileID,*)
  END DO
  CLOSE (FileID)
  END SUBROUTINE WriteASCFile

!============================================================
!
!============================================================
  SUBROUTINE WriteASCFile_Int(FileName, IntMat, nCols, nRows, &
    XLLCorner,YLLCorner, CellSize, NoData_Value,bIsError)
  IMPLICIT NONE
  CHARACTER(*):: FileName
  INTEGER :: lC, lR, nCols, nRows
  INTEGER :: IntMat(0:nCols-1, 0:nRows-1)
  DOUBLE PRECISION :: XLLCorner, YLLCorner
  DOUBLE PRECISION :: CellSize, NoData_Value
  INTEGER :: FileID
  LOGICAL :: bIsError

  bIsError=.FALSE.
  CALL XXWGetFreeFile(FileID)

  OPEN(FileID,FILE=TRIM(FileName),FORM='FORMATTED')
  WRITE(FileID,"('nCols 		')",ADVANCE='no')
  WRITE(FileID,"(i8)")nCols
  WRITE(FileID,"('nRows 		')",ADVANCE='no')
  WRITE(FileID,"(i8)")nRows
  WRITE(FileID,"('XLLCorner 	')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")XLLCorner
  WRITE(FileID,"('YLLCorner 	')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")YLLCorner
  WRITE(FileID,"('cellSize 	 ')",ADVANCE='no')
  WRITE(FileID,"(f11.6)")CellSize
  WRITE(FileID,"('NODATA_value ')",ADVANCE='no')
  WRITE(FileID,"(f11.0)")NoData_Value
  DO lR=0, nRows-1
    DO lC=0, nCols-1
      WRITE(FileID,'(I10,A)',ADVANCE='no')IntMat(lC,lR),' ';
    END DO
    WRITE(FileID,*)
  END DO
  CLOSE (FileID)
  END SUBROUTINE WriteASCFile_Int

!============================================================
!
!============================================================
  SUBROUTINE ReadBBIFFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  REAL(KIND=4) :: XLLCorner_R,YLLCorner_R
  REAL(KIND=4) :: CellSize_R,NoData_Value_R

  LOGICAL :: bIsError

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*4
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R

  ! ALL SwapInt4(NCols)
  ! ALL SwapInt4(NRows)
  !
  ! ALL SwapSng(XLLCorner_R)
  ! ALL SwapSng(YLLCorner_R)
  !
  ! ALL SwapSng(CellSize_R)
  ! ALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R

  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  !WRITE(*,*)"--Test--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value
  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadBBIFFileHeader


!============================================================
!
!============================================================
  SUBROUTINE ReadBBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  REAL(KIND=4) :: VOut2(0:NCols-1,0:NRows-1), &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R,test
  LOGICAL :: bIsError
  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*4+NCols*NRows*4 !6*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, &
    ((VOut2(j,i), j=0, NCols-1), i=0, NRows-1)

  !WRITE(*,*)TRIM(FileName)
  ! ALL SwapInt4(NCols)
  ! ALL SwapInt4(NRows)
  !
  ! ALL SwapSng(XLLCorner_R)
  ! ALL SwapSng(YLLCorner_R)
  !
  ! ALL SwapSng(CellSize_R)
  ! ALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R
  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  ! DO i=0, NRows-1
  !DO j=0, NCols-1
  ! CALL SwapSng(VOut2(j,i))
  !END DO
  ! END DO

  dblMat=VOut2


  !test=VOut2(174,0)
  ! ALL SwapSng(test)
  !WRITE(*,*)"--Test-2--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value


  !WRITE(*,*) (VOut2(j,0), j=0, NCols-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadBBIFFile



!============================================================
!
!============================================================
  SUBROUTINE WriteBBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  LOGICAL :: bIsError

  REAL(KIND=4) :: dblMat2(0:NCols-1,0:NRows-1), &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)
  XLLCorner_R=XLLCorner
  YLLCorner_R=YLLCorner
  CellSize_R=CellSize
  NoData_Value_R=NoData_Value
  dblMat2=dblMat

  blocklength=2*4+4*4+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), STATUS="REPLACE", &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)

  WRITE(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, &
    ((dblMat2(j,i), j=0, NCols-1), i=0, NRows-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE WriteBBIFFile


!============================================================
!
!============================================================
  SUBROUTINE ReadDBIFFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp

  LOGICAL :: bIsError

  bIsError=.TRUE.

  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*8
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner,YLLCorner, &
    CellSize,NoData_Value


  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadDBIFFileHeader

!============================================================
!
!============================================================
  SUBROUTINE ReadDBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*8+26+NCols*NRows*8
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner,YLLCorner, &
    CellSize,NoData_Value, EmptySpace, &
    ((dblMat(j,i), j=0, NCols-1), i=0, NRows-1)
  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadDBIFFile

!============================================================
!
!============================================================
  SUBROUTINE WriteDBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*8+26+NCols*NRows*8
  OPEN(fileid, FILE=TRIM(FileName), STATUS="REPLACE", &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)

  WRITE(fileid, rec=1)NCols, NRows, XLLCorner,YLLCorner, &
    CellSize,NoData_Value, EmptySpace, &
    ((dblMat(j,i), j=0, NCols-1), i=0, NRows-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE WriteDBIFFile


!============================================================
!
!============================================================
  SUBROUTINE ReadBIFFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  REAL(KIND=4) :: XLLCorner_R,YLLCorner_R
  REAL(KIND=4) :: CellSize_R,NoData_Value_R

  LOGICAL :: bIsError

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*4
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R

  ! ALL SwapInt4(NCols)
  ! ALL SwapInt4(NRows)
  !
  ! ALL SwapSng(XLLCorner_R)
  ! ALL SwapSng(YLLCorner_R)
  !
  ! ALL SwapSng(CellSize_R)
  ! ALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R

  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  !WRITE(*,*)"--Test--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value
  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadBIFFileHeader

!============================================================
!
!============================================================
  SUBROUTINE ReadBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  REAL(KIND=4) :: VOut2(0:NCols-1,0:NRows-1), &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R,test
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError
  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*4+26+NCols*NRows*4 !6*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, EmptySpace, &
    ((VOut2(j,i), j=0, NCols-1), i=0, NRows-1)

  !WRITE(*,*)TRIM(FileName)
  ! ALL SwapInt4(NCols)
  ! ALL SwapInt4(NRows)
  !
  ! ALL SwapSng(XLLCorner_R)
  ! ALL SwapSng(YLLCorner_R)
  !
  ! ALL SwapSng(CellSize_R)
  ! ALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R
  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  ! DO i=0, NRows-1
  !DO j=0, NCols-1
  ! CALL SwapSng(VOut2(j,i))
  !END DO
  ! END DO

  dblMat=VOut2


  !test=VOut2(174,0)
  ! ALL SwapSng(test)
  !WRITE(*,*)"--Test-2--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value


  !WRITE(*,*) (VOut2(j,0), j=0, NCols-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadBIFFile

!============================================================
!
!============================================================
  SUBROUTINE ReadBIFFile_Int(FileName,intMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  INTEGER :: intMat(0:NCols-1,0:NRows-1)
  REAL(KIND=4) :: &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError
  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  blocklength=2*4+4*4+26+NCols*NRows*4 !6*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength, &
    STATUS='OLD',action='READ')

  READ(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, EmptySpace, &
    ((intMat(j,i), j=0, NCols-1), i=0, NRows-1)
  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R
  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadBIFFile_Int

!============================================================
!
!============================================================
  SUBROUTINE WriteBIFFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError

  REAL(KIND=4) :: dblMat2(0:NCols-1,0:NRows-1), &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)
  XLLCorner_R=XLLCorner
  YLLCorner_R=YLLCorner
  CellSize_R=CellSize
  NoData_Value_R=NoData_Value
  dblMat2=dblMat

  blocklength=2*4+4*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), STATUS="REPLACE", &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)

  WRITE(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, EmptySpace, &
    ((dblMat2(j,i), j=0, NCols-1), i=0, NRows-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE WriteBIFFile

!============================================================
!
!============================================================
  SUBROUTINE WriteBIFFile_Int(FileName,intMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  INTEGER :: intMat(0:NCols-1,0:NRows-1)
  CHARACTER*26:: EmptySpace
  LOGICAL :: bIsError

  REAL(KIND=4) :: &
    XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)
  XLLCorner_R=XLLCorner
  YLLCorner_R=YLLCorner
  CellSize_R=CellSize
  NoData_Value_R=NoData_Value

  blocklength=2*4+4*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), STATUS="REPLACE", &
    ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength, &
    action='WRITE')

  WRITE(fileid, rec=1)NCols, NRows, XLLCorner_R,YLLCorner_R, &
    CellSize_R,NoData_Value_R, EmptySpace, &
    ((intMat(j,i), j=0, NCols-1), i=0, NRows-1)

  CLOSE (fileid)
  RETURN
  END SUBROUTINE WriteBIFFile_Int

!============================================================
!
!============================================================
  SUBROUTINE ReadNMQBINFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  INTEGER*4:: Year, Month, Day, Hour, Minute, Second
  INTEGER*4:: souNZ

  LOGICAL :: bIsError


  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  OPEN (fileid, FILE=TRIM(FileName), ACCESS='DIRECT', &
    FORM='UNFORMATTED', RECL=4)
  READ(fileid, rec=1) Year
  READ(fileid, rec=2) Month
  READ(fileid, rec=3) Day
  READ(fileid, rec=4) Hour
  READ(fileid, rec=5) Minute
  READ(fileid, rec=6) Second
  READ(fileid, rec=7) NCols
  READ(fileid, rec=8) NRows
  READ(fileid, rec=9) souNZ
  CLOSE(fileid)

  RETURN
  END SUBROUTINE ReadNMQBINFileHeader

!============================================================
!
!============================================================
  SUBROUTINE ReadNMQBINFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)

  INTEGER*2, ALLOCATABLE:: ShortData(:,:)

  LOGICAL :: bIsError
  INTEGER*4:: Year, Month, Day, Hour, Minute, Second
  INTEGER*4:: souNZ
  INTEGER*4:: mapProject, truLat1, truLat2, truLon, xy_scaLE, &
    map_scaLE_temp, nw_lon_temp, nw_lat_temp
  INTEGER*4:: dx, dy, dxy_scaLE, z_height, z_scaLE, filLEr, &
    varScaLETemp, varMissingTemp, numRad
  CHARACTER*1:: varName(20), varUnit(6), radName(4)
  REAL*4:: nw_lat, nw_lon, map_scaLE, varScaLE, varMissing
  REAL*4:: souCS, souND

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  BlockLenGTh = (NCols * NRows * 2) + (34 * 4) + 30 + 4
  ALLOCATE(ShortData(1:NCols,1:NRows))

  OPEN (fileid, FILE=TRIM(FileName), ACCESS='DIRECT', &
    FORM='UNFORMATTED', RECL=blockLenGTh)
  READ(fileid, rec=1) Year, Month, Day, Hour, Minute, Second, &
    NCols, NRows, souNZ, mapProject, map_scaLE_temp, &
    truLat1, truLat2, truLon, nw_lon_temp, nw_lat_temp, xy_scaLE, &
    dx, dy, dxy_scaLE, z_height, z_scaLE, filLEr, filLEr, filLEr, &
    filLEr, filLEr, filLEr, filLEr, filLEr, filLEr, filLEr, &
    varName, varUnit, varScaLETemp, varMissingTemp, &
    numRad, radName, ShortData

  CLOSE(fileid)

  ! onvert our ints to floats and scale things properly
  map_scaLE = REAL(map_scaLE_temp)
  nw_lon = REAL(nw_lon_temp) / map_scaLE
  nw_lat = REAL(nw_lat_temp) / map_scaLE
  varScaLE = REAL(varScaLETemp)
  varMissing = REAL(varMissingTemp)
  DO i = 0, NRows-1
    DO j = 0, NCols-1
      dblMat(j, i) = REAL(ShortData(i+1,(NRows - j)))/varScaLE
      IF (dblMat(j, i) .LT. 0) THEN
        dblMat(j, i) = varMissing
      ENDIF
    END DO
  END DO
  XLLCorner = nw_lon
  YLLCorner = nw_lat-((REAL(dy)/REAL(dxy_scaLE)) * REAL(NRows))
  CellSize = REAL(dx) / REAL(dxy_scaLE)
  NoData_Value = varMissing
  RETURN
  END SUBROUTINE ReadNMQBINFile

!+ygan
!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE ReadDiamond4FileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp

  INTEGER(4) :: year,month,day,hour,min,time,dmdcode
  INTEGER(4) :: smoothv,heavierv

  INTEGER(4) irec,iost,tmpI4

  DOUBLE PRECISION :: intv,minv,maxv
  DOUBLE PRECISION :: xRes,yRes
  DOUBLE PRECISION :: stX,endX,stY,endY
  CHARACTER(LEN=10)::dmdstr
  CHARACTER(LEN=80)::discrip

  LOGICAL :: bIsError


  bIsError=.FALSE.

  NCols=0
  NRows=0
  NoData_Value = -9999

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),STATUS='OLD', IOSTAT=iost, ERR=99)

  READ(fileid,*,ERR=98,IOSTAT=iost) dmdstr,dmdcode,discrip
  READ(fileid,*,ERR=98,IOSTAT=iost) year,month,day,hour,min,time, &
  xRes,yRes,stX,endX,stY,endY,NCols,NRows, &
  intv,minv,maxv,smoothv,heavierv

  !Check the order of rows and columes saved in file
  CellSize = xRes
  IF(xRes.LT.0.) THEN
    XLLCorner = endX
  ELSE
    XLLCorner = stX
  ENDIF
  IF(yRes.LT.0.) THEN
    YLLCorner = endY
  ELSE
    YLLCorner = stY
  ENDIF

  IF(DABS(xRes).NE.DABS(yRes)) THEN
    WRITE(*,*) " X gridsize and Y gridsize are different!"
    GOTO 98
  ENDIF

98 CLOSE (fileid)

  !WRITE(*,*)"--Test--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value

   RETURN

99 WRITE(*,*) 'Open Diamond4 file error, ERR=', TRIM(FileName)
  RETURN

  END SUBROUTINE ReadDiamond4FileHeader



!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE ReadDiamond4File(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)

  INTEGER(4) :: year,month,day,hour,min,time,dmdcode
  INTEGER(4) :: smoothv,heavierv

  INTEGER(4) irec,iost,tmpI4

  DOUBLE PRECISION :: intv,minv,maxv
  DOUBLE PRECISION :: xRes,yRes
  DOUBLE PRECISION :: stX,endX,stY,endY
  CHARACTER(LEN=10)::dmdstr
  CHARACTER(LEN=80)::discrip

  LOGICAL :: bIsError


  bIsError=.FALSE.

  NoData_Value = -9999

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),STATUS='OLD', IOSTAT=iost, ERR=99)

  READ(fileid,*,ERR=98,IOSTAT=iost) dmdstr,dmdcode,discrip
  READ(fileid,*,ERR=98,IOSTAT=iost) year,month,day,hour,min,time, &
  xRes,yRes,stX,endX,stY,endY,NCols,NRows, &
  intv,minv,maxv,smoothv,heavierv

  !Check the order of rows and columes saved in file
  CellSize = xRes
  IF(xRes.LT.0.) THEN
    XLLCorner = endX
  ELSE
    XLLCorner = stX
  ENDIF
  IF(yRes.LT.0.) THEN
    YLLCorner = endY
  ELSE
    YLLCorner = stY
  ENDIF

  IF(DABS(xRes).NE.DABS(yRes)) THEN
    WRITE(*,*) " X gridsize and Y gridsize are different!"
    GOTO 98
  ENDIF

  dblMat=NoData_Value

  !READ(fileid,*,ERR=98,IOSTAT=iost)((dblMat(j,i),j=0,NCols-1),i=NRows-1,0,-1)
  READ(fileid,*,ERR=98,IOSTAT=iost)((dblMat(j,i),j=0,NCols-1),i=0,NRows-1)

  WHERE(dblMat<0)
    dblMat=NoData_Value
  END WHERE

98 CLOSE (fileid)

  !WRITE(*,*)"--Test--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value

  RETURN

99 WRITE(*,*) 'Open Diamond4 file error, ERR=', TRIM(FileName)
  RETURN

  END SUBROUTINE ReadDiamond4File

!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE WriteDiamond4File(FileName,dblMat, &
  NCols, NRows,XLLCorner,YLLCorner, &
  CellSize,NoData_Value,bIsError,strDate,strProdDescrip)
  IMPLICIT NONE

  CHARACTER(*)::strDate,strProdDescrip !for diamond
  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)

  INTEGER(4) :: year,month,day,hour,min,time,dmdcode
  INTEGER(4) :: smoothv,heavierv

  INTEGER(4) irec,iost,tmpI4
  DOUBLE PRECISION :: intv,minv,maxv
  DOUBLE PRECISION :: xRes,yRes
  DOUBLE PRECISION :: stX,endX,stY,endY
  CHARACTER(LEN=10)::dmdstr
  CHARACTER(LEN=80)::discrip

  CHARACTER(10) :: cDataFmt
  CHARACTER(5) :: cHeadFmt
  CHARACTER(62) :: cDescFmt
  CHARACTER(256):: cFullFileName

  LOGICAL :: bIsError

  bIsError=.FALSE.

  NoData_Value = -9999
  dmdstr="diamond 4  "
  discrip=strDate//strProdDescrip

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),STATUS='UNKNOWN',FORM='FORMATTED', IOSTAT=iost, ERR=99)

  xRes = CellSize != xRes
  stX = XLLCorner != stX
  endX= XLLCorner+CellSize*(NCols-1)
  yRes= -xRes
  stY = YLLCorner+CellSize*(NRows-1)
  endY= YLLCorner

  intv=0.1
  minv=0.
  maxv=0.
  smoothv=0
  heavierv=0

  DO i=0,NRows-1
    DO j=0,NCols-1
      !IF(dblMat(j,i).EQ.NoData_Value) dblMat(j,i)=0.
      IF(maxv.LT.dblMat(j,i)) maxv=dblMat(j,i)
    ENDDO
  ENDDO

  cHeadFmt = '(A,A)'
  WRITE(fileid,cHeadFmt,ERR=98,IOSTAT=iost)'diamond 4 ', TRIM(discrip)
  cDescFmt = '(A4,4(" ",A2)," 0 ",2f10.5,4f10.4,2i8,3f8.2,2(i4))'
  WRITE(fileid,TRIM(cDescFmt),ERR=98,IOSTAT=iost) strDate(1:4),strDate(5:6),strDate(7:8),strDate(9:10),' 0', &
  xRes,yRes,stX,endX,stY,endY,NCols,NRows, &
  intv,minv,maxv,smoothv,heavierv

  !DO i=NRows-1,0,-1
  DO i=0,NRows-1
    DO j=0,NCols-1
      WRITE(fileid,'(F12.2,A)',ADVANCE='no') dblMat(j,i),' ';
      IF(MOD(j+1,20).EQ.0)WRITE(fileid,*)
    END DO
    WRITE(fileid,*)
  END DO

98 CLOSE (fileid)
  RETURN

99 WRITE(*,*) 'Open Diamond4 file error, ERR=', TRIM(FileName)
  RETURN

  END SUBROUTINE WriteDiamond4File

!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE WriteDiamond4File_Int(FileName,intMat, &
  NCols, NRows,XLLCorner,YLLCorner, &
  CellSize,NoData_Value,bIsError,strDate,strProdDescrip)
  IMPLICIT NONE

  CHARACTER(*)::strDate,strProdDescrip !for diamond
  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  INTEGER :: intMat(0:NCols-1,0:NRows-1)

  INTEGER(4) :: year,month,day,hour,min,time,dmdcode
  INTEGER(4) :: smoothv,heavierv
  INTEGER(4) :: irec,iost,tmpI4

  DOUBLE PRECISION :: intv,minv,maxv
  DOUBLE PRECISION :: xRes,yRes
  DOUBLE PRECISION :: stX,endX,stY,endY
  CHARACTER(LEN=10):: dmdstr
  CHARACTER(LEN=80):: discrip
  LOGICAL :: bIsError

  CHARACTER(10) :: cDataFmt
  CHARACTER(5) :: cHeadFmt
  CHARACTER(62) :: cDescFmt
  CHARACTER(256):: cFullFileName

  !LOGICAL :: bIsError

  bIsError=.FALSE.
  NoData_Value = -9999
  dmdstr="diamond 4  "
  discrip=TRIM(strDate)//TRIM(strProdDescrip)

  CALL XXWGetFreeFile(fileid)

  OPEN(fileid,FILE=TRIM(FileName),STATUS='UNKNOWN',FORM='FORMATTED', IOSTAT=iost, ERR=99)
  xRes = CellSize != xRes
  stX = XLLCorner != stX
  endX= XLLCorner+CellSize*(NCols-1)
  yRes= -xRes
  stY = YLLCorner+CellSize*(NRows-1)
  endY= YLLCorner

  intv=0.1
  minv=0.
  maxv=0.
  smoothv=0
  heavierv=0

  DO i=0,NRows-1
    DO j=0,NCols-1
      !IF(dblMat(j,i).EQ.NoData_Value) dblMat(j,i)=0.
      IF(maxv.LT.intMat(j,i)) maxv=intMat(j,i)
    ENDDO
  ENDDO

  cHeadFmt = '(A,A)'
  WRITE(fileid,cHeadFmt,ERR=98,IOSTAT=iost)'diamond 4 ', TRIM(discrip)
  cDescFmt = '(A4,4(" ",A2)," 0 ",2f10.5,4f10.4,2i8,3f8.2,2(i4))'
  WRITE(fileid,TRIM(cDescFmt),ERR=98,IOSTAT=iost) strDate(1:4),strDate(5:6),strDate(7:8),strDate(9:10),' 0', &
  xRes,yRes,stX,endX,stY,endY,NCols,NRows, &
  intv,minv,maxv,smoothv,heavierv

  !cDataFmt = '(2(f10.3))'
  DO i=NRows-1, 0,-1
    DO j=0, NCols-1
      !IF(NoData_Value.EQ.dblMat(j,i)) intMat(j,i)=0.
      WRITE(fileid,'(i12,A)',ADVANCE='no') intMat(j,i),' ';
      IF(MOD(j+1,20).EQ.0)WRITE(fileid,*)
    END DO
    WRITE(fileid,*)
  END DO

98 CLOSE (fileid)

  RETURN
99 WRITE(*,*) 'Open Diamond4 file error, ERR=', TRIM(FileName)
  RETURN
  END SUBROUTINE WriteDiamond4File_Int

!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE ReadDiamond131FileHeader(FileName, &
  NCols, NRows, XLLCorner, YLLCorner, &
  CellSize,NoData_Value,bIsError)
  USE CREST_Diamond
  IMPLICIT NONE
  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  REAL(KIND=4) :: XLLCorner_R,YLLCorner_R
  REAL(KIND=4) :: CellSize_R,NoData_Value_R

  TYPE (diamond131_header)::mheader

  LOGICAL :: bIsError

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  blocklength=1024
  OPEN(fileid, FILE=TRIM(FileName), &
  ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid,rec=1) mheader

  IF(SCAN(mheader%Flag,'SWAN')<=0 .AND. SCAN(mheader%Flag,'swan')<=0)THEN
    WRITE(*,*) " Not SWAN Mosaic QP Product!"
    GOTO 98
  ENDIF

  IF(mheader%YReso.NE.mheader%XReso) THEN
    WRITE(*,*) " X gridsize and Y gridsize are different!"
    GOTO 98
  ENDIF

  !year   =mheader%year
  !month  =mheader%month
  !day    =mheader%day
  !hour   =mheader%hour
  !min    =mheader%minute

  CellSize_R =mheader%XReso
  !CALL SwapSng(CellSize_R)
  NCols =mheader%XNumGrids
  NRows =mheader%YNumGrids
  XLLCorner_R =mheader%StartLon
  YLLCorner_R =mheader%StartLat-NRows*CellSize_R

  NoData_Value_R=-9999

  !CALL SwapInt4(NCols)
  !CALL SwapInt4(NRows)
  !
  !CALL SwapSng(XLLCorner_R)
  !CALL SwapSng(YLLCorner_R)
  !
  !CALL SwapSng(CellSize_R)
  !CALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R

  CellSize=CellSize_R
  NoData_Value=NoData_Value_R

  !WRITE(*,*)"--Test--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value
98 CONTINUE

  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadDiamond131FileHeader


!============================================================
!
!============================================================
! Author: Yanjun Gan, 04/10/2015

  SUBROUTINE ReadDiamond131File(FileName,dblMat, &
  NCols, NRows,XLLCorner,YLLCorner, &
  CellSize,NoData_Value,bIsError)
  USE CREST_Diamond
  IMPLICIT NONE

  TYPE (diamond131_header):: mheader

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j,i1
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  INTEGER(2),ALLOCATABLE::i2qpe_1h(:,:) !~0:NCols-1,0:NRows-1)
  !REAL(KIND=4) :: VOut2(0:NCols-1,0:NRows-1),
  REAL(KIND=4) :: XLLCorner_R,YLLCorner_R,CellSize_R,NoData_Value_R,test
  LOGICAL :: bIsError

  bIsError=.FALSE.

  CALL XXWGetFreeFile(fileid)

  ALLOCATE(i2qpe_1h(0:NCols-1,0:NRows-1))

  blocklength=1024+NCols*NRows*2 !6*4+26+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), &
  ACCESS='DIRECT', FORM='UNFORMATTED', RECL=blocklength)
  READ(fileid, rec=1) mheader,((i2qpe_1h(j,i), j=0, NCols-1), i=0, NRows-1)
  !NCols, NRows, XLLCorner_R,YLLCorner_R, &
  ! CellSize_R,NoData_Value_R, &
  ! ((VOut2(j,i), j=0, NCols-1), i=0, NRows-1)

  IF(SCAN(mheader%Flag,'SWAN')<=0 .AND. SCAN(mheader%Flag,'swan')<=0)THEN
    WRITE(*,*) " Not SWAN Mosaic QP Product!"
    GOTO 98
  ENDIF

  IF(mheader%YReso.NE.mheader%XReso) THEN
    WRITE(*,*) " X gridsize and Y gridsize are different!"
    GOTO 98
  ENDIF

  CellSize_R =mheader%XReso
  !CALL SwapSng(CellSize_R)
  NCols =mheader%XNumGrids
  NRows =mheader%YNumGrids
  XLLCorner_R =mheader%StartLon
  YLLCorner_R =mheader%StartLat-NRows*CellSize_R

  !NoData_Value_R=-9999

  !CALL SwapInt4(NCols)
  !CALL SwapInt4(NRows)
  !
  !CALL SwapSng(XLLCorner_R)
  !CALL SwapSng(YLLCorner_R)
  !
  !CALL SwapSng(CellSize_R)
  !CALL SwapSng(NoData_Value_R)

  XLLCorner=XLLCorner_R
  YLLCorner=YLLCorner_R

  CellSize=CellSize_R
  !NoData_Value=NoData_Value_R

  !DO i=0, NRows-1
  !  DO j=0, NCols-1
  !    CALL SwapSng(VOut2(j,i))
  !  END DO
  !END DO

  dblMat=NoData_Value !VOut2
  !copy DATA for saving
  DO i=0,NRows-1
    i1= i !NRows-1-i
    DO j=0,NCols-1
      IF(i2qpe_1h(j,i1)>diam131_no_value .AND. i2qpe_1h(j,i1)<diam131_max_qpe) THEN
        dblMat(j,i)=i2qpe_1h(j,i1)*0.1 !0.1mm GrdData(0:NCols-1,i)=tempGrd(0:NCols-1,i1)*0.1 !0.1mm
      ENDIF
    ENDDO
  ENDDO

  DEALLOCATE(i2qpe_1h)
  !test=VOut2(174,0)
  !CALL SwapSng(test)
  !WRITE(*,*)"--Test-2--",NCols, NRows, XLLCorner,YLLCorner, &
  ! CellSize,NoData_Value


  !WRITE(*,*) (VOut2(j,0), j=0, NCols-1)

98 CONTINUE

  CLOSE (fileid)
  RETURN
  END SUBROUTINE ReadDiamond131File
!.ygan

!============================================================
!
!============================================================
  SUBROUTINE ReadTRMMRTFileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp

  LOGICAL :: bIsError

  bIsError=.FALSE.
  NRows=480
  NCols=1440
  NoData_Value=-9999
  XLLCorner=-180.
  YLLCorner=-60.
  CellSize=0.25

  RETURN
  END SUBROUTINE ReadTRMMRTFileHeader



!============================================================
!
!============================================================
  SUBROUTINE ReadTRMMRTFile(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER*2880:: EmptySpace
  LOGICAL :: bIsError

  INTEGER*4:: gC

  INTEGER*2, ALLOCATABLE:: Idir(:,:)
  REAL*4, ALLOCATABLE:: temS(:,:)
  ALLOCATE(Idir(0:NCols-1, 0:NRows-1), temS(0:NCols-1, 0:NRows-1))

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  BlockLenGTh=2880+nCols*nRows*2

  OPEN(fileid, FILE=TRIM(FileName), ACCESS='DIRECT', &
    FORM='UNFORMATTED', RECL=BlockLenGTh)
  READ(fileid, rec=1)EmptySpace, &
    ((Idir(j,i), j=0, NCols-1), i=0, NRows-1)
  CLOSE(fileid)

  DO i=0, NRows-1
    DO j=0, NCols-1
      CALL SwapINT(Idir(j,i))
    END DO
  END DO

  temS=NoData_Value
  DO i=0, NRows-1
    DO j=0, NCols-1
      IF(Idir(j, i)>=0)THEN
        temS(j, i) = REAL(Idir(j, i)) / 100.0
      ELSEIF (Idir(j, i)==-1)THEN
        temS(j, i) = 0.0
      ELSEIF (Idir(j, i)==-31999)THEN
        temS(j, i) = NoData_Value
      ELSE
        temS(j, i)=-REAL(Idir(j, i))/100.0
      END IF
    ENDDO
  ENDDO

  DO j = 0, NCols-1
    IF (j < 720) THEN
      gC = j + 720
    ELSE
      gC = j - 720
    END IF
    DO i = 0, NRows-1
      dblMat(gC, i) = temS(j, i)
    ENDDO
  ENDDO

  RETURN
  END SUBROUTINE ReadTRMMRTFile

!============================================================
!
!============================================================
  SUBROUTINE SwapINT(Vsou)
  IMPLICIT NONE
  INTEGER*2:: Vsou, temI
  CHARACTER*1:: cha2(1:2), cha
  EquivaLENce (temI, cha2)
  temI=Vsou
  cha=cha2(1)
  cha2(1)=cha2(2)
  cha2(2)=cha
  Vsou=temI
  END SUBROUTINE SwapINT

!============================================================
!
!============================================================
  SUBROUTINE ReadTRMMV6FileHeader(FileName, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,blocklength,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp

  LOGICAL :: bIsError

  bIsError=.FALSE.
  NRows=400
  NCols=1440
  NoData_Value=-9999
  XLLCorner=-180.
  YLLCorner=-50.
  CellSize=0.25

  RETURN
  END SUBROUTINE ReadTRMMV6FileHeader



!============================================================
!
!============================================================
  SUBROUTINE ReadTRMMV6File(FileName,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError)
  IMPLICIT NONE

  CHARACTER(*):: FileName
  INTEGER :: fileid
  INTEGER :: NCols, NRows,BlockLenGTh,i,j
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value
  CHARACTER(20):: sTemp
  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER*298:: EmptySpace
  LOGICAL :: bIsError
  REAL(4), ALLOCATABLE:: temS(:,:),temSTemp(:)

  bIsError=.FALSE.
  CALL XXWGetFreeFile(fileid)

  ALLOCATE(temS(0:NRows-1, 0:NCols-1))
  BlockLenGTh=298+NCols*NRows*4
  OPEN(fileid, FILE=TRIM(FileName), ACCESS='DIRECT', &
    FORM='UNFORMATTED', RECL=BlockLenGTh)

  ! Developed by Tiger
  ! READ(fileid, rec=1)EmptySpace,
  ! & ((temS(i,j), i=0, NRows-1), j=0, NCols-1) !Tiger
  ! Developed by Tiger

  ! Developed by Xianwu xue 2010.12.26
  ALLOCATE(temSTemp(0:NCols-1))
  READ(fileid, rec=1)EmptySpace, &
    ((temS(i,j), i=1, NRows-1),temSTemp(j), j=0, NCols-1)

  DO j=0,NCols-1
    temS(0,j)=temSTemp(j)
  END DO
  ! Developed by Xianwu xue 2010.12.26

  CLOSE(fileid)

  DO i=0,NRows-1
    DO j=0,NCols-1
      CALL SwapSng(temS(i,j))
      IF (temS(i,j)<0) temS(i,j)=NoData_Value
    ENDDO
  ENDDO

  DO i=0,NRows-1
    DO j=0,NCols-1
      dblMat(j,NRows-1-i)=temS(i,j)
    ENDDO
  ENDDO

  RETURN
  END SUBROUTINE ReadTRMMV6File

!============================================================
!
!============================================================
  SUBROUTINE SwapSng(Vsou)
  IMPLICIT NONE
  REAL(4):: Vsou, temS
  CHARACTER*1:: cha2(1:4), cha
  EquivaLENce (temS, cha2)
  temS=Vsou
  cha=cha2(1)
  cha2(1)=cha2(4)
  cha2(4)=cha
  cha=cha2(2)
  cha2(2)=cha2(3)
  cha2(3)=cha
  Vsou=temS
  END SUBROUTINE SwapSng

!============================================================
!
!============================================================
  SUBROUTINE SwapInt4(Vsou)
  IMPLICIT NONE
  INTEGER*4:: Vsou, temS
  CHARACTER*1:: cha2(1:4), cha
  EquivaLENce (temS, cha2)
  temS=Vsou
  cha=cha2(1)
  cha2(1)=cha2(4)
  cha2(4)=cha
  cha=cha2(2)
  cha2(2)=cha2(3)
  cha2(3)=cha
  Vsou=temS
  END SUBROUTINE SwapInt4

!============================================================
!
!============================================================
  SUBROUTINE WriteMatrixFile(FileName_In,dblMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError,sFileFormat_In)
  IMPLICIT NONE

  CHARACTER(*):: FileName_In,sFileFormat_In
  CHARACTER(LEN=200):: FileName

  LOGICAL :: bIsError

  INTEGER :: NCols, NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value

  DOUBLE PRECISION :: dblMat(0:NCols-1,0:NRows-1)
  CHARACTER(LEN=10) :: sFileFormat

  sFileFormat=TRIM(sFileFormat_In)
  CALL UPCASE(sFileFormat)

  bIsError=.FALSE.

  SELECT CASE (TRIM(sFileFormat))
  CASE ("ASC")
    FileName=TRIM(FileName_In) // ".asc"
    CALL WriteASCFile(TRIM(FileName), dblMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("TXT")
    FileName=TRIM(FileName_In) // ".txt"
    CALL WriteASCFile(TRIM(FileName), dblMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("BIFFIT")
    FileName=TRIM(FileName_In) // ".bif"
    CALL WriteBIFFile(TRIM(FileName), dblMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("DBIF")
    FileName=TRIM(FileName_In) // ".dbif"
    CALL WriteDBIFFile(TRIM(FileName), dblMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)

  CASE ("BBIF")
    FileName=TRIM(FileName_In) // ".bbif"
    CALL WriteBBIFFile(TRIM(FileName), dblMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)

    CASE DEFAULT


  END SELECT


  RETURN
  END SUBROUTINE WriteMatrixFile


!============================================================
!
!============================================================
  SUBROUTINE WriteMatrixFile_Int(FileName_In,IntMat, &
    NCols, NRows,XLLCorner,YLLCorner, &
    CellSize,NoData_Value,bIsError,sFileFormat_In)
  IMPLICIT NONE

  CHARACTER(*):: FileName_In,sFileFormat_In
  CHARACTER(LEN=200):: FileName

  LOGICAL :: bIsError

  INTEGER :: NCols, NRows
  DOUBLE PRECISION :: XLLCorner,YLLCorner,CellSize,NoData_Value

  INTEGER :: IntMat(0:NCols-1,0:NRows-1)
  CHARACTER(LEN=10) :: sFileFormat

  sFileFormat=sFileFormat_In
  CALL UPCASE(sFileFormat)

  bIsError=.FALSE.

  SELECT CASE (TRIM(sFileFormat))
  CASE ("ASC")
    FileName=TRIM(FileName_In) // ".asc"
    CALL WriteASCFile_Int(TRIM(FileName), IntMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("TXT")
    FileName=TRIM(FileName_In) // ".txt"
    CALL WriteASCFile_Int(TRIM(FileName), IntMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("BIFFIT")
    FileName=TRIM(FileName_In) // ".bif"
    CALL WriteBIFFile_Int(TRIM(FileName), IntMat, &
      NCols, NRows,XLLCorner,YLLCorner, &
      CellSize,NoData_Value,bIsError)
  CASE ("DBIF")
  CASE ("BBIF")

  CASE DEFAULT


  END SELECT


  RETURN
  END SUBROUTINE WriteMatrixFile_Int

!============================================================
!
!============================================================
  SUBROUTINE XXWGetFreeFile(fileid)
  IMPLICIT NONE
  INTEGER::i,fileid
  LOGICAL::bOpen
  bOpen = .TRUE.

  DO i = 15, 100000
    INQUIRE(UNIT=i, OPENED=bOpen)
    IF (.NOT. bOpen) THEN
      fileid = i
      RETURN
    END IF
  END DO
  fileid = -1
  END SUBROUTINE XXWGetFreeFile



!============================================================
!
!============================================================
! Revised: Yanjun Gan, 04/10/2015

  SUBROUTINE InquireMatrixFile(FileName_In,bIsError,sFileFormat_In, &
    strDate_In)
  IMPLICIT NONE

  CHARACTER(*):: FileName_In,sFileFormat_In
  CHARACTER(*) :: strDate_In
  CHARACTER(LEN=200):: FileName
  CHARACTER(14):: strDate,strTemp
  INTEGER :: intTemp
  LOGICAL :: bIsError,fExist
  CHARACTER(LEN=10) :: sFileFormat

  sFileFormat=TRIM(ADJUSTL(sFileFormat_In))
  strDate=TRIM(strDate_In)

  CALL UPCASE(sFileFormat)

  bIsError=.FALSE.

  SELECT CASE (TRIM(sFileFormat))
    ! #########################################################################
  CASE ("ASBIMO") !For PET--------------------------------------------
    FileName=TRIM(FileName_In)//"PET025."  &
      //TRIM(strDate(5:6))// ".txt"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    ! #########################################################################
  CASE ("BIBIMO") !For PET--------------------------------------------
    FileName=TRIM(FileName_In)//"PET025." &
      //TRIM(strDate(5:6))// ".bif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    ! #########################################################################
  CASE ("ASC","TXT") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)// ".asc"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In)//TRIM(strDate)// ".txt"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF

    ! #########################################################################
  CASE ("DBIF") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)// ".dbif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    ! #########################################################################
  CASE ("BIFFIT") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)//".bif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    ! #########################################################################
  CASE ("BBIF") !--------------------------------------------
    FileName=TRIM(FileName_In)//TRIM(strDate)//".bbif"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

    ! #########################################################################
  CASE ("TRMMRT") !--------------------------------------------
    FileName=TRIM(FileName_In) // "3B42RT."  &
      //TRIM(strDate(1:10))//".7R2.bin"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In) // "3B42RT."  &
        //TRIM(strDate(1:10))//".7.bin"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        FileName=TRIM(FileName_In) // "3B42RT." &
          //TRIM(strDate(1:10))// ".6A.bin"
        INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
        IF(fExist .EQV. .FALSE.)THEN
          FileName=TRIM(FileName_In) // "3B42RT."  &
            //TRIM(strDate(1:10))// ".6.bin"
          INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
          IF(fExist .EQV. .FALSE.)THEN
            FileName=TRIM(FileName_In) // "3B42RT."  &
              //TRIM(strDate(1:10))// ".bin"
            INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
            IF(fExist .EQV. .FALSE.)THEN
              bIsError=.TRUE.
              RETURN
            END IF
          END IF
        END IF
      END IF
    END IF

    ! #########################################################################
  CASE ("TRMMV6") !--------------------------------------------
    strTemp=TRIM(strDate(9:10))
    READ(strTemp,*) intTemp
    IF(intTemp<10)THEN
      WRITE(strTemp,"(I1)")intTemp
    ELSE
      WRITE(strTemp,"(I2)")intTemp
    END IF

    FileName=TRIM(FileName_In) // "3B42." &
      //TRIM(strDate(3:8))//"."//TRIM(strTemp)//".6.HDF"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      FileName=TRIM(FileName_In) // "3B42."  &
        //TRIM(strDate(3:8))//"."//TRIM(strTemp)// ".6A.HDF"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
        bIsError=.TRUE.
        RETURN
      END IF
    END IF


    ! #########################################################################
  CASE ("NMQBIN") !------XXW: Not Validate--------------------------------------
    strTemp=TRIM(strDate(9:10))
    FileName=TRIM(FileName_In) // "1HRAD.HSR."   &
      //TRIM(strDate(1:8))//"."//TRIM(strTemp)//"00"
    INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
    IF(fExist .EQV. .FALSE.)THEN
      bIsError=.TRUE.
      RETURN
    END IF

!+ygan
     !##########################################################################
    CASE ("DIAMOND4") !--------------------------------------------
      FileName=TRIM(FileName_In)//TRIM(strDate)//".000.TMP"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
         bIsError=.TRUE.
         RETURN
      END IF

      !##########################################################################
    CASE ("DIAMOND131") !-------------------------------------------- SWAN QPE QPR
      FileName=TRIM(FileName_In)//TRIM(strDate)//"0000.BIN"
      INQUIRE(FILE=TRIM(FileName), EXIST=fExist)
      IF(fExist .EQV. .FALSE.)THEN
         bIsError=.TRUE.
         RETURN
      END IF
!.ygan

    CASE DEFAULT
    WRITE(*,*) "ERROR!!! This version DO not USE this Format!"
    bIsError=.TRUE.
    RETURN
  END SELECT

  RETURN
  END SUBROUTINE InquireMatrixFile

!============================================================
!
!============================================================
  SUBROUTINE SaveStates(strDate,W0,SS0,SI0)
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE
  CHARACTER*(*) :: strDate
  LOGICAL :: bIsError
  DOUBLE PRECISION :: W0(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: SS0(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: SI0(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: dblTemp(0:g_NCols-1,0:g_NRows-1)

  ! Save the State Data
  IF(g_SaveState .EQV. .TRUE.)THEN
    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=W0
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE
    CALL WriteMatrixFile(TRIM(g_StatePath) &
      // "State_"//TRIM(strDate)// "_W0",  &
      dblTemp,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_StateFormat)

    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=SS0
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE
    CALL WriteMatrixFile(TRIM(g_StatePath) &
      // "State_"//TRIM(strDate)// "_SS0", &
      dblTemp,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_StateFormat)

    WHERE(g_Mask/=g_NoData_Value)
      dblTemp=SI0
    ELSEWHERE
      dblTemp=g_NoData_Value
    END WHERE
    CALL WriteMatrixFile(TRIM(g_StatePath) &
      // "State_"//TRIM(strDate)// "_SI0",  &
      dblTemp,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
      g_CellSize,g_NoData_Value,bIsError,g_StateFormat)
  END IF

  RETURN
  END SUBROUTINE SaveStates

!============================================================
!
!============================================================
  SUBROUTINE LoadStates(strDate,W0,SS0,SI0,bIsError)
  USE CREST_Project
  USE CREST_Basic

  IMPLICIT NONE
  CHARACTER*(*) :: strDate
  LOGICAL :: bIsError
  DOUBLE PRECISION :: W0(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: SS0(0:g_NCols-1,0:g_NRows-1)
  DOUBLE PRECISION :: SI0(0:g_NCols-1,0:g_NRows-1)

  CALL ReadMatrixFile(TRIM(g_StatePath) &
    // "State_"//TRIM(strDate)// "_W0", &
    W0,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_StateFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN ! File does not exist!
    RETURN
  END IF

  CALL ReadMatrixFile(TRIM(g_StatePath) &
    // "State_"//TRIM(strDate)// "_SS0", &
    SS0,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_StateFormat,"")

  IF(bIsError .EQV. .TRUE.)THEN ! File does not exist!
    RETURN
  END IF


  CALL ReadMatrixFile(TRIM(g_StatePath) &
    // "State_"//TRIM(strDate)// "_SI0", &
    SI0,g_NCols, g_NRows,g_XLLCorner,g_YLLCorner, &
    g_CellSize,g_NoData_Value,bIsError,g_StateFormat,"")
  IF(bIsError .EQV. .TRUE.)THEN ! File does not exist!
    RETURN
  END IF

  bIsError=.FALSE.
  RETURN
  END SUBROUTINE LoadStates

!============================================================
! THIS IS THE MAIN PROGRAM WHICH RUNS SCEUA PROGRAMS ON
! THE SIXPAR MODEL - BY QINGYUN DUAN, APRIL 1992

! Xianwu Xue modified it to a subroutine named as "SCEUSMain"
! in Oct. 12 2010
! Email: xuexianwu@ou.edu or xuexianwuqd@gmail.com
!============================================================
  SUBROUTINE GeneralSCEUA2(a_In,bl_In,bu_In,nopt_In,maxn_In, &
    kstop_In,pcento_In,iseed_In,ngs_In,idfunc_In, &
    sPath,iOutputID)
  USE CREST_Project

  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER(ijk=6000,ldt=24)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  COMMON /iopar/ ipr
  REAL*8 aa
  COMMON /ablu/ bl(50), bu(50)
  DIMENSION a(50),jseed(10)
  DIMENSION a_In(50),bl_In(50), bu_In(50)
  DATA jseed/2,3,5,7,11,13,17,19,23,29/
  CHARACTER*(*) sPath

  !
  WRITE (*,"(/,28X,A)") ' ENTER THE MAIN PROGRAM SCE-UA '

  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,"(/,28X,A)") &
    ' ENTER THE MAIN PROGRAM SCE-UA '
  WRITE(g_CREST_LogFileID,*)

  ipr=iOutputID
  idfunc=idfunc_In
  a(1:nopt_In)=a_In
  bl(1:nopt_In)=bl_In
  bu(1:nopt_In)=bu_In
  nopt=nopt_In
  maxn=maxn_In
  kstop=kstop_In
  pcento=pcento_In
  iseed=iseed_In
  ngs=ngs_In
  ! CALL XXWGetFreeFile(ipr)
  ! OPEN(UNIT=ipr,FILE=sPath // 'SCEUAOut.dat',STATUS='UNKNOWN')
  !
  ! Read input variable for SCE-UA
  CALL SCEIN(a,bl,bu,nopt,maxn,kstop,pcento,iseed, &
    ngs,npg,nps,nspl,mings,iniflg,iprint,sPath)

  !print*,npg,nps,nspl,mings,iniflg,iprint

  iniflg=1
  SELECT CASE (idfunc)
  CASE (8)
    !CALL sxpinp(nopt,a,sPath)

  CASE (1:7)
    !CALL otherf(a,50,nopt,idfunc,loc,bl,bu)
  CASE (55)

  CASE DEFAULT
    WRITE(*,*) "Something wrong in your parameters--idfunc"
  END SELECT

  IF (iseed .GT. 0) THEN
    nrun = min(iseed,10)
  ELSE
    nrun = 1
  END IF

  DO i=1, nrun
    IF (nrun .NE. 1) iseed = jseed(i)
    WRITE (*,*) '@ SCE-UA Run Number',i, &
      '    Random Seed Value',iseed

    WRITE(g_CREST_LogFileID,*)
    WRITE(g_CREST_LogFileID,*) &
      '@ SCE-UA Run Number',i, &
      '    Random Seed Value',iseed
    WRITE(g_CREST_LogFileID,*)

    CALL SCEUA(a,bl,bu,nopt,maxn,kstop,pcento,iseed, &
      ngs,npg,nps,nspl,mings,iniflg,iprint)

  END DO

  RETURN
  END SUBROUTINE GeneralSCEUA2

!============================================================
!
!============================================================
  SUBROUTINE XXWCCE(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed)
  USE CREST_Project
  USE CREST_Calib_SCEUA
  ! ALGORITHM GENERATE A NEW POINT(S) FROM A SUB-COMPLEX
  !
  ! SUB-COMPLEX VARIABLES
  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER (c1=0.8,c2=0.4)
  DIMENSION s(50,16),sf(50),bu(16),bl(16),xnstd(16)


  DIMENSION sw(16),sb(16),ce(16),snew(16)
  !
  ! EQUIVALENCE OF VARIABLES FOR READABILTY OF CODE
  n = nps
  m = nopt
  alpha = 1.0
  beta = 0.5
  !
  ! IDENTIFY THE WORST POINT wo OF THE SUB-COMPLEX s
  ! COMPUTE THE CENTROID ce OF THE REMAINING POINTS
  ! COMPUTE step, THE VECTOR BETWEEN wo AND ce
  ! IDENTIFY THE WORST FUNCTION VALUE fw
  DO j = 1, m
    sb(j) = s(1,j)
    sw(j) = s(n,j)
    ce(j) = 0.0
    DO i = 1, n-1
      ce(j) = ce(j) + s(i,j)
    END DO
    ce(j) = ce(j)/DBLE(n-1)
  END DO
  fw = sf(n)
  !
  ! COMPUTE THE NEW POINT snew
  !
  ! FIRST TRY A REFLECTION STEP
  DO j = 1, m
    snew(j) = ce(j) + alpha * (ce(j) - sw(j))
  END DO
  !
  ! CHECK IF snew SATISFIES ALL CONSTRAINTS
  CALL chkcst(nopt,snew,bl,bu,ibound)
  !
  !
  ! snew IS OUTSIDE THE BOUND,
  ! CHOOSE A POINT AT RANDOM WITHIN FEASIBLE REGION ACCORDING TO
  ! A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
  ! AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
  IF (ibound .GE. 1) CALL Getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
  !
  !
  ! COMPUTE THE FUNCTION VALUE AT snew
  fnew = XXWCRESTFunctn(nopt,snew)

  icall = icall + 1

  WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
    "Trials: ",icall,fnew,"Region Number: ", &
    g_tCalibSta(g_RegNum)%Value
  !
  ! COMPARE fnew WITH THE WORST FUNCTION VALUE fw
  !
  ! fnew IS LESS THAN fw, ACCEPT THE NEW POINT snew AND RETURN
  IF (fnew .LE. fw) GO TO 2000
  IF (icall .GE. maxn) GO TO 3000
  !
  !
  ! fnew IS GREATER THAN fw, SO TRY A CONTRACTION STEP
  DO j = 1, m
    snew(j) = ce(j) - beta * (ce(j) - sw(j))
  END DO
  !
  ! COMPUTE THE FUNCTION VALUE OF THE CONTRACTED POINT
  ! fnew = functn(nopt,snew,idfunc) !XXW
  fnew = XXWCRESTFunctn(nopt,snew)
  icall = icall + 1

  WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
    "Trials: ",icall,fnew,"Region Number: ", &
    g_tCalibSta(g_RegNum)%Value
  !
  ! COMPARE fnew TO THE WORST VALUE fw
  ! IF fnew IS LESS THAN OR EQUAL TO fw, THEN ACCEPT THE POINT AND RETURN
  IF (fnew .LE. fw) GO TO 2000
  IF (icall .GE. maxn) GO TO 3000
  !
  !
  ! IF BOTH REFLECTION AND CONTRACTION FAIL, CHOOSE ANOTHER POINT
  ! ACCORDING TO A NORMAL DISTRIBUTION WITH BEST POINT OF THE SUB-COMPLEX
  ! AS MEAN AND STANDARD DEVIATION OF THE POPULATION AS STD
1000 CALL Getpnt(nopt,2,iseed,snew,bl,bu,xnstd,sb)
  !
  ! COMPUTE THE FUNCTION VALUE AT THE RANDOM POINT
  ! fnew = functn(nopt,snew,idfunc)
  fnew = XXWCRESTFunctn(nopt,snew)
  icall = icall + 1

  WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
    "Trials: ",icall,fnew,"Region Number: ", &
    g_tCalibSta(g_RegNum)%Value
  !
  !
  ! REPLACE THE WORST POINT BY THE NEW POINT
2000 CONTINUE
  DO j = 1, m
    s(n,j) = snew(j)
  END DO
  sf(n) = fnew
3000 CONTINUE
  !
  ! END OF SUBROUTINE CCE
  RETURN
  END SUBROUTINE XXWCCE

!============================================================
!
!============================================================
  SUBROUTINE chkcst(nopt,x,bl,bu,ibound)
  !
  ! This SUBROUTINE check IF the trial point satisfies all
  ! constraints.
  !
  ! ibound - violation indicator
  ! = -1 initial value
  ! = 0  no violation
  ! = 1  violation
  ! nopt = number of optimizing variabLEs
  ! ii = the ii'th variable of the arrays x, bl, and bu
  !
  IMPLICIT REAL*8 (a-h,o-z)
  DIMENSION x(nopt),bl(nopt),bu(nopt)
  !
  ibound = -1
  !
  ! Check if explicit constraints are violated
  !
  DO ii=1, nopt
    IF (x(ii) .LT. bl(ii) .OR. x(ii) .GT. bu(ii)) GO TO 10
  END DO
  ! IF (nopt .EQ. 1) GO TO 9
  !
  ! Check IF IMPLICIT constraints are violated
  !
  ! c1 = 5.0*x(1) + 2.0*x(3) - 100.0
  ! c2 = x(3) + 50.*x(4) - 50.0
  ! c3 = 10.*x(5) - x(6) - 2.0
  ! c4 = -10.*x(5) + x(6) - 2.0
  ! IF (c1 .GT. 0.) GO TO 10
  ! IF (c2 .GT. 0.) GO TO 10
  ! IF (c3 .GT. 0.) GO TO 10
  ! IF (c4 .GT. 0.) GO TO 10
  !
  ! No constraints are violated
  !
9 ibound = 0
  RETURN
10 ibound = 1
  RETURN
  END SUBROUTINE chkcst

!============================================================
!
!============================================================
  SUBROUTINE comp(n,npt,ngs1,ngs2,npg,a,af,b,bf)
  !
  !
  ! THIS SUBROUTINE REDUCE INPUT MATRIX a(n,ngs2*npg) TO MATRIX
  ! b(n,ngs1*npg) AND VECTOR af(ngs2*npg) TO VECTOR bf(ngs1*npg)
  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER(ijk=6000,np=50)
  DIMENSION a(ijk,np),af(ijk),b(ijk,np),bf(ijk)
  DO igs=1, ngs1
    DO ipg=1, npg
      k1=(ipg-1)*ngs2 + igs
      k2=(ipg-1)*ngs1 + igs
      DO i=1, n
        b(k2,i) = a(k1,i)
      END DO
      bf(k2) = af(k1)
    END DO
  END DO
  !
  DO j=1, npt
    DO i=1, n
      a(j,i) = b(j,i)
    END DO
    af(j) = bf(j)
  END DO
  !
  ! END OF SUBROUTINE COMP
  RETURN
  END SUBROUTINE comp !
  !
  !

!============================================================
!
!============================================================
  REAL*8 FUNCTION gasdev(idum)
  !
  IMPLICIT REAL*8 (a-h,o-z)
  COMMON /gasblk/ iset
  save
  DATA iset / 0 /
  IF (iset .EQ. 0) THEN
1   v1 = (2. * ran1(idum)) - 1.
    v2 = (2. * ran1(idum)) - 1.
    r = (v1 ** 2) + (v2 ** 2)
    IF (r .GE. 1.) GOTO 1
    fac = sqrt(- ((2. * log(r)) / r))
    gset = v1 * fac
    gasdev = v2 * fac
    iset = 1
  ELSE
    gasdev = gset
    iset = 0
  END IF
  !
  ! END OF SUBROUTINE GASDEV
  RETURN
  END FUNCTION gasdev

!============================================================
!
!============================================================
  SUBROUTINE Getpnt(nopt,idist,iseed,x,bl,bu,std,xi)
  !
  !
  IMPLICIT REAL*8 (a-h,o-z)
  DIMENSION x(16),bl(16),bu(16),std(16),xi(16)
  !
1 DO j=1, nopt

2   IF(bu(j)-bl(j)/=0.0)THEN
      IF (idist .EQ. 1) rand = ran1(iseed)
      IF (idist .EQ. 2) rand = gasdev(iseed)

      x(j) = xi(j) + std(j) * rand * (bu(j) - bl(j))
    ELSE
      x(j)=bl(j)
    END IF
    ! IF(rand.LT.0.0000001)WRITE(*,*)rand
    !
    ! Check explicit constraints
    !
    CALL chkcst(1,x(j),bl(j),bu(j),ibound)
    IF (ibound .GE. 1) GO TO 2
  END DO
  !
  ! Check IMPLICIT constraints
  !
  CALL chkcst(nopt,x,bl,bu,ibound)
  IF (ibound .GE. 1) GO TO 1
  !
  RETURN
  END SUBROUTINE Getpnt

!============================================================
!
!============================================================
  SUBROUTINE indexx(n, arrin, indx)

  ! THIS SUBROUTINE IS FROM "NUMERICAL RECIPES" BY PRESS ET AL.
  IMPLICIT REAL*8 (a-h,o-z)
  INTEGER n
  DOUBLE PRECISION arrin
  INTEGER index
  DIMENSION arrin(n), indx(n)

  INTEGER i,j,l,ir
  !
  DO 11 j = 1, n
    indx(j) = j
11 CONTINUE
  l = (n / 2) + 1
  ir = n
10 CONTINUE
  IF (l .GT. 1) THEN
    l = l - 1
    indxt = indx(l)
    q = arrin(indxt)
  ELSE
    indxt = indx(ir)
    q = arrin(indxt)
    indx(ir) = indx(1)
    ir = ir - 1
    IF (ir .EQ. 1) THEN
      indx(1) = indxt
      RETURN
    END IF
  END IF
  i = l
  j = l + l
20 IF (j .LE. ir) THEN
    IF (j .LT. ir) THEN
      IF (arrin(indx(j)) .LT. arrin(indx(j + 1))) j = j + 1
    END IF
    IF (q .LT. arrin(indx(j))) THEN
      indx(i) = indx(j)
      i = j
      j = j + j
    ELSE
      j = ir + 1
    END IF
    GOTO 20
  END IF
  indx(i) = indxt
  GOTO 10
  !
  ! END OF SUBROUTINE INDEXX
  END SUBROUTINE indexx

!============================================================
!
!============================================================
  SUBROUTINE parstt(npt,nopt,x,xnstd,bound,gnrng,ipcnvg)
  !
  ! SUBROUTINE CHECKING FOR PARAMETER CONVERGENCE
  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER(ijk=6000,np=50)
  DIMENSION x(ijk,np),xmax(np),xmin(np)
  DIMENSION xmean(np),xnstd(np),bound(np)
  PARAMETER (deLTa = 1.0d-20,peps=1.0d-3)
  !
  ! COMPUTE MAXIMUM, MINIMUM AND STANDARD DEVIATION OF PARAMETER VALUES
  gsum = 0.d0
  DO k = 1, nopt
    xmax(k) = -1.0d+20
    xmin(k) = 1.0d+20
    xsum1 = 0.d0
    xsum2 = 0.d0
    DO i = 1, npt
      xmax(k) = dmax1(x(i,k), xmax(k))
      xmin(k) = dmin1(x(i,k), xmin(k))
      xsum1 = xsum1 + x(i,k)
      xsum2 = xsum2 + x(i,k)*x(i,k)
    END DO
    xmean(k) = xsum1 / DBLE(npt)
    xnstd(k) = (xsum2 / DBLE(npt) - xmean(k)*xmean(k))
    IF (xnstd(k) .LE. deLTa) xnstd(k) = deLTa
    xnstd(k) = DSQRT(xnstd(k))
    xnstd(k) = xnstd(k) / bound(k)
    gsum = gsum + DLOG( deLTa + (xmax(k)-xmin(k))/bound(k) )
  END DO
  gnrng = DEXP(gsum/DBLE(nopt))
  !
  ! CHECK IF NORMALIZED STANDARD DEVIATION OF PARAMETER IS <= eps
  ipcnvg = 0
  IF (gnrng .LE. peps) THEN
    ipcnvg = 1
  END IF
  !
  ! END OF SUBROUTINE PARSTT
  RETURN
  END SUBROUTINE parstt
  !
  !
  !
!============================================================
!
!============================================================
  REAL*8 FUNCTION ran1(idum)
  !
  IMPLICIT REAL*8 (a-h,o-z)
  DIMENSION r(97)
  PARAMETER (m1 = 259200, ia1 = 7141, ic1 = 54773, &
    rm1 =3.8580247e-6)
  PARAMETER (m2 = 134456, ia2 = 8121, ic2 = 28411, &
    rm2 =7.4373773e-6)
  PARAMETER (m3 = 243000, ia3 = 4561, ic3 = 51349)
  save
  DATA iff / 0 /
  IF ((idum .LT. 0) .OR. (iff .EQ. 0)) THEN
    iff = 1
    ix1 = MOD(ic1 - idum,m1)
    ix1 = MOD((ia1 * ix1) + ic1,m1)
    ix2 = MOD(ix1,m2)
    ix1 = MOD((ia1 * ix1) + ic1,m1)
    ix3 = MOD(ix1,m3)
    DO 11 j = 1, 97
      ix1 = MOD((ia1 * ix1) + ic1,m1)
      ix2 = MOD((ia2 * ix2) + ic2,m2)
      r(j) = (DBLE(ix1) + (DBLE(ix2) * rm2)) * rm1
11  CONTINUE
    idum = 1
  END IF

  ix1 = MOD((ia1 * ix1) + ic1,m1)
  ix2 = MOD((ia2 * ix2) + ic2,m2)
  ix3 = MOD((ia3 * ix3) + ic3,m3)
  j = 1 + ((97 * ix3) / m3)
  !IF ((j .GT. 97) .OR. (j .LT. 1)) paUSE
  ran1 = r(j)
  r(j) = (DBLE(ix1) + (DBLE(ix2) * rm2)) * rm1
  !
  ! END OF SUBROUTINE RAN1
  RETURN
  END FUNCTION ran1

!============================================================
! scein.for
!============================================================
  SUBROUTINE SCEIN(a,bl,bu,  &
    nopt,maxn,kstop,pcento,iseed,  &
    ngs,npg,nps,nspl,mings,iniflg,iprint, &
    sPath)
  !
  ! THIS SUBROUTINE READS AND PRINTS THE INPUT VARIABLES FOR
  ! SHUFFLED COMPLEX EVOLUTION METHOD FOR GLOBAL OPTIMIZATION
  ! -- Version 2.2
  !
  ! WRITTEN BY QINGYUN DUAN - UNIVERSITY OF ARIZONA, APRIL 1992
  !
  !
  IMPLICIT REAL*8 (a-h,o-z)
  ! ############Parameters

  DOUBLE PRECISION a,bl,bu
  DIMENSION a(50),bl(50),bu(50)
  CHARACTER*(*)sPath
  INTEGER nopt,maxn,kstop,iseed
  INTEGER ngs,npg,nps,nspl,mings,iniflg,iprint
  DOUBLE PRECISION pcento
  ! ############Global Variables
  INTEGER ipr !Block iopar

  CHARACTER*10 pcntrl,deflt,usrsp
  CHARACTER*4 reduc,initl,ysflg,noflg,xname(50)
  PARAMETER(ijk=6000,ldt=24)

  ! ############Local Variables

  INTEGER in !Block iopar

  ! COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  ! COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  ! COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  ! COMMON /routing/ckm,cxm,nm
  COMMON /iopar/ ipr
  ! COMMON /idid/iduh,area

  DATA deflt/' DEFAULT  '/
  DATA usrsp/'USER SPEC.'/
  DATA ysflg/'YES '/
  DATA noflg/'NO  '/
  DATA xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7', &
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16', &
    ' X17',' X18',' X19',' X20',' X21',' X22',' X23',' X24',' X25', &
    ' X26',' X27',' X28',' X29',' X30',' X31',' X32',' X33',' X34', &
    ' X35',' X36',' X37',' X38',' X39',' X40',' X41',' X42',' X43', &
    ' X44',' X45',' X46',' X47',' X48',' X49',' X50'/
  !
  WRITE (*,"(30X,A)") ' ENTER THE SCEIN SUBROUTINE'
  !
  !
  ! INITIALIZE I/O VARIABLES
  in = 12
  IF(idfunc/=55)THEN
    OPEN(UNIT=in,FILE=sPath // 'scein.dat',STATUS='OLD')
  END IF

  ierror = 0
  iwarn = 0
  WRITE(ipr,700)
700 FORMAT(10X,'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION', &
    /,10x,46(1H=))
  !
  !
  ! READ THE SCE CONTROL PARAMETERS
  ideflt = 0
  IF(idfunc/=55)THEN
    READ(in,*)idfunc
    READ(in,*) nopt,maxn,kstop,pcento,ngs,iseed,ideflt
  END IF
  IF (iseed .EQ. 0) iseed = 1969
  !
  !
  ! IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS
  IF (ideflt .EQ. 1) THEN
    IF(idfunc/=55)THEN
      READ(in,*) npg,nps,nspl,mings,iniflg,iprint
    END IF
    pcntrl = usrsp
  ELSE
    ! READ(in,*)
    pcntrl = deflt
  END IF

  IF(idfunc.NE.8) GO TO 555

  DO 654 i=1,nopt
    IF(idfunc/=55)THEN
      READ(in,*) a(i),bl(i),bu(i)
    END IF
    loc(i)=i
654 CONTINUE

555 CONTINUE
  !
  ! IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO
  ! THE DEFAULT VALUES
  IF (ideflt .EQ. 0) THEN
    npg = 2*nopt + 1
    nps = nopt + 1
    nspl = npg
    mings = ngs
    iniflg = 0
    iprint = 0
  END IF
  !
  !
  ! CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
  IF (ngs .LT. 1 .OR. ngs .GE. 1320) THEN
    WRITE(ipr,900) ngs
900 FORMAT(//,1X,'**ERROR** NUMBER OF COMPLEXES IN INITIAL ', &
      ' POPULATION ',I5,' IS NOT A VALID CHOICE')
    ierror = ierror + 1
  END IF
  !
  IF (kstop .LT. 0 .OR. kstop .GE. 20) THEN
    WRITE(ipr,901) kstop
901 FORMAT(//,1X,'**WARNING** THE NUMBER OF SHUFFLING LOOPS IN', &
      ' WHICH THE CRITERION VALUE MUST CHANGE ',/,13X,'SHOULD BE', &
      ' GREATER THAN 0 AND LESS THAN 10.  ','kstop = ',i2, &
      ' WAS SPECIFIED.'/,13X,'BUT kstop = 5 WILL BE USED INSTEAD.')
    iwarn = iwarn + 1
    kstop=5
  END IF
  !
  IF (mings .LT. 1 .OR. mings .GT. ngs) THEN
    WRITE(ipr,902) mings
902 FORMAT(//,1X,'**WARNING** THE MINIMUM NUMBER OF COMPLEXES ', &
      I2,' IS NOT A VALID CHOICE. SET IT TO DEFAULT')
    iwarn = iwarn + 1
    mings = ngs
  END IF
  !
  IF (npg .LT. 2 .OR. npg .GT. 1320/MAX(ngs,1)) THEN
    WRITE(ipr,903) npg
903 FORMAT(//,1X,'**WARNING** THE NUMBER OF POINTS IN A COMPLEX ', &
      I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
    iwarn = iwarn + 1
    npg = 2*nopt+1
  END IF
  !
  IF (nps.LT.2 .OR. nps.GT.npg .OR. nps.GT.50) THEN
    WRITE(ipr,904) nps
904 FORMAT(//,1X,'**WARNING** THE NUMBER OF POINTS IN A SUB-', &
      'COMPLEX ',I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
    iwarn = iwarn + 1
    nps = nopt + 1
  END IF
  !
  IF (nspl .LT. 1) THEN
    WRITE(ipr,905) nspl
905 FORMAT(//,1X,'**WARNING** THE NUMBER OF EVOLUTION STEPS ', &
      'TAKEN IN EACH COMPLEX BEFORE SHUFFLING ',I4,/,13X, &
      'IS NOT A VALID CHOICE, SET IT TO DEFAULT')
    iwarn = iwarn + 1
    nspl = npg
  END IF
  !
  ! COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
  npt = ngs * npg
  !
  IF (npt .GT. 1320) THEN
    WRITE(ipr,906) npt
906 FORMAT(//,1X,'**WARNING** THE NUMBER OF POINTS IN INITIAL ', &
      'POPULATION ',I5,' EXCEED THE POPULATION LIMIT,',/,13X, &
      'SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS')
    iwarn = iwarn + 1
    ngs = 2
    npg = 2*nopt + 1
    nps = nopt + 1
    nspl = npg
  END IF
  !
  ! PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
  IF (ierror .GE. 1) WRITE(ipr,907) ierror
907 FORMAT(//,1X,'*** TOTAL NUMBER OF ERROR MESSAGES IS ',I2)
  !
  IF (iwarn .GE. 1) WRITE(ipr,908) iwarn
908 FORMAT(//,1X,'*** TOTAL NUMBER OF WARNING MESSAGES IS ',I2)
  !
  IF (mings .LT. ngs) THEN
    reduc = ysflg
  ELSE
    reduc = noflg
  END IF
  !
  IF (iniflg .NE. 0) THEN
    initl = ysflg
  ELSE
    initl = noflg
  END IF
  !
  !
  ! PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
104 WRITE(ipr,910)
910 FORMAT(//,2X,'SCE CONTROL',5X,'MAX TRIALS',5X, &
    'REQUIRED IMPROVEMENT',5X,'RANDOM',/,3X,'PARAMETER',8X, &
    'ALLOWED',6X,'PERCENT',4X,'NO. LOOPS',6X,'SEED',/, &
    2X,11(1H-),5X,10(1H-),5X,7(1H-),4X,9(1H-),5X,6(1H-))
  !
  pcenta=pcento*100.
  WRITE(ipr,912) pcntrl,maxn,pcenta,kstop,iseed
912 FORMAT(3X,A10,7X,I5,10X,F5.1,9X,I2,9X,I5)
  WRITE(ipr,914) ngs,npg,npt,nps,nspl
914 FORMAT(//,18X,'SCE ALGORITHM CONTROL PARAMETERS',/,18X,32(1H=), &
    //,2X,'NUMBER OF',5X,'POINTS PER',5X,'POINTS IN',6X,'POINTS PER', &
    4X,'EVOL. STEPS',/,2X,'COMPLEXES',6X,'COMPLEX',6X,'INI. POPUL.', &
    5X,'SUB-COMPLX',4X,'PER COMPLEX',/,2X,9(1H-),5X,10(1H-),4X, &
    11(1H-),5X,10(1H-),4X,11(1H-),5X,/,2X,5(I5,10X))
  WRITE(ipr,915) reduc,mings,initl
915 FORMAT(//,15X,'COMPLX NO.',5X,'MIN COMPLEX',5X,'INI. POINT',/, &
    15X,'REDUCTION',6X,'NO. ALLOWED',6X,'INCLUDED',/, &
    15X,10(1H-),5X,11(1H-),5X,10(1H-),/,18X,A4,6X,I8,13X,A4)
  WRITE(ipr,916)
916 FORMAT(//,8X,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/, &
    8X,45(1H=),//,2X,'PARAMETER',5X,'INITIAL VALUE',5X &
    'LOWER BOUND',5X,'UPPER BOUND',/,2X,9(1H-),5X,13(1H-),5X, &
    11(1H-),5X,11(1H-))
  DO 920 i = 1, nopt
    WRITE(ipr,918) xname(i),a(i),bl(i),bu(i)
918 FORMAT(4X,A4,4X,3(6X,F10.3))
920 CONTINUE
  IF (ierror .GE. 1) THEN
    WRITE(ipr,922)
922 FORMAT(//,'*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE', &
      ' OF INPUT DATA ERROR ***')
    stop
  END IF

  RETURN
  END SUBROUTINE SCEIN


!============================================================
!  sceua.for
!============================================================
  SUBROUTINE SCEUA(a,bl,bu,nopt,maxn,kstop,pcento,iseed, &
    ngs,npg,nps,nspl,mings,iniflg,iprint)
  USE CREST_Project
  USE CREST_Calib_SCEUA

  !
  IMPLICIT REAL*8 (a-h,o-z)
  !
  ! ARRAYS FROM THE INPUT DATA
  DIMENSION a(16),bl(16),bu(16)
  !
  ! LOCAL ARRAYS
  PARAMETER(ijk=6000,np=50)
  DIMENSION x(ijk,np),xx(np),bestx(np),worstx(np),xf(ijk)
  DIMENSION s(50,np),sf(50),lcs(50),cx(ijk,np),cf(ijk)
  DIMENSION xnstd(np),bound(np),criter(20),unit(np)

  DOUBLE PRECISION :: XXWBestValue !XXW
  LOGICAL :: fExist
  !
  COMMON /iopar/ ipr
  !
  CHARACTER*4 xname(16)
  DATA xname /'  X1','  X2','  X3','  X4','  X5','  X6','  X7', &
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16'/
  !
  WRITE (*,"(30X,A,/)") ' ENTER THE SCEUA SUBROUTINE '
  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,"(30X,A,/)") &
    ' ENTER THE SCEUA SUBROUTINE '
  WRITE(g_CREST_LogFileID,*)
  !
  ! INITIALIZE VARIABLES
  nloop = 0
  loop = 0
  igs = 0
  nopt1 = 8
  IF (nopt.LT.8) nopt1 = nopt
  nopt2 = 12
  IF (nopt.LT.12) nopt2 = nopt
  !
  ! INITIALIZE RANDOM SEED TO A NEGATIVE INTEGER
  iseed1 = -ABS(iseed)
  !
  ! COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPUALTION
  npt = ngs * npg
  ngs1 = ngs
  npt1 = npt
  !
  WRITE(ipr,400)
  WRITE (*,*) ' ***  Evolution Loop Number ',nloop
  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,*) &
    ' ***  Evolution Loop Number ',nloop
  !
  ! COMPUTE THE BOUND FOR PARAMETERS BEING OPTIMIZED
  DO j = 1, nopt
    bound(j) = bu(j) - bl(j)
    unit(j) = 1.
    !print*,bu(j),bl(j),bound(j)
  END DO
  !
  ! COMPUTE THE FUNCTION VALUE OF THE INITIAL POINT
  !print*,nopt,a
  ! fa = functn(nopt,a,idfunc) !XXW
  fa = XXWCRESTFunctn(nopt,a)
  WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
    "Trials: ",icall,fa,"Region Number: ", &
    g_tCalibSta(g_RegNum)%Value

  !To stop for anytime when "XXW.CREST.STOP" EXISTs
  INQUIRE(FILE=TRIM(g_ProjectWS)//TRIM("XXW.CREST.STOP"), &
    EXIST=fExist)
  IF(fExist .EQV. .TRUE.)THEN
    GO TO 9500
  END IF
  !
  ! PRINT THE INITIAL POINT AND ITS CRITERION VALUE
  WRITE(ipr,500)
  WRITE(ipr,510) (xname(j),j=1,nopt2)
  WRITE(ipr,520) fa,(a(j),j=1,nopt2)
  IF (nopt.LT.13) GO TO 101
  WRITE(ipr,530) (xname(j),j=13,nopt)
  WRITE(ipr,540) (a(j),j=13,nopt)
101 CONTINUE
  !
  ! GENERATE AN INITIAL SET OF npt1 POINTS IN THE PARAMETER SPACE
  ! IF iniflg IS EQUAL TO 1, SET x(1,.) TO INITIAL POINT a(.)
  IF (iniflg .EQ. 1) THEN
    DO j = 1, nopt
      x(1,j) = a(j)
    END DO
    xf(1) = fa
  !
  ! ELSE, GENERATE A POINT RANDOMLY AND SET IT EQUAL TO x(1,.)
  ELSE
    CALL Getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)
    DO j=1, nopt
      x(1,j) = xx(j)
    END DO
    ! xf(1) = functn(nopt,xx,idfunc) !XXW
    xf(1) = XXWCRESTFunctn(nopt,xx)

    icall = 1
    WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
      "Trials: ",icall,xf(1),"Region Number: ", &
      g_tCalibSta(g_RegNum)%Value

  END IF
  icall = 1
  IF (icall .GE. maxn) GO TO 9000
  !
  ! GENERATE npt1-1 RANDOM POINTS DISTRIBUTED UNIFORMLY IN THE PARAMETER
  ! SPACE, AND COMPUTE THE CORRESPONDING FUNCTION VALUES
  DO i = 2, npt1
    CALL Getpnt(nopt,1,iseed1,xx,bl,bu,unit,bl)
    DO j = 1, nopt
      x(i,j) = xx(j)
    END DO
    ! xf(i) = functn(nopt,xx,idfunc)
    xf(i) = XXWCRESTFunctn(nopt,xx)

    icall = icall + 1

    WRITE(g_CREST_LogFileID,"(4X,A8,I5,F20.15,3X,A15,I3)") &
      "Trials: ",icall,xf(i),"Region Number: ", &
      g_tCalibSta(g_RegNum)%Value

    IF (icall .GE. maxn) THEN
      nopt1 = i
      GO TO 45
    END IF
  END DO
  !
  ! ARRANGE THE POINTS IN ORDER OF INCREASING FUNCTION VALUE
45 CALL sort(npt1,nopt,x,xf)
  !
  ! RECORD THE BEST AND WORST POINTS
  DO j = 1, nopt
    bestx(j) = x(1,j)
    worstx(j) = x(npt1,j)
  END DO
  bestf = xf(1)
  worstf = xf(npt1)
  !
  ! COMPUTE THE PARAMETER RANGE FOR THE INITIAL POPULATION
  CALL parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)
  !
  ! PRINT THE RESULTS FOR THE INITIAL POPULATION
  WRITE(ipr,600)
  WRITE(ipr,610) (xname(j),j=1,nopt1)
  IF (nopt .LT. 9) GO TO 201
  WRITE(ipr,620) (xname(j),j=9,nopt)
201 CONTINUE
  WRITE(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng, &
    (bestx(j),j=1,nopt1)
  IF (nopt .LT. 9) GO TO 301
  WRITE(ipr,640) (bestx(j),j=9,nopt)
301 CONTINUE
  ! iprint=1
  IF (iprint .EQ. 1) THEN
    WRITE(ipr,650) nloop
    DO i = 1, npt1
      WRITE(ipr,660) xf(i),(x(i,j),j=1,nopt1)
      IF (nopt .LT. 9) GO TO 401
      WRITE(ipr,640) (x(i,j),j=9,nopt)
401 END DO
  END IF
  !
  IF (icall .GE. maxn) GO TO 9000
  IF (ipcnvg .EQ. 1) GO TO 9200
  !
  ! BEGIN THE MAIN LOOP ----------------
1000 CONTINUE


  INQUIRE(FILE=TRIM(g_ProjectWS)//TRIM("XXW.CREST.STOP"),  &
    EXIST=fExist)
  IF(fExist .EQV. .TRUE.)THEN
    GO TO 9500
  END IF
  nloop = nloop + 1
  !
  WRITE (*,*) ' ***  Evolution Loop Number ',nloop
  WRITE(g_CREST_LogFileID,*) &
    ' ***  Evolution Loop Number ',nloop
  !
  ! BEGIN LOOP ON COMPLEXES
  DO 3000 igs = 1, ngs1
    !
    ! ASSIGN POINTS INTO COMPLEXES
    DO k1 = 1, npg
      k2 = (k1-1) * ngs1 + igs
      DO j = 1, nopt
        cx(k1,j) = x(k2,j)
      END DO
      cf(k1) = xf(k2)
    END DO
    !
    ! BEGIN INNER LOOP - RANDOM SELECTION OF SUB-COMPLEXES ---------------
    DO 2000 loop = 1, nspl
      !
      ! CHOOSE A SUB-COMPLEX (nps points) ACCORDING TO A LINEAR
      ! PROBABILITY DISTRIBUTION
      IF (nps .EQ. npg) THEN
        DO k = 1, nps
          lcs(k) = k !lcs(.) = indices locating position of s(.,.)
        END DO
        GO TO 85
      END IF
      !
      rand = ran1(iseed1)


      lcs(1) = 1 + DINT(npg + 0.5 - DSQRT( (npg+.5)**2 - &
        npg * (npg+1) * rand ))
      DO k = 2, nps
60      rand = ran1(iseed1)
        lpos = 1 + DINT(npg + 0.5 - DSQRT((npg+.5)**2 - &
          npg * (npg+1) * rand ))
        DO k1 = 1, k-1
          IF (lpos .EQ. lcs(k1)) GO TO 60
        END DO
        lcs(k) = lpos
      END DO
      !
      ! ARRANGE THE SUB-COMPLEX IN ORDER OF INCEASING FUNCTION VALUE
      CALL sort1(nps,lcs)
      !
      ! CREATE THE SUB-COMPLEX ARRAYS
85    DO k = 1, nps
        DO j = 1, nopt
          s(k,j) = cx(lcs(k),j)
        END DO
        sf(k) = cf(lcs(k))
      END DO
      !
      ! USE THE SUB-COMPLEX TO GENERATE NEW POINT(S)
      CALL XXWCCE(nopt,nps,s,sf,bl,bu,xnstd,icall,maxn,iseed1)
      !
      ! IF THE SUB-COMPLEX IS ACCEPTED, REPLACE THE NEW SUB-COMPLEX
      ! INTO THE COMPLEX
      DO k = 1, nps
        DO j = 1, nopt
          cx(lcs(k),j) = s(k,j)
        END DO
        cf(lcs(k)) = sf(k)
      END DO
      !
      ! SORT THE POINTS
      CALL sort(npg,nopt,cx,cf)
      !
      ! IF MAXIMUM NUMBER OF RUNS EXCEEDED, BREAK OUT OF THE LOOP
      IF (icall .GE. maxn) GO TO 2222
      !
      ! END OF INNER LOOP ------------
2000 CONTINUE
2222 CONTINUE
    !
    ! REPLACE THE NEW COMPLEX INTO ORIGINAL ARRAY x(.,.)
    DO k1 = 1, npg
      k2 = (k1-1) * ngs1 + igs
      DO j = 1, nopt
        x(k2,j) = cx(k1,j)
      END DO
      xf(k2) = cf(k1)
    END DO
    IF (icall .GE. maxn) GO TO 3333
    !
    ! END LOOP ON COMPLEXES
3000 CONTINUE
  !
  ! RE-SORT THE POINTS
3333 CALL sort(npt1,nopt,x,xf)
  !
  ! RECORD THE BEST AND WORST POINTS
  DO j = 1, nopt
    bestx(j) = x(1,j)
    worstx(j) = x(npt1,j)
  END DO
  bestf = xf(1)
  worstf = xf(npt1)
  !
  ! TEST THE POPULATION FOR PARAMETER CONVERGENCE
  CALL parstt(npt1,nopt,x,xnstd,bound,gnrng,ipcnvg)
  ! WRITE(*,*)"   Nloop:",nloop
  !
  ! PRINT THE RESULTS FOR CURRENT POPULATION
  IF (MOD(nloop,5) .NE. 0) GO TO 501
  WRITE(ipr,610) (xname(j),j=1,nopt1)
  IF (nopt .LT. 9) GO TO 501
  WRITE(ipr,620) (xname(j),j=9,nopt)
501 CONTINUE
  WRITE(ipr,630) nloop,icall,ngs1,bestf,worstf,gnrng, &
    (bestx(j),j=1,nopt1)
  IF (nopt.LT.9) GO TO 601
  WRITE(ipr,640) (bestx(j),j=9,nopt)
601 CONTINUE
  IF (iprint .EQ. 1) THEN
    WRITE(ipr,650) nloop
    DO i = 1, npt1
      WRITE(ipr,660) xf(i),(x(i,j),j=1,nopt1)
      IF (nopt .LT. 9) GO TO 701
      WRITE(ipr,640) (x(i,j),j=9,nopt)
701 END DO
  END IF
  !
  ! TEST IF MAXIMUM NUMBER OF FUNCTION EVALUATIONS EXCEEDED
  IF (icall .GE. maxn) GO TO 9000
  !
  ! COMPUTE THE COUNT ON SUCCESSIVE LOOPS W/O FUNCTION IMPROVEMENT
  criter(20) = bestf
  IF (nloop .LT. (kstop+1)) GO TO 132
  denomi = DABS(criter(20-kstop) + criter(20)) / 2.
  timeou = DABS(criter(20-kstop) - criter(20)) / denomi
  IF (timeou .LT. pcento) GO TO 9100
132 CONTINUE
  DO l = 1, 19
    criter(l) = criter(l+1)
  END DO
  !
  ! IF POPULATION IS CONVERGED INTO A SUFFICIENTLY SMALL SPACE
  IF (ipcnvg .EQ. 1) GO TO 9200
  !
  ! NONE OF THE STOPPING CRITERIA IS SATISFIED, CONTINUE SEARCH
  !
  ! CHECK FOR COMPLEX NUMBER REDUCTION
  IF (ngs1 .GT. mings) THEN
    ngs2 = ngs1
    ngs1 = ngs1 - 1
    npt1 = ngs1 * npg
    CALL comp(nopt,npt1,ngs1,ngs2,npg,x,xf,cx,cf)
  END IF
  !
  ! END OF MAIN LOOP -----------
  GO TO 1000
  !
  ! SEARCH TERMINATED
9000 CONTINUE
  WRITE(ipr,800) maxn,loop,igs,nloop
  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,800) maxn,loop,igs,nloop
  GO TO 9999
9100 CONTINUE
  WRITE(ipr,810) pcento*100.,kstop
  WRITE(g_CREST_LogFileID,*)
  WRITE(g_CREST_LogFileID,810) pcento*100.,kstop
  GO TO 9999
9200 CONTINUE
  WRITE(g_CREST_LogFileID,*)
  WRITE(ipr,820) gnrng*100.
  WRITE(g_CREST_LogFileID,820) gnrng*100.
  GO TO 9999
9500 CONTINUE
  WRITE(g_CREST_LogFileID,*)
  WRITE(ipr,850) pcento*100.,kstop !XXW
  WRITE(g_CREST_LogFileID,850) pcento*100.,kstop !XXW
9999 CONTINUE
  !
  ! PRINT THE FINAL PARAMETER ESTIMATE AND ITS FUNCTION VALUE
  WRITE(ipr,830)
  WRITE(ipr,510) (xname(j),j=1,nopt2)
  WRITE(ipr,520) bestf,(bestx(j),j=1,nopt2)
  IF (nopt.LT.13) GO TO 801
  WRITE(ipr,530) (xname(j),j=13,nopt)
  WRITE(ipr,540) (bestx(j),j=13,nopt)

801 CONTINUE

  XXWBestValue = XXWCRESTFunctn(nopt,bestx)

  WRITE(g_CREST_LogFileID,"(4X,A12,F20.15,3X,A15,I3)") &
    "Best Value: ",XXWBestValue,"Region Number: ", &
    g_tCalibSta(g_RegNum)%Value
  !
  ! END OF SUBROUTINE SCEUA
  RETURN
400 FORMAT(//,2X,50(1H=),/,2X,'ENTER THE SHUFFLED COMPLEX EVOLUTION', &
    ' GLOBAL SEARCH',/,2X,50(1H=))

500 FORMAT(//,'*** PRINT THE INITIAL POINT AND ITS CRITERION ',   &
    'VALUE ***')
510 FORMAT(/,' CRITERION',12(6X,A4),/1X,60(1H-))
520 FORMAT(G10.3,12F10.3)
530 FORMAT(10X,12(6X,A4))
540 FORMAT(10X,12F10.3)

600 FORMAT(//,1X,'*** PRINT THE RESULTS OF THE SCE SEARCH ***')
610 FORMAT(/,1X,'LOOP',1X,'TRIALS',1X,'COMPLXS',2X,'BEST F',3X, &
    'WORST F',3X,'PAR RNG',1X,8(6X,A4))
620 FORMAT(49X,8(6X,A4))
630 FORMAT(I5,1X,I5,3X,I5,3G10.3,8(F10.3))
640 FORMAT(49X,8(F10.3))
650 FORMAT(/,1X,'POPULATION AT LOOP ',I3,/,1X,22(1H-))
660 FORMAT(15X,G10.3,20X,8(F10.3))

800 FORMAT(//,1X,'*** OPTIMIZATION SEARCH TERMINATED BECAUSE THE',  &
    ' LIMIT ON THE MAXIMUM',/,5X,'NUMBER OF TRIALS ',I5, &
    ' EXCEEDED.  SEARCH WAS STOPPED AT',/,5X,'SUB-COMPLEX ', &
    I3,' OF COMPLEX ',I3,' IN SHUFFLING LOOP ',I3,' ***')
810 FORMAT(//,1X,'*** OPTIMIZATION TERMINATED BECAUSE THE CRITERION',  &
    ' VALUE HAS NOT CHANGED ',/,5X,F5.2,' PERCENT IN',I3, &
    ' SHUFFLING LOOPS ***')
820 FORMAT(//,1X,'*** OPTIMIZATION TERMINATED BECAUSE THE POPULATION',  &
    ' HAS CONVERGED INTO ',/,4X,F5.2,' PERCENT OF THE', &
    ' FEASIBLE SPACE ***')
830 FORMAT(//,'*** PRINT THE FINAL PARAMETER ESTIMATE AND ITS',  &
    ' CRITERION VALUE ***')

850 FORMAT(//,1X,'*** OPTIMIZATION TERMINATED BECAUSE DR. XIANWU XUE',  &
    ' WANT TO STOP IT! ',/,5X,F5.2,' PERCENT IN',I3, &
    ' SHUFFLING LOOPS ***')

  END SUBROUTINE SCEUA

!============================================================
!
!============================================================
  SUBROUTINE sort(n,m,rb,ra)

  !
  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER(ijk=6000)
  DIMENSION ra(ijk),rb(ijk,16),wk(ijk,16),iwk(ijk)
  !
  CALL indexx(n, ra, iwk)
  DO 11 i = 1, n
    wk(i,1) = ra(i)
11 CONTINUE
  DO 12 i = 1, n
    ra(i) = wk(iwk(i),1)
12 CONTINUE
  DO 14 j = 1, m
    DO 13 i = 1, n
      wk(i,j) = rb(i,j)
13 CONTINUE
14 CONTINUE
  DO 16 j = 1, m
    DO 15 i = 1, n
      rb(i,j) = wk(iwk(i),j)
15 CONTINUE
16 CONTINUE
  !
  ! END OF SUBROUTINE SORT
  RETURN
  END SUBROUTINE sort

!============================================================
!
!============================================================
  SUBROUTINE sort1(n,ra)
  !
  IMPLICIT REAL*8 (a-h,o-z)
  DIMENSION ra(n)
  !
  INTEGER ra, rra
  !
  l = (n / 2) + 1
  ir = n
10 CONTINUE
  IF (l .GT. 1) THEN
    l = l - 1
    rra = ra(l)
  ELSE
    rra = ra(ir)
    ra(ir) = ra(1)
    ir = ir - 1
    IF (ir .EQ. 1) THEN
      ra(1) = rra
      RETURN
    END IF
  END IF
  i = l
  j = l + l
20 IF (j .LE. ir) THEN
    IF (j .LT. ir) THEN
      IF (ra(j) .LT. ra(j + 1)) j = j + 1
    END IF
    IF (rra .LT. ra(j)) THEN
      ra(i) = ra(j)
      i = j
      j = j + j
    ELSE
      j = ir + 1
    END IF
    GOTO 20
  END IF
  ra(i) = rra
  GOTO 10
  !
  ! END OF SUBROUTINE SORT1
  END SUBROUTINE sort1

!============================================================
!
!============================================================
  FUNCTION functn(nopt,x,idfunc)
  IMPLICIT NONE
  DOUBLE PRECISION functn,functOthers,functnsc
  INTEGER nopt
  DOUBLE PRECISION x
  DIMENSION x(nopt)

  INTEGER idfunc

  IF ((idfunc .GE. 1).AND.(idfunc .LE. 7)) THEN
    functn = functOthers(nopt,x)
  END IF
  ! This is the SAC model
  IF (idfunc .EQ. 8) THEN
    functn = functnsc(nopt,x)
  END IF

  ! END OF FUNCTION FUNCTN
  RETURN
  END FUNCTION functn

!============================================================
!
!============================================================
  FUNCTION functOthers(nopt,x)
  IMPLICIT REAL*8 (a-h,o-z)

  PARAMETER(ijk=6000,ldt=24,bigm = 1.0d10)

  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  ! COMMON /iopar/ ipr
  DIMENSION x(nopt), par(50)
  COMMON /ablu/ bl(50), bu(50)


  DIMENSION a1(10,4), c1(10), a2(4,6), c2(4), p2(4,6), a3(4,3), c3(4), p3(4,3)
  ! Data for Skekel FUNCTION coefficients (n=4, m=10)
  DATA a1 /4.,1.,8.,6.,3.,2.,5.,8.,6.,7., 4.,1.,8.,6.,7.,9.,5., &
    1.,2.,3.6, 4.,1.,8.,6.,3.,2.,3.,8.,6.,7., &
    4.,1.,8.,6.,7.,9.,3.,1.,2.,3.6/
  DATA c1 /.1,.2,.2,.4,.4,.6,.3,.7,.5,.5/
  ! Data for Hartman FUNCTION coefficients (6-D)
  DATA a2 /10.,0.05,3.,17.,3.,10.,3.5,8.,17.,17.,1.7,.05, 3.5, &
    0.1,10.,10.,1.7,8.,17.,.1,8.,14.,8.,14./
  DATA c2 /1.,1.2,3.,3.2/
  DATA p2 /.1312,.2329,.2348,.4047,.1696,.4135,.1451,.8828, .5569, &
    .8307,.3522,.8732,.0124,.3736,.2883,.5743, &
    .8283,.1004,.3047,.1091,.5886,.9991,.6650,.0381/
  ! Data for Hartman FUNCTION coefficient (3-D)
  DATA a3 /3.,.1,3.,.1,10.,10.,10.,10.,30.,35.,30.,35./
  DATA c3 /1.,1.2,3.,3.2/
  DATA p3 /.3689,.4699,.1091,.03815, .1170,.4387,.8732,.5743, &
    .2673,.7470,.5547,.8828/

  ! IF PARAMETER IS OUTSIDE BOUND, ASSIGN bigm to FUNCTION VALUE
  DO 1 j = 1, nopt
    IF (x(j) .LT. bl(loc(j)) .OR. x(j) .GT. bu(loc(j))) THEN
      functOthers = bigm
      RETURN
    END IF
1 CONTINUE
  ! CHECK nopt AND PASS x(1) AND x(2) TO x1 AND x2 FOR IFunc = 1, 4
  IF (idfunc .GE. 1 .AND. idfunc .LE. 4) THEN
    IF (nopt .NE. 2) THEN
      WRITE(*,*) ' Error in nopt'
      STOP
    END IF
    x1 = x(1)
    x2 = x(2)
  END IF

  ! This is the Goldstein-Price FUNCTION G:3,(0,-1)
  ! (True G:3.0,(0.0,-1.0)
  IF (idfunc .EQ. 1) THEN
    u1 = (x1 + x2 + 1.0)**2
    u2 = 19. - 14.*x1 + 3.*x1**2 - 14.*x2 + 6.*x1*x2 +3.*x2**2
    u3 = (2.*x1 - 3.*x2)**2
    u4 = 18. - 32.*x1 + 12.*x1**2 + 48.*x2 -36.*x1*x2 + 27.*x2**2
    u5 = u1 * u2
    u6 = u3 * u4
    functOthers = (1. + u5) * (30. + u6)
    functOthers = functOthers - 3.0
  END IF
  ! This is the Rosenbrock FUNCTION with Parameter (a) G:0,(1,1)
  IF(idfunc .EQ. 2) THEN
    a = 100.0
    functOthers = a * (x2 - x1**2)**2 + (1 - x1)**2
  END IF
  ! This is the Six-hump Camelback FUNCTION.
  ! (True optima = (-0.08983,0.7126), (0.08983,-0.7126)
  IF (idfunc .EQ. 3) THEN
    functOthers = (4. - 2.1*x1**2 + x1**4/3.)*x1**2 + x1*x2 + &
      (-4. + 4.*x2**2)*x2**2
    functOthers = functOthers + 1.031628453489877
  END IF
  ! This is the Rastrigin FUNCTION G:-2,(0,0)
  IF (idfunc .EQ. 4) THEN
    functOthers = x1**2 + x2**2 - cos(18.0*x1) - cos(18.0*x2)
    functOthers = 2.0 + functOthers
  END IF
  ! This is the Griewank FUNCTION G:0,(0,0), 2-D or 10-D
  IF (idfunc .EQ. 5) THEN
    IF (nopt .EQ. 2) THEN
      d = 200.0
    ELSE
      d = 4000.0
    END IF
    u1 = 0.0
    u2 = 1.0
    DO 10 j = 1, nopt
      u1 = u1 + x(j)**2 / d
      u2 = u2 * cos(x(j)/sqrt(DBLE(j)))
10 CONTINUE
    functOthers = u1 -u2 + 1
  END IF
  ! This is the Shekel FUNCTION G:-10.5364098252,(4,4), 4-D
  IF (idfunc .EQ. 6) THEN
    functOthers = 11.0309996713 - 0.494589854607
    DO 30 i = 1, 10
      u = 0.0
      DO 20 j = 1, nopt
        u = u + (x(j) - a1(i,j))**2
20 CONTINUE
      u = 1.0 / (u + c1(i))
      functOthers = functOthers - u
30 CONTINUE
  END IF
  ! This is the Hartman FUNCTION
  ! G:-3.322368011415515,(0.201,0.150,0.477,0.275,0.311,0.657)
  IF (idfunc .EQ. 7) THEN
    functOthers = 3.322368011415515
    DO 50 i = 1, 4
      u = 0.0
      DO 40 j = 1, nopt
        IF (nopt .EQ. 3) THEN
          a2(i,j) = a3(i,j)
          p2(i,j) = p3(i,j)
        END IF
        u = u + a2(i,j) * (x(j) - p2(i,j))**2
40    CONTINUE
      IF (nopt .EQ. 3) c2(i) = c3(i)
      functOthers = functOthers - c2(i) * EXP(-u)
50  CONTINUE
  END IF

  ! END OF FUNCTION functOthers
  RETURN
  END FUNCTION functOthers

!============================================================
! ccc ------- other FUNCTION conditions
!============================================================
  SUBROUTINE otherf(a,n,nopt,idfunc,loc,bl,bu)
  ! IMPLICIT REAL*8 (a-h,o-z)
  IMPLICIT NONE
  INTEGER n,nopt,idfunc
  DOUBLE PRECISION a,bl,bu
  DIMENSION a(n),bl(50), bu(50)

  INTEGER loc
  DIMENSION loc(50)

  INTEGER i

  DO 111 i=1,50
111 loc(i)=i

  IF(idfunc.EQ.1) THEN
    nopt=2
    a(1)=0
    a(2)=0
    bl(1)=-2
    bu(1)=2
    bl(2)=-2
    bu(2)=2
  ELSE IF(idfunc.EQ.2) THEN
    nopt=2
    a(1)=0
    a(2)=0
    bl(1)=-5
    bu(1)=5
    bl(2)=-2
    bu(2)=8
  ELSE IF(idfunc.EQ.3) THEN
    nopt=2
    a(1)=0
    a(2)=0
    bl(1)=-5
    bu(1)=5
    bl(2)=-5
    bu(2)=5
  ELSE IF(idfunc.EQ.4) THEN
    nopt=2
    a(1)=0
    a(2)=0
    bl(1)=-1
    bu(1)=1
    bl(2)=-1
    bu(2)=1
  ELSE IF(idfunc.EQ.5) THEN
    nopt=10
    DO 10 i=1,nopt
      a(i)=0
      bl(i)=-600
      bu(i)=600
10  CONTINUE
  ELSE IF(idfunc.EQ.6) THEN
    nopt=4
    DO 11 i=1,nopt
      a(i)=5
      bl(i)=0
      bu(i)=10
11  CONTINUE
  ELSE IF(idfunc.EQ.7) THEN
    nopt=6
    DO 12 i=1,nopt
      a(i)=0.5
      bl(i)=0
      bu(i)=1
12  CONTINUE
  END IF

  END SUBROUTINE otherf

!============================================================
!
!============================================================
  SUBROUTINE ex1(par,n)
  ! ******************************************************************
  ! This is the execution SUBROUTINE for the Sacramento
  ! soil moisture accounting operation as modIFied for research
  ! on the University of Arizona campus. No derivative computations
  ! are included in this SUBROUTINE. Last revised 4/25/88 wtw.
  ! ******************************************************************
  !
  IMPLICIT REAL*8 (a-h,o-z)
  REAL*8 lztwc,lzfsc,lzfpc,lztwm,lzfpm,lzfsm,lzpk,lzsk


  COMMON /fsums1/ srot,simpvt,srodt,srost,sintft,sgwfp,sgwfs, &
    srecht,sett,se1,se3,se4,se5

  PARAMETER(ijk=6000,ldt=24)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  ! COMMON /block5/ pt(50),idfunc
  ! COMMON /idid/iduh,area
  COMMON /incntnt/ uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc
  DIMENSION par(n)
  COMMON /pval/ uztwm,uzfwm,uzk,pctim,adimp,riva,zperc,rexp, &
    lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,rserv,side

  ! Evapotranspiration distribution (EPDIST) is unIForm each day.

  epdist = 1.0/ndt

  uztwm=par(1)
  uzfwm=par(2)
  lztwm=par(3)
  lzfpm=par(4)
  lzfsm=par(5)
  adimp=par(6)
  uzk=par(7)
  lzpk=par(8)
  lzsk=par(9)
  pctim=par(10)
  zperc=par(11)
  rexp=par(12)
  pfree=par(13)
  riva=par(14)
  rserv=par(15)
  side=par(16)


  !
  ! Initialize sums and carryover.
  uztwc = state(1)
  uzfwc = state(2)
  lztwc = state(3)
  lzfsc = state(4)
  lzfpc = state(5)
  adimc = state(6)
  IF (uztwc .GT. uztwm) uztwc = uztwm
  IF (uzfwc .GT. uzfwm) uzfwc = uzfwm
  IF (lztwc .GT. lztwm) lztwc = lztwm
  IF (lzfsc .GT. lzfsm) lzfsc = lzfsm
  IF (lzfpc .GT. lzfpm) lzfpc = lzfpm
  ! IF (adimc .GT. adimp) adimc = adimp
  !
  srot = 0.
  simpvt = 0.
  srodt = 0.
  srost = 0.
  sintft = 0.
  sgwfp = 0.
  sgwfs = 0.
  srecht = 0.
  sett = 0.
  sprt = 0.
  spet = 0.
  se1 = 0.
  se3 = 0.
  se4 = 0.
  se5 = 0.
  sro = 0.
  srech = 0.
  set = 0.
  spr = 0.
  !
  ! DT is the length of each time interval in days.
  dt = 1.0/ndt
  !
  ! Begin day and hour loops for ET and PXV calculations.
  DO 105 i = 1, nday
    ep = pet(i)
    ! spet = spet + ep
    DO 100 j = 1, ndt
      pxv = px(i,j)
      ! PerFORM soil moisture accounting operations.
      ! sprt = sprt + pxv

      CALL fland1(dt, pxv, ep, epdist, tlci)

      tci((i - 1) * ndt + j) = tlci
100 CONTINUE
105 CONTINUE
  RETURN
  END SUBROUTINE ex1

!============================================================
!
!============================================================
  SUBROUTINE ex2 !(area)
  ! ***************************************************************
  ! This is a revised subroutine ex2. It converts total channel
  ! inflow to simulated instantaneous streamflow in cms, then
  ! converts simulated instantaneous streamflow to simulated
  ! mean daily streamflow in cms. This subroutine combines
  ! the functions of subroutines ex2 and ex6 of the NWSRFS
  ! model as modified for research on the University of
  ! Arizona campus. Last revised 4/25/88 wtw.
  ! ***************************************************************
  !
  IMPLICIT REAL*8 (a-h,o-z)
  ! DOUBLE PRECISION area
  PARAMETER(ijk=6000,ldt=24)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/ px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  ! COMMON /block5/ pt(50),idfunc
  COMMON /routing/ ckm,cxm,nm
  COMMON /idid/ iduh,area

  !
  ! Convert total channel inflow (TCI, in mm), to simulated
  ! instantaneous streamflow (SQIN, in cms).
  !
  IF(iduh.EQ.1) THEN
    DO 205 i = 1, nday * ndt
      sqin(i) = 0.
      nno = nord
      IF (i .LT. nord) nno = i
      DO 200 j = 1, nno
        sqin(i) = sqin(i) + (tci(i + 1 - j) * uhg(j))
200   CONTINUE
205 CONTINUE
    ! ELSE IF(iduh.EQ.2)THEN
    ! CALL msk(k,x,n,dt)
    !
    ! ELSE IF(iduh.EQ.3) THEN
    ! CALL kinm()
  ELSE
    DO 322 i=1,nday*ndt
322 sqin(i)=tci(i)*area/(3.6*24/ndt)
  END IF

  !
  ! Convert simulated instantaneous streamflow to simulated mean
  ! daily streamflow (SMQ, in cms).
  !
  DO 215 i = 1, nday
    IF (i .EQ. 1) THEN
      smq(i)=0.0
      DO 216 j=1,ndt-1
        smq(i) = smq(i)+sqin(j)
216   CONTINUE
      smq(i)=smq(i)/ndt+sqin(ndt)/(2.0*ndt)
    END IF

    IF (i .NE. 1) THEN
      j = ndt * (i - 1)
      smq(i) = (sqin(j) + sqin(j + ndt)) /(2.0*ndt)
      DO 210 k = 1, ndt-1
        smq(i) = smq(i) + (sqin(j + k) / ndt)
210   CONTINUE
    END IF
215 CONTINUE
  RETURN
  END SUBROUTINE ex2

!============================================================
!
!============================================================
  SUBROUTINE fland1(dt, pxv, ep, epdist, tci)
  ! .......................................
  ! THIS SUBROUTINE EXECUTES THE 'SAC-SMA ' OPERATION FOR ONE TIME
  ! PERIOD.
  ! .......................................
  ! SUBROUTINE INITIALLY WRITTEN BY. . .
  ! ERIC ANDERSON - HRL APRIL 1979 VERSION 1
  ! MODIFIED FEB 1988 FOR RESEARCH ON THE U OF A CAMPUS
  ! Last revised 4/25/88 wtw.
  ! ............................................................
  IMPLICIT REAL*8 (a-h,o-z)

  REAL*8 lztwc,lzfsc,lzfpc,lztwm,lzfpm,lzfsm,lzpk,lzsk

  DIMENSION rsum(7)
  COMMON /incntnt/ uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc
  COMMON /fsums1/ srot,simpvt,srodt,srost,sintft,sgwfp,sgwfs, &
    srecht,sett,se1,se3,se4,se5
  COMMON /pval/ uztwm,uzfwm,uzk,pctim,adimp,riva,zperc,rexp, &
    lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,rserv,side

  ! .......................................
  !
  ! REPARAMETERIZATION OPTION
  ! izflag=2, rePARAMETERized percolation,
  ! otherwise, original percolation
  ! izflag = 1
  zp = zperc
  ! IF (izflag .EQ. 2) aperc = zperc
  !
  ! COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL.
  ! EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL
  !
  ! COMPUTE ET FROM UPPER ZONE.
  edmnd = ep * epdist
  ! RED IS RESIDUAL EVAP DEMAND
  e1 = edmnd * uztwc / uztwm
  red = edmnd - e1
  uztwc = uztwc - e1
  IF (ABS(uztwc) .LT. 0.00001) uztwc = 0.0
  e2 = 0.0
  IF (uztwc .GE. 0.) GOTO 220
  ! E1 CAN NOT EXCEED UZTWC
  e1 = e1 + uztwc
  uztwc = 0.0
  red = edmnd - e1
  IF (uzfwc .GE. red) GOTO 221
  ! E2 IS EVAP FROM UZFWC.
  e2 = uzfwc
  uzfwc = 0.0
  red = red - e2
  GOTO 225
221 e2 = red
  uzfwc = uzfwc - e2
  red = 0.0
220 IF ((uztwc / uztwm) .GE. (uzfwc / uzfwm)) GOTO 225
  ! UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE
  ! TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION
  uzrat = (uztwc + uzfwc) / (uztwm + uzfwm)
  uztwc = uztwm * uzrat
  uzfwc = uzfwm * uzrat
  !
  ! COMPUTE ET FROM THE LOWER ZONE.
  ! COMPUTE ET FROM LZTWC (E3)
225 e3 = red * lztwc / (uztwm + lztwm)
  lztwc = lztwc - e3
  IF (ABS(lztwc) .LT. 0.00001) lztwc = 0.0
  IF (lztwc .GE. 0.0) GOTO 226
  ! E3 CAN NOT EXCEED LZTWC
  e3 = e3 + lztwc
  lztwc = 0.0
226 ratlzt = lztwc / lztwm
  saved = rserv * (lzfpm + lzfsm)

  ratlz = (lztwc+lzfpc+lzfsc-saved) / (lztwm+lzfpm+lzfsm-saved)
  IF (ratlzt .GE. ratlz) GOTO 230
  ! RESUPPLY LOWER ZONE TENSION WATER FROM LOWER
  ! ZONE FREE WATER IF MORE WATER AVAILABLE THERE.
  del = (ratlz - ratlzt) * lztwm
  ! TRANSFER FROM LZFSC TO LZTWC.
  lztwc = lztwc + del
  lzfsc = lzfsc - del
  IF (lzfsc .GE. 0.0) GOTO 230
  ! IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC
  lzfpc = lzfpc + lzfsc
  lzfsc = 0.0
  !
  ! COMPUTE ET FROM ADIMP AREA.-E5
230 e5 = e1 + (red+e2) * (adimc-e1-uztwc) / (uztwm+lztwm)
  ! ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR EVAPORATION.
  adimc = adimc - e5
  IF (ABS(adimc) .LT. 0.00001) adimc = 0.0
  IF (adimc .GE. 0.0) GOTO 231
  ! E5 CAN NOT EXCEED ADIMC.
  e5 = e5 + adimc
  adimc = 0.0
231 e5 = e5 * adimp
  ! E5 IS ET FROM THE AREA ADIMP.
  ! .......................................
  ! COMPUTE PERCOLATION AND RUNOFF AMOUNTS.
  twx = pxv + uztwc - uztwm
  ! TWX IS THE TIME INTERVAL AVAILABLE MOISTURE IN EXCESS
  ! OF UZTW REQUIREMENTS.
  IF (twx .GE. 0.0) GOTO 232
  ! ALL MOISTURE HELD IN UZTW--NO EXCESS.
  uztwc = uztwc + pxv
  twx = 0.0
  GOTO 233
  ! MOISTURE AVAILABLE IN EXCESS OF UZTW STORAGE.
232 uztwc = uztwm
233 adimc = adimc + pxv - twx
  !
  ! COMPUTE IMPERVIOUS AREA RUNOFF.
  roimp = pxv * pctim
  ! ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA.
  simpvt = simpvt + roimp
  !
  ! INITIALIZE TIME INTERVAL SUMS.
  sbf = 0.0
  ssur = 0.0
  sif = 0.0
  sperc = 0.0
  sdro = 0.0
  spbf = 0.0
  !
  ! DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC TIME
  ! INTERVAL
  ! ninc = 1.0 + 0.2 * uzfwc + twx
  ninc = 1.0 + 0.2 * (uzfwc + twx)
  ! NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
  ! IS DIVIDED INTO FOR FURTHER
  ! SOIL-MOISTURE ACCOUNTING. NO ONE INCREMENT
  ! WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV

  dinc = 1.0 / ninc * dt
  ! DINC=LENGTH OF EACH INCREMENT IN DAYS.
  pinc = twx / ninc
  ! PINC=AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT.
  ! COMPUTE FREE WATER DEPLETION FRACTIONS FOR
  ! THE TIME INCREMENT BEING USED-BASIC DEPLETIONS
  ! ARE FOR ONE DAY
  IF (uzk .GT. 1.d0) WRITE(*,*)'uzk = ', uzk
  IF (lzpk .GT. 1.d0) WRITE(*,*)'lzpk = ', lzpk
  IF (lzsk .GT. 1.d0) WRITE(*,*)'lzsk = ', lzsk
  duz = 1.0 - (1.0 - uzk) ** dinc
  dlzp = 1.0 - (1.0 - lzpk) ** dinc
  dlzs = 1.0 - (1.0 - lzsk) ** dinc
  parea = 1.0 - adimp - pctim
  ! .......................................
  ! START INCREMENTAL DO LOOP FOR THE TIME INTERVAL.
  ! .......................................
  DO 240 i = 1, ninc
    adsur = 0.0
    ! COMPUTE DIRECT RUNOFF (FROM ADIMP AREA).
    ratio = (adimc - uztwc) / lztwm
    addro = pinc * ratio ** 2
    sdro = sdro + addro * adimp
    ! ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM
    ! THE AREA ADIMP-SDRO IS THE SIX HOUR SUMMATION
    !
    ! COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
    bf = lzfpc * dlzp
    lzfpc = lzfpc - bf
    IF (lzfpc .GT. 0.0001) GOTO 234
    bf = bf + lzfpc
    lzfpc = 0.0
234 sbf = sbf + bf
    spbf = spbf + bf
    bf = lzfsc * dlzs
    lzfsc = lzfsc - bf
    IF (lzfsc .GT. 0.0001) GOTO 235
    bf = bf + lzfsc
    lzfsc = 0.0
235 sbf = sbf + bf
    !
    ! COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP
    IF ((pinc + uzfwc) .GT. 0.01) GOTO 251
    uzfwc = uzfwc + pinc
    GOTO 249
251 percm = lzfpm * dlzp + lzfsm * dlzs
    ! IF (izflag .EQ. 2) zp = (uzfwm - percm) / (percm * aperc**rexp)
    IF (zp .LT. 0.0) zp = 0.0
    perc = percm * uzfwc / uzfwm

    defr = 1.0 - (lztwc+lzfpc+lzfsc) / (lztwm+lzfpm+lzfsm)

    IF (defr .LT. 0.d0) THEN
      WRITE(*,*)'defr = ', defr
      WRITE(*,*) lztwc, lzfpc, lzfsc
      WRITE(*,*) lztwm, lzfpm, lzfsm
      STOP
    END IF
    ! DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO
    uzdefr = 1.0 - (uztwc + uzfwc) / (uztwm + uzfwm)
239 perc = perc * (1.0 + zp * defr ** rexp)
    ! NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED.
    IF (perc .LT. uzfwc) GOTO 241
    ! PERCOLATION RATE EXCEEDS UZFWC.
    perc = uzfwc
    ! PERCOLATION RATE IS LESS THAN UZFWC.
241 uzfwc = uzfwc - perc
    ! CHECK TO SEE IF PERCOLATION EXCEEDS LOWER ZONE DEFICIENCY.
    check = lztwc+lzfpc+lzfsc+perc-lztwm-lzfpm-lzfsm
    IF (check .LE. 0.0) GOTO 242
    perc = perc - check
    uzfwc = uzfwc + check
242 sperc = sperc + perc
    ! SPERC IS THE TIME INTERVAL SUMMATION OF PERC
    !
    ! COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
    ! NOTE...PINC HAS NOT YET BEEN ADDED
    del = uzfwc * duz
    sif = sif + del
    uzfwc = uzfwc - del
    ! DISTRIBE PERCOLATED WATER INTO THE LOWER ZONES
    ! TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE PFREE AREA.
    ! PERCT IS PERCOLATION TO TENSION WATER AND PERCF IS PERCOLATION
    ! GOING TO FREE WATER.
    perct = perc * (1.0 - pfree)
    IF ((perct + lztwc) .GT. lztwm) GOTO 243
    lztwc = lztwc + perct
    percf = 0.0
    GOTO 244
243 percf = perct + lztwc - lztwm
    lztwc = lztwm
    !
    ! DISTRIBUTE PERCOLATION IN EXCESS OF TENSION
    ! REQUIREMENTS AMONG THE FREE WATER STORAGES.
244 percf = percf + (perc * pfree)
    IF (percf .EQ. 0.0) GOTO 245
    hpl = lzfpm / (lzfpm + lzfsm)
    ! HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE
    ! AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE.
    ratlp = lzfpc / lzfpm
    ratls = lzfsc / lzfsm
    ! RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR
    ! IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE

    fracp = hpl * 2.0 * (1.0-ratlp) / (1.0-ratlp+1.0-ratls)
    ! FRACP IS THE FRACTION GOING TO PRIMARY.
    IF (fracp .GT. 1.0) fracp = 1.0
    percp = percf * fracp
    percs = percf - percp
    ! PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS
    ! PERCOLATION GOING TO PRIMARY AND SUPPLEMENTAL
    ! STORGES,RESPECTIVELY.
    lzfsc = lzfsc + percs
    IF (lzfsc .LE. lzfsm) GOTO 246
    percs = percs - lzfsc + lzfsm
    lzfsc = lzfsm
246 lzfpc = lzfpc + percf - percs
    ! CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM.
    IF (lzfpc .LE. lzfpm) GOTO 245
    excess = lzfpc - lzfpm
    lztwc = lztwc + excess
    lzfpc = lzfpm
    !
    ! DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF.
245 IF (pinc .EQ. 0.0) GOTO 249
    ! CHECK IF PINC EXCEEDS UZFWM
    IF ((pinc + uzfwc) .GT. uzfwm) GOTO 248
    ! NO SURFACE RUNOFF
    uzfwc = uzfwc + pinc
    GOTO 249
    !
    ! COMPUTE SURFACE RUNOFF (SUR) AND KEEP TRACK OF TIME INTERVAL SUM.
248 sur = pinc + uzfwc - uzfwm
    uzfwc = uzfwm
    ssur = ssur + (sur * parea)

    adsur = sur * (1.0 - addro / pinc)
    ! ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES
    ! FROM THAT PORTION OF ADIMP WHICH IS NOT
    ! CURRENTLY GENERATING DIRECT RUNOFF. ADDRO/PINC
    ! IS THE FRACTION OF ADIMP CURRENTLY GENERATING
    ! DIRECT RUNOFF.
    ssur = ssur + adsur * adimp
249 adimc = adimc + pinc - addro - adsur
    IF (adimc .LE. (uztwm + lztwm)) GOTO 240
    addro = (addro + adimc) - (uztwm + lztwm)
    adimc = uztwm + lztwm
240 CONTINUE
  ! .......................................
  ! END OF INCREMENTAL DO LOOP.
  ! .......................................
  ! COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER
  ! WHICH THEY ARE GENERATED.
  eused = e1 + e2 + e3
  ! EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM
  sif = sif * parea
  !
  ! SEPARATE CHANNEL COMPONENT OF BASEFLOW
  ! FROM THE NON-CHANNEL COMPONENT
  tbf = sbf * parea
  ! TBF IS TOTAL BASEFLOW

  bfcc = tbf * (1.0 / (1.0 + side))
  ! BFCC IS BASEFLOW, CHANNEL COMPONENT
  bfp = (spbf * parea) / (1.0 + side)
  bfs = bfcc - bfp
  IF (bfs .LT. 0.0) bfs = 0.0
  bfncc = tbf - bfcc
  ! BFNCC IS BASEFLOW,NON-CHANNEL COMPONENT
  !
  ! ADD TO MONTHLY SUMS.
  sintft = sintft + sif
  sgwfp = sgwfp + bfp
  sgwfs = sgwfs + bfs
  srecht = srecht + bfncc
  srost = srost + ssur
  srodt = srodt + sdro
  !
  ! COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL.
  tci = roimp + sdro + ssur + sif + bfcc
  !
  ! COMPUTE E4-ET FROM RIPARIAN VEGETATION.
  e4 = (edmnd - eused) * riva
  !
  ! SUBTRACT E4 FROM CHANNEL INFLOW
  tci = tci - e4
  IF (tci .GE. 0.0) GOTO 250
  e4 = e4 + tci
  tci = 0.0
250 srot = srot + tci
  !
  ! COMPUTE TOTAL EVAPOTRANSPIRATION-TET
  eused = eused * parea
  tet = eused + e5 + e4
  sett = sett + tet
  se1 = se1 + e1 * parea
  se3 = se3 + e3 * parea
  se4 = se4 + e4
  se5 = se5 + e5
  ! CHECK THAT ADIMC IS GT UZTWC
  !
  IF (adimc .LT. uztwc) adimc = uztwc
  !
  ! ADD TO SUMS OF RUNOFF COMPONENTS.

  ! .......................................
  RETURN
  END SUBROUTINE fland1

!============================================================
!
!============================================================
  REAL*8 FUNCTION functnsc(nopt,x)
  !
  !
  ! THIS IS THE SIXPAR MODEL
  ! Bound: PAR(1)=UM=[0,50], PAR(2)=UK=[0,1], PAR(3)=BM=[0,50],
  ! PAR(4)=BK=[0,1], PAR(5)=A=[0,1], PAR(6)=X=[0,10]
  ! Global Optimum: 0, [10, 0.5, 20, 0.2, 0.31, 3]
  !
  IMPLICIT REAL*8 (a-h,o-z)
  PARAMETER(ijk=6000,ldt=24,bigm = 1.0d10)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  COMMON /iopar/ ipr
  COMMON /ablu/ bl(50), bu(50)
  DIMENSION x(nopt), par(50)

  ! IF PARAMETER IS OUTSIDE BOUND, ASSIGN bigm to FUNCTION VALUE
  DO 1 j = 1, nopt
    IF (x(j) .LT. bl(loc(j)) .OR. x(j) .GT. bu(loc(j))) THEN
      functnsc = bigm
      RETURN
    END IF
1 CONTINUE
  !
  IF (nopt .GT. 16) THEN
    WRITE (*,*) ' ???? CHECK THE INPUT FILE ???? '
    STOP
  END IF
  !
  ! ASSIGN x TO APPROPRIATE PARAMETERS IN SIXPAR MODEL
  DO 10 j = 1, npar
    par(j) = pt(j)
10 CONTINUE
  DO 20 k = 1, nopt
    par(loc(k)) = x(k)
20 CONTINUE


  !
  ! COMPUTE THE STREAMFLOW FOR PARAMETER SET par
  CALL model(par,npar)
  !
  ! COMPUTE THE SIMPLE LEAST SQUARE OBJECTIVE FUNCTION VALUE
  IF (iobj .EQ. 1) THEN
    functnsc = 0.d0
    DO 40 i = 1, ndata
      functnsc = functnsc + (qobs(i) - smq(i))**2
40  CONTINUE
  END IF
  !
  ! COMPUTE THE HMLE VALUE
  IF (iobj .EQ. 2) functnsc = hmle()
  !
  ! END OF FUNCTION functnsc
  RETURN
  END FUNCTION functnsc


!============================================================
!
!============================================================
  REAL*8 FUNCTION hmle()
  IMPLICIT REAL*8 (a-h,o-z)

  COMMON /fnblk/ rlamda, ad
  PARAMETER(ijk=6000,ldt=24)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  COMMON /iopar/ ipr

  DIMENSION ra(2)
  PARAMETER (eps=5.d-02, del=5.d-02)
  DATA iflag /0/
  !
  ! COMPUTE THE MEAN OF LOGARITHM OF OBSERVED FLOWS
  IF (iflag .EQ. 0) THEN
    ad = 0.d0
    DO 10 i = 1, ndata
      ad = ad + DLOG(qobs(i))
10  CONTINUE
    ad = ad / DBLE(ndata)
    rlamda = 1.d0
    iflag = 1
  END IF
  !
  ! ESTIMATE THE LAMDA VALUE
  lcount = 0
  ict = 1
  ra(1) = 0.d0
  ra(2) = 0.d0
25 CONTINUE
  lcount = lcount + 1
  IF(lcount .GT. 40) THEN
    WRITE(*,*) 'LAMDA ITERATION GO OVER 40', rlamda, ra(1), ra(2)
    GO TO 50
  END IF
  rd = 0.d0
  rn = 0.d0
  DO 30 i = 1, ndata
    a = DLOG(qobs(i)) / ad
    w = qobs(i)**(2*(rlamda-1.d0))
    rd = rd + w*(qobs(i) - smq(i))**2
    rn = rn + w*(qobs(i) - smq(i))**2 * a
30 CONTINUE
  ra(ict) = rn / rd - 1.d0
  IF (DABS(ra(ict)) .LE. eps) GO TO 50
  isign = -1
  IF (ra(ict) .LT. 0.d0) isign = 1
  rlamda = rlamda + isign * del
  IF (ict .EQ. 2) GO TO 35
  ict = 2
  GO TO 25
  !
35 CONTINUE
  IF (ra(1)*ra(2) .LT. 0.d0) GO TO 40
  ra(1) = ra(2)
  GO TO 25
  !
40 CONTINUE
  rlamda = rlamda - isign * del / 2.d0
  !
  ! COMPUTE HMLE
50 CONTINUE
  hmle = 0.d0
  ex = 2. * (rlamda - 1.)
  DO 60 i = 1, ndata
    hmle = hmle + qobs(i)**ex * (qobs(i) -smq(i))**2
60 CONTINUE
  hmle = hmle / DBLE(ndata)
  hmle = hmle / DEXP(ex * ad)
  !
  RETURN
  END FUNCTION hmle


!============================================================
!
!============================================================
  SUBROUTINE model(par,n)
  ! **********************************************************************
  ! This is the execution controlLEr for the SAC-SMA model. It controls
  ! operation of the streamflow simulation model for each set of
  ! PARAMETER values. Parameter values are set in the NWS main program
  ! IF no PARAMETER estimation is perFORMed. SeLEcted PARAMETER values
  ! are varied in the OARSCH SUBROUTINE IF PARAMETER estimation is
  ! perFORMed. Last revised 2/1/89 wtw.
  ! **********************************************************************
  !
  IMPLICIT REAL*8 (a-h,o-z)

  PARAMETER(ijk=6000,ldt=24)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
!-ygan
!  COMMON /idid/iduh,area
!.ygan
  DIMENSION par(n)

  ! Calculate total channel inflow (TCI) in mm.

  CALL ex1(par,n)

  CONTINUE
  !
  ! Convert channel inflow in mm to simulated instantaneous streamflow
  ! (SQIN) and simulated mean daily streamflow (SMQ) in cms
  ! CALL ex2(tci, uhg, sqin, smq)

  CALL ex2 !(area)
  CONTINUE
  RETURN
  END SUBROUTINE model

!============================================================
!
!============================================================
  SUBROUTINE sxpinp(nopt,a,sPath)
  !
  !
  ! THIS SUBROUTINE READS THE INPUT DATA FOR SIXPAR MODEL
  !
  ! SIXPAR MODEL VARIABLES AND PARAMETERS
  !      ndata = number of data points
  !      idata = flag on what data set to use for calibration
  !            = 0, use existing data
  !            = 1, generate data within the program
  !      iobj = flag on which objective function to use
  !            = 1, simple least square (SLS)
  !            = 2, heteroscedastic maximum likelihood estimator (HMLE)
  !      pt(i) =.true.value for ith parameter, i=1,2,...,6
  !      loc(i) = index on the order of PARAMETERs to be optimized
  !      ns = number of state variabLEs
  !      p(.) = precipitation data
  !      obsq(.) = observed streamflow data
  !      s0(.) = initial states
  !
  IMPLICIT REAL*8 (a-h,o-z)

  PARAMETER(ijk=6000,ldt=24)
  COMMON /block1/ ndata,ns,iobj,nday,nord,ndt
  COMMON /block2/px(ijk,ldt),pet(ijk),qobs(ijk),uhg(50),state(20)
  COMMON /block3/ smq(ijk),sqin(ijk*ldt),tci(ijk*ldt)
  COMMON /block4/ npar,loc(50)
  COMMON /block5/ pt(50),idfunc
  COMMON /idid/iduh,area
  COMMON /iopar/ ipr
  INTEGER in
  CHARACTER*4 name(50) ! haracter*3 name(50)
  CHARACTER*4 objty,obj1,obj2
  CHARACTER*9 dataty,data1,data2
  DATA obj1/'SLS '/
  DATA obj2/'HMLE'/
  DATA data1/'OBSERVED '/
  DATA data2/'SYNTHETIC'/
  DATA name /'  X1','  X2','  X3','  X4','  X5','  X6','  X7',&
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X16',&
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X25', &
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X34', &
    '  X8','  X9',' X10',' X11',' X12',' X13',' X14',' X15',' X43', &
    '  X8','  X9',' X10',' X11',' X12',' X13',' X50'/
  DIMENSION a(50)
  CHARACTER*(*)sPath

  WRITE (*,*) ' ENTER SUBROUTINE SXPINP --- '

  WRITE(ipr,800)
800 FORMAT(//,10X,'READ AND WRITE THE INPUT INFORMATION FOR SIXPAR', &
    ' MODEL',/,10X,53(1H=))

  in=7

  OPEN(in, FILE=sPath // 'sxpinp.dat')

  READ(in,*) npar,ns,iobj,idata,ndt

  DO 111 i=1,nopt
111 pt(i)=a(i)

  READ(in,*) (pt(j),j=nopt+1,npar)

  READ(in,*)iduh,area

  READ(in,*) (state(j),j=1,ns)

  ! routing by UH
  IF (iduh.EQ.1) THEN
    READ(in,*)nord
    READ(in,*)(uhg(i),i=1,nord)
    ! By MushK---
  ELSEIF (iduh.EQ.2) THEN
    READ(in,*)ckm,cxm,nm
    ! by Kinematic method----
  ELSEIF(iduh.EQ.3) THEN
    READ(in,*)c
  END IF

  ! READ THE PRE-SELECTED DATA SET

  i=0
777 i=i+1
  READ(in,*,END=778) iid,iim,iiy,qobs(i),pet(i),(px(i,j),j=1,ndt)

  ! qobs(i)=qobs(i)*24*3.6/area iid,

  GO TO 777
778 ndata=i-1


  ! 50 FORMAT(i5,2f10.4)
100 CONTINUE

  nday=ndata

  ! IF idata IS EQUAL TO 1, GENERATE ERROR-FREE SYNTHETIC DATA
  ! BY CALLING THE SIXPAR MODEL USING PARAMETER SET pt
  IF (idata .EQ. 1) THEN
    ft = functn(nopt,pt,idfunc)
    DO 200 i = 1, ndata
      qobs(i) = smq(i)
200 CONTINUE
    dataty = data2
  ELSE
    dataty = data1
  END IF

  IF (iobj .EQ. 1) objty = obj1
  IF (iobj .EQ. 2) objty = obj2

  ! PRINT THE INPUT INFORMATION FOR THE SIXPAR MODEL
  WRITE(ipr,810) npar, ndata, objty, dataty
810 FORMAT(//,2X,'PARAMETERS',4X,'DATA POINTS',4X,'CRITERION',4X,  &
    'DATA TYPE',/,2X,10(1H-),4X,11(1H-),4X,9(1H-),4X,9(1H-), &
    /,6X,i2,10X,i4,12X,a4,6X,a9)
  WRITE(ipr,820) (pt(j), j = 1, npar)
820 FORMAT(//,2X,'TRUE PARAMETER VALUES:',//,5X, &
    'UM',8X,'UK',8X,'BM',8X,'BK',8X,'A ',8X,'X ',/,6(F8.2,2X))
  WRITE(ipr,830) (name(loc(j)),j = 1, nopt)
830 FORMAT(//,2X,'THE PARAMETERS TO BE OPTIMIZED:',//,1X,6(4X,A2,4X))
  !
  ! CLOSE THE SIXPAR MODEL INPUT FILE
  IF(idfunc/=55)THEN
    CLOSE(in)
  END IF
  !
  ! END OF SUBROUTINE SXPINP
  RETURN
  END SUBROUTINE sxpinp
