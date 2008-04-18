MODULE deck_control_block

  USE shared_data
  USE strings

SAVE 
  INTEGER,PARAMETER :: ControlBlockElements = 26
  LOGICAL, DIMENSION(ControlBlockElements) :: ControlBlockDone =.FALSE.
  CHARACTER(len=30),DIMENSION(ControlBlockElements) :: ControlBlockName = (/"nx","ny","nprocx","nprocy","x_stretch","y_stretch",&
       "nsteps","t_end","dt_snapshots","x_start","x_end","y_start","y_end",&
         "visc1","visc2","visc3","resistivemhd","j_max","eta0","eta_background","dt_multiplier","gamma",&
         "rke","data_dir","restart","restart_snapshot"/)

CONTAINS

  FUNCTION HandleControlDeck(Element,Value)

    CHARACTER(*),INTENT(IN) :: Element,Value
    INTEGER :: HandleControlDeck
    INTEGER :: loop,elementselected
    HandleControlDeck=ERR_UNKNOWN_ELEMENT

    elementselected=0

    DO loop=1,ControlBlockElements
       IF(StrCmp(Element,TRIM(ADJUSTL(ControlBlockName(loop))))) THEN
          elementselected=loop
          EXIT
       ENDIF
    ENDDO

    IF (elementselected .EQ. 0) RETURN
    IF (ControlBlockDone(elementselected)) THEN
       HandleControlDeck=ERR_PRESET_ELEMENT
       RETURN
    ENDIF
    ControlBlockDone(elementselected)=.TRUE.
    HandleControlDeck=ERR_NONE

    SELECT CASE (elementselected)
    CASE(1)
       !nx
       nx_global=AsInteger(Value,HandleControlDeck)
    CASE(2)
       !ny
       ny_global=AsInteger(Value,HandleControlDeck)
    CASE(3)
       !nprocx
       nprocx=AsInteger(Value,HandleControlDeck)
    CASE(4)
       !nprocy
       nprocy=AsInteger(Value,HandleControlDeck)
    CASE(5)
       x_stretch=AsLogical(Value,HandleControlDeck)
    CASE(6)
       y_stretch=AsLogical(Value,HandleControlDeck)
    CASE(7)
       nsteps=AsInteger(Value,HandleControlDeck)
    CASE(8)
       t_end=AsReal(Value,HandleControlDeck)
    CASE(9)
       dt_snapshots=AsReal(Value,HandleControlDeck)
    CASE(10)
       x_start=AsReal(Value,HandleControlDeck)
    CASE(11)
       x_end=AsReal(Value,HandleControlDeck)
    CASE(12)
       y_start=AsReal(Value,HandleControlDeck)
    CASE(13)
       y_end=AsReal(Value,HandleControlDeck)
    CASE(14)
       visc1=AsReal(Value,HandleControlDeck)
    CASE(15)
       visc2=AsReal(Value,HandleControlDeck)
    CASE(16)
       visc3=AsReal(Value,HandleControlDeck)
    CASE(17)
       ResistiveMHD=AsLogical(Value,HandleControlDeck)
    CASE(18)
       j_max=AsReal(Value,HandleControlDeck)
    CASE(19)
       eta0=AsReal(Value,HandleControlDeck)
    CASE(20)
       eta_background=AsReal(Value,HandleControlDeck)
    CASE(21)
       dt_multiplier=AsReal(Value,HandleControlDeck)
    CASE(22)
       gamma=AsReal(Value,HandleControlDeck)
    CASE(23)
       rke=AsLogical(Value,HandleControlDeck)
    CASE(24)
       data_dir=Value(1:MIN(LEN(Value),Data_Dir_Max_Length))
    CASE(25)
       restart=AsLogical(Value,HandleControlDeck)
    CASE(26)
       restart_snapshot=AsInteger(Value,HandleControlDeck)
    END SELECT

  END FUNCTION HandleControlDeck

END MODULE
