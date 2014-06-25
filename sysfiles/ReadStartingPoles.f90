! If <useActiveDamping>, reads <gamma> values from FN_GAMMA
! If file could not be opened, Message is called and gamma is set to zero
! If read error occurs, Message is called and program is terminated

SUBROUTINE ReadStartingPoles
	USE Declare
	USE FilesModule
	USE MessageModule
	IMPLICIT NONE
	INTEGER ios
	CHARACTER(12) IntToString

	IF (useExtStartingPoles) THEN
		ALLOCATE (StartingPoles(n), STAT = err)
		IF (err /= 0) CALL AllocationError

		StartingPoles = 0.06d0
		OPEN (FH_POLES, FILE = FN_POLES, IOSTAT = ios)
	
		IF (ios /=0) THEN
!			CALL Message (MessageText(POLESFILEOPENERROR), MESSAGETITLE(POLESFILEOPENERROR)) 
			WRITE (*,*) 'StartingPoles-file open error'
		ELSE
			DO i=1, n
				READ (FH_POLES, *, IOSTAT = ios) StartingPoles(i)
				IF (ios /=0) THEN
!					messageboxQQresult = MESSAGEBOXQQ ("Only " // TRIM(IntToString(i-1)) // " startingpole-values could be read, while "//TRIM(IntToString(n)) // " values\nare needed (one for every section). The startingpoles for the last sections\nwill be set to 0.06"C, "Error reading "//FN_POLES // ""C, MB$OK .OR. MB$ICONINFORMATION) 
					WRITE (*,*) 'Error reading StartingPoles-values !'
					EXIT
				ENDIF
			ENDDO
			CLOSE (FH_POLES)
		ENDIF

!Read in the Vthresholds that go with the startingpoles to have a fixed nonlinearity threshold across the model
		ALLOCATE (Vthresholds(n), STAT = err)
		IF (err /= 0) CALL AllocationError

		Vthresholds = 4.3652d-6
		OPEN (FH_VTHS, FILE = FN_VTHS, IOSTAT = ios)
	
		IF (ios /=0) THEN
!			CALL Message (MessageText(POLESFILEOPENERROR), MESSAGETITLE(POLESFILEOPENERROR)) 
			WRITE (*,*) 'StartingPoles-file open error'
		ELSE
			DO i=1, n
				READ (FH_VTHS, *, IOSTAT = ios) Vthresholds(i)
				IF (ios /=0) THEN
!					messageboxQQresult = MESSAGEBOXQQ ("Only " // TRIM(IntToString(i-1)) // " Vthresholds-values could be read, while "//TRIM(IntToString(n)) // " values\nare needed (one for every section). The Vthresholds for the last sections\nwill be set to default"C, "Error reading "//FN_POLES // ""C, MB$OK .OR. MB$ICONINFORMATION) 
					WRITE (*,*) 'Error reading StartingPoles-values !'
					EXIT
				ENDIF
			ENDDO
			CLOSE (FH_VTHS)
		ENDIF
	ENDIF

END SUBROUTINE ReadStartingPoles
