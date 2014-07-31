      PROGRAM FINDIT_NEW


      INTEGER i
      CHARACTER*10 ID
      CHARACTER*9 TIM
c     REAL*8 TIM
      CHARACTER*3 INFILE
      CHARACTER*11 PSLW,SSLW,PPSLW,DUM
c     REAL*8 PSLW,SSLW,PPSLW,DUM
      DIMENSION ID(8)
      DIMENSION TIM(8)
      DIMENSION DUM(8)
     
      PSLW=' -12345'
      SSLW=' -12345'
      PPSLW=' -12345'

      INFILE='dat'
      OPEN(UNIT=7,FILE=INFILE,STATUS='OLD')

c initialize all variables to -12345

      DO 33 I=1,8
	id(i)=' -12345'
	tim(i)=' -12345'
	dum(i)=' -12345'
33    CONTINUE
      
      DO 10 I=1,8

      READ(7,1000,END=99) ID(I),TIM(I),DUM(I)
1000   FORMAT(8x,a10,2x,a9,2x,a11)
      IF (ID(I).EQ.'  P       ') THEN
      PSLW=DUM(I)
      ELSE IF (ID(I).EQ.'  S       ') THEN
      SSLW=DUM(I)
      ELSE IF (ID(I).EQ.'  PP      ') THEN
      PPSLW=DUM(I)
      ELSE IF (ID(I).EQ."  P'P'    ") THEN
      ID(I)='  PPdash  '
      ELSE IF (ID(I).EQ."  P'P'df  ") THEN
      ID(I)='  PPdashdf' 
      ELSE IF (ID(I).EQ."  P'P'ab  ") THEN
      ID(I)='  PPdashab'
      ELSE IF (ID(I).EQ."  P'P'bc  ") THEN
      ID(I)='  PPdashbc'
      ENDIF

10    CONTINUE
      
99    CLOSE(7)

      OPEN(UNIT=2,FILE='dat1',STATUS='NEW')
      write(2,*) id(1)
      write(2,*) id(2)
      write(2,*) id(3)
      write(2,*) id(4)
      write(2,*) id(5)
      write(2,*) id(6)
      write(2,*) id(7)
      write(2,*) id(8)
      write(2,*) tim(1)
      write(2,*) tim(2)
      write(2,*) tim(3)
      write(2,*) tim(4)
      write(2,*) tim(5)
      write(2,*) tim(6)
      write(2,*) tim(7)
      write(2,*) tim(8)
      write(2,*) pslw
      write(2,*) sslw
      write(2,*) ppslw

      close(2)

      END
