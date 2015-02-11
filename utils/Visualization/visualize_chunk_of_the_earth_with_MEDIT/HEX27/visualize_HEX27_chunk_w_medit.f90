PROGRAM visualize_HEX27_chunk_w_medit

  implicit none

  integer   :: np, nhex, nquad, ios, i, nsph, nq1, nq2, nq3, nq4, nq5, nq6, nhexdecoup, nqdecoup
  integer   :: s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, &
               s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27
  integer   :: sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
  real*8    :: p1, p2, p3

  open(unit=10,file='nodes_coords_file',status='unknown',iostat=ios)
  open(unit=11,file='mesh_file',status='unknown',iostat=ios)
  open(unit=31,file='absorbing_surface_file_bottom',status='unknown',iostat=ios)
  open(unit=32,file='absorbing_surface_file_xmax',status='unknown',iostat=ios)
  open(unit=33,file='absorbing_surface_file_xmin',status='unknown',iostat=ios)
  open(unit=34,file='absorbing_surface_file_ymax',status='unknown',iostat=ios)
  open(unit=35,file='absorbing_surface_file_ymin',status='unknown',iostat=ios)
  open(unit=36,file='free_surface',status='unknown',iostat=ios)

  open(unit=20,file='test_out_HEX27.mesh',status='unknown',iostat=ios)

  READ(10,*) np
  READ(11,*) nhex
  READ(31,*) nq1
  READ(32,*) nq2
  READ(33,*) nq3
  READ(34,*) nq4
  READ(35,*) nq5
  READ(36,*) nq6

  nhexdecoup  = 8*nhex
  nqdecoup    = 4*(nq1 + nq2 + nq3 + nq4 + nq5 + nq6)

  WRITE(20,*) 'MeshVersionFormatted 1'
  WRITE(20,*) 'Dimension 3'
  WRITE(20,*) ' '
  WRITE(20,*) 'Vertices'
  WRITE(20,*) np

  DO i=1,np
    READ(10,*) nsph, p1, p2, p3
    WRITE(20,*) p1, p2, p3, 1
  enddo

  WRITE(20,*) ' '
  WRITE(20,*) 'Quadrilaterals'
  WRITE(20,*) nqdecoup


  DO i=1,nq1
    READ(31,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 2
    WRITE(20,*) sf5, sf2, sf6, sf9, 2
    WRITE(20,*) sf9, sf6, sf3, sf7, 2
    WRITE(20,*) sf8, sf9, sf7, sf4, 2
  enddo
  DO i=1,nq2
    READ(32,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 2
    WRITE(20,*) sf5, sf2, sf6, sf9, 2
    WRITE(20,*) sf9, sf6, sf3, sf7, 2
    WRITE(20,*) sf8, sf9, sf7, sf4, 2
  enddo
  DO i=1,nq3
    READ(33,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 2
    WRITE(20,*) sf5, sf2, sf6, sf9, 2
    WRITE(20,*) sf9, sf6, sf3, sf7, 2
    WRITE(20,*) sf8, sf9, sf7, sf4, 2
  enddo
  DO i=1,nq4
    READ(34,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 2
    WRITE(20,*) sf5, sf2, sf6, sf9, 2
    WRITE(20,*) sf9, sf6, sf3, sf7, 2
    WRITE(20,*) sf8, sf9, sf7, sf4, 2
  enddo
  DO i=1,nq5
    READ(35,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 2
    WRITE(20,*) sf5, sf2, sf6, sf9, 2
    WRITE(20,*) sf9, sf6, sf3, sf7, 2
    WRITE(20,*) sf8, sf9, sf7, sf4, 2
  enddo
  DO i=1,nq6
    READ(36,*) nsph, sf1, sf2, sf3, sf4, sf5, sf6, sf7, sf8, sf9
    WRITE(20,*) sf1, sf5, sf9, sf8, 3
    WRITE(20,*) sf5, sf2, sf6, sf9, 3
    WRITE(20,*) sf9, sf6, sf3, sf7, 3
    WRITE(20,*) sf8, sf9, sf7, sf4, 3
  enddo

  WRITE(20,*) ' '
  WRITE(20,*) 'Hexahedra'
  WRITE(20,*) nhexdecoup

  DO i=1,nhex
    READ(11,*) nsph, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, &
               s15, s16, s17, s18, s19, s20, s21, s22, s23, s24, s25, s26, s27
!
    WRITE(20,*) s1 , s9 , s21, s12, s13, s22, s27, s25, 1
    WRITE(20,*) s9 , s2 , s10, s21, s22, s14, s23, s27, 1
    WRITE(20,*) s21, s10, s3 , s11, s27, s23, s15, s24, 1
    WRITE(20,*) s12, s21, s11, s4 , s25, s27, s24, s16, 1
    WRITE(20,*) s13, s22, s27, s25, s5 , s17, s26, s20, 1
    WRITE(20,*) s22, s14, s23, s27, s17, s6 , s18, s26, 1
    WRITE(20,*) s27, s23, s15, s24, s26, s18, s7 , s19, 1
    WRITE(20,*) s25, s27, s24, s16, s20, s26, s19, s8 , 1
  enddo

  WRITE(20,*) ' '
  WRITE(20,*) 'End'

  close(10)
  close(11)
  close(12)
  close(20)

END PROGRAM visualize_HEX27_chunk_w_medit
