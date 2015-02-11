PROGRAM visualize_HEX8_chunk_w_medit

  implicit none

  integer            :: np, nhex, nquad, ios, i, nsph, nq1, nq2, nq3, nq4, nq5, nq6
  integer            :: s1, s2, s3, s4, s5, s6, s7, s8, sf1, sf2, sf3, sf4
  real*8             :: p1, p2, p3

  open(unit=10,file='nodes_coords_file',status='unknown',iostat=ios)
  open(unit=11,file='mesh_file',status='unknown',iostat=ios)
  open(unit=31,file='absorbing_surface_file_bottom',status='unknown',iostat=ios)
  open(unit=32,file='absorbing_surface_file_xmax',status='unknown',iostat=ios)
  open(unit=33,file='absorbing_surface_file_xmin',status='unknown',iostat=ios)
  open(unit=34,file='absorbing_surface_file_ymax',status='unknown',iostat=ios)
  open(unit=35,file='absorbing_surface_file_ymin',status='unknown',iostat=ios)
  open(unit=36,file='free_surface',status='unknown',iostat=ios)

  open(unit=20,file='test_out_HEX8.mesh',status='unknown',iostat=ios)

  READ(10,*) np
  READ(11,*) nhex

  READ(31,*) nq1
  READ(32,*) nq2
  READ(33,*) nq3
  READ(34,*) nq4
  READ(35,*) nq5
  READ(36,*) nq6

  nquad = nq1 + nq2 + nq3 + nq4 + nq5 + nq6

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
  WRITE(20,*) nquad


  DO i=1,nq1
    READ(31,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 2
  enddo
  DO i=1,nq2
    READ(32,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 2
  enddo
  DO i=1,nq3
    READ(33,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 2
  enddo
  DO i=1,nq4
    READ(34,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 2
  enddo
  DO i=1,nq5
    READ(35,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 2
  enddo
  DO i=1,nq6
    READ(36,*) nsph, sf1, sf2, sf3, sf4
    WRITE(20,*) sf1, sf2, sf3, sf4, 3
  enddo

  WRITE(20,*) ' '
  WRITE(20,*) 'Hexahedra'
  WRITE(20,*) nhex

  DO i=1,nhex
    READ(11,*) nsph, s1, s2, s3, s4, s5, s6, s7, s8
    WRITE(20,*) s1, s2, s3, s4, s5, s6, s7, s8, 1
  enddo

  WRITE(20,*) ' '
  WRITE(20,*) 'End'

  close(10)
  close(11)
  close(12)
  close(20)

END PROGRAM visualize_HEX8_chunk_w_medit
