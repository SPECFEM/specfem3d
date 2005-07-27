
!! Dimitri Komatitsch, Univ. of Pau, April 2004

!! draw polarograms for 3D seismograms in copper crystal for study with Jose Carcione

  program plot_polarograms

!! DK DK UGLY regarding rotation for polarization below:
!! DK DK UGLY check if this is really what we want to do in an anisotropic medium

  implicit none

  include "constants.h"

! normalize displacement globally or for each trace
  logical, parameter :: NORMALIZE_GLOBALLY = .false.

! half-size of window around the source in which we mute displacement
! window extends from - NEX_MUTE to + NEX_MUTE around the source
  integer, parameter :: NEX_MUTE = 5

! generate *.ps or *.eps
  logical, parameter :: GENERATE_EPS = .true.

! total length of the time axis in centimeters
  double precision, parameter :: LENGTH_TIME_AXIS = 24.5d0
  double precision, parameter :: additional_shift_arrow = 0.5d0

! max length to which all the vectors are normalized globally in cm
  double precision, parameter :: SIZEVECT = 5.8d0

! PostScript parameters
  double precision, parameter :: orig_x = 2.2d0
  double precision, parameter :: orig_z = 9.4d0
  double precision, parameter :: length_arrow_cm = 0.35d0

! size of A4 page for PostScript display
  double precision, parameter :: sizex = 21.d0
  double precision, parameter :: sizez = 29.7d0

  character(len=150) filename

! to store seismograms everywhere on vertical face for polarograms
! use single precision to store this array
  real(kind=4), allocatable, dimension(:,:,:) :: displ2D_all_time_steps
  logical, allocatable, dimension(:) :: elem2D_has_been_muted

  integer it,ispec2D,ifigure_latex
  integer ispec2D_xi,ispec2D_gamma,ispec2D_xi_source,ispec2D_gamma_source
  integer NSPEC2D

  double precision xval,zval,xpos,zpos,xend,zend,tscale,maxnorm
  double precision x_diff_source,z_diff_source,dist_diff_source,angle_diff_source_rad,angle_diff_source_deg

! parameters read from parameter file
  integer NEX_GAMMA,NEX_ETA,NEX_XI,NSEIS,NSTEP,SOURCE_TIME_FUNCTION

  double precision UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK
  double precision DT,LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX,FACTOR_SOURCE,SOURCE_DOMINANT_FREQ

  logical POSTSCRIPT_SNAPSHOTS

  character(len=150) LOCAL_PATH

!--------- program starts here

  print *,'**** reading seismograms to construct polarograms ****'
  print *

! read the parameter file
  call read_parameter_file(LAT_MIN,LAT_MAX,LONG_MIN,LONG_MAX, &
        UTM_X_MIN,UTM_X_MAX,UTM_Y_MIN,UTM_Y_MAX,Z_DEPTH_BLOCK, &
        NEX_GAMMA,NEX_ETA,NEX_XI,NSEIS,NSTEP,DT, &
        LOCAL_PATH,FACTOR_SOURCE,POSTSCRIPT_SNAPSHOTS,SOURCE_TIME_FUNCTION,SOURCE_DOMINANT_FREQ)

! number of spectral elements on vertical face Ymin
  NSPEC2D = NEX_XI*NEX_GAMMA

  print *,'total number of time samples for seismograms = ',NSTEP
  print *,'subsampling rate for polarograms = ',ISAMP_POLAROGRAMS
  print *,'total number of time samples for polarograms = ',NSTEP / ISAMP_POLAROGRAMS
  print *
  print *,'number of spectral elements on vertical face Ymin = ',NSPEC2D
  print *,'total number of polarograms = ',NSPEC2D
  print *

! check that even number of elements were used
  if(mod(NEX_XI,2) /= 0) stop 'NEX_XI must be even'
  if(mod(NEX_GAMMA,2) /= 0) stop 'NEX_GAMMA must be even'

  allocate(displ2D_all_time_steps(3,NSPEC2D,NSTEP))
  allocate(elem2D_has_been_muted(NSPEC2D))

! read seismograms everywhere on vertical face for polarograms
! use unformatted file format for faster disk access
  open(unit=IOUT,file='seismos/displ2D_all_time_steps.bin',status='old',form='unformatted')
  read(IOUT) displ2D_all_time_steps
  close(IOUT)

!! DK DK UGLY use analytical function for tests
! do ispec2D=1,NSPEC2D
!   do it=1,NSTEP
!     displ2D_all_time_steps(1,ispec2D,it) = 50.*cos(2.*3.141592*it/1000.)
!     displ2D_all_time_steps(2,ispec2D,it) = 35.*sin(2.*3.141592*it/800.)
!     displ2D_all_time_steps(3,ispec2D,it) = 20.*sin(2.*3.141592*it/800.)
!   enddo
! enddo

! print information about displacement values read
  maxnorm = maxval(sqrt(displ2D_all_time_steps(1,:,:)**2 + displ2D_all_time_steps(2,:,:)**2 + displ2D_all_time_steps(3,:,:)**2))

  print *
  print *,'Amplitude before muting near source and edges:'
  print *
  print *,'minval maxval displ_1 read = ',minval(displ2D_all_time_steps(1,:,:)),maxval(displ2D_all_time_steps(1,:,:))
  print *,'minval maxval displ_2 read = ',minval(displ2D_all_time_steps(2,:,:)),maxval(displ2D_all_time_steps(2,:,:))
  print *,'minval maxval displ_3 read = ',minval(displ2D_all_time_steps(3,:,:)),maxval(displ2D_all_time_steps(3,:,:))
  print *
  print *,'maxval norm displacement read = ',maxnorm
  print *

  print *,'muting displacement around the source and the edges from ',-NEX_MUTE,' to ',+NEX_MUTE
  print *

! muting displacement in elements very close to the source to avoid saturation
! also exclude elements very close to edges of the crystal
  ispec2D = 0
  elem2D_has_been_muted(:) = .false.

  do ispec2D_gamma=1,NEX_GAMMA
    do ispec2D_xi=1,NEX_XI

    ispec2D = ispec2D + 1

    if((ispec2D_xi >= NEX_XI/2 + 1 - NEX_MUTE .and. &
       ispec2D_xi <= NEX_XI/2 + NEX_MUTE .and. &
       ispec2D_gamma >= NEX_GAMMA/2 + 1 - NEX_MUTE .and. &
       ispec2D_gamma <= NEX_GAMMA/2 + NEX_MUTE) &
       .or. ispec2D_xi >= NEX_XI - NEX_MUTE + 1 &
       .or. ispec2D_gamma >= NEX_GAMMA - NEX_MUTE + 1) then

         displ2D_all_time_steps(:,ispec2D,:) = 0.
         elem2D_has_been_muted(ispec2D) = .true.

    endif

    enddo
  enddo

! print information about displacement values after muting
  maxnorm = maxval(sqrt(displ2D_all_time_steps(1,:,:)**2 + displ2D_all_time_steps(2,:,:)**2 + displ2D_all_time_steps(3,:,:)**2))

  print *
  print *,'Amplitude after muting near source and edges:'
  print *
  print *,'minval maxval displ_1 read = ',minval(displ2D_all_time_steps(1,:,:)),maxval(displ2D_all_time_steps(1,:,:))
  print *,'minval maxval displ_2 read = ',minval(displ2D_all_time_steps(2,:,:)),maxval(displ2D_all_time_steps(2,:,:))
  print *,'minval maxval displ_3 read = ',minval(displ2D_all_time_steps(3,:,:)),maxval(displ2D_all_time_steps(3,:,:))
  print *
  print *,'maxval norm displacement read = ',maxnorm
  print *


! normalize displacement globally after muting
  if(NORMALIZE_GLOBALLY) displ2D_all_time_steps(:,:,:) = displ2D_all_time_steps(:,:,:) / maxnorm

! source is located in lower-left corner of middle element
  ispec2D_xi_source = NEX_XI/2 + 1
  ispec2D_gamma_source = NEX_GAMMA/2 + 1

! create LaTeX file to assemble all the figures
  ifigure_latex = 0
  open(unit=IOUT,file='assemble_all_figures.tex',status='unknown')
  write(IOUT,*) '\documentclass[12pt]{article}'
  write(IOUT,*) '%'
  write(IOUT,*) '\usepackage[dvips]{epsfig}'
  write(IOUT,*) '%'
  write(IOUT,*) '\textwidth 17cm'
  write(IOUT,*) '\textheight 24.8cm'
  write(IOUT,*) '\oddsidemargin -4.5mm'
  write(IOUT,*) '\evensidemargin -4.5mm'
  write(IOUT,*) '\topmargin -20mm'
  write(IOUT,*) '%'
  write(IOUT,*) '\begin{document}'
  write(IOUT,*) '%'
  write(IOUT,*) '\begin{center}'
  write(IOUT,*) '\parindent 0pt'

!----
!----  main loop over all the spectral elements on vertical face Ymin
!----

  ispec2D = 0

  do ispec2D_gamma=1,NEX_GAMMA
  do ispec2D_xi=1,NEX_XI

  ispec2D = ispec2D + 1

  print *
  print *,'spectral element ',ispec2D,' out of ',NSPEC2D

! use upper-right corner only, ignore other elements
  if(ispec2D_xi <= NEX_XI/2 .or. ispec2D_gamma <= NEX_GAMMA/2) then

    print *,'element not located in upper-right corner, ignored'
    cycle

  else

    print *,'element located in upper-right corner, creating polarogram'

! ouverture fichier PostScript
  if(GENERATE_EPS) then
    write(filename,72) ispec2D_xi,ispec2D_gamma
  else
    write(filename,70) ispec2D_xi,ispec2D_gamma
  endif

 70 format('polarograms/polarogram_ixi_',i3.3,'_igamma_',i3.3,'.ps')
 72 format('polarograms/polarogram_ixi_',i3.3,'_igamma_',i3.3,'.eps')

  open(unit=24,file=filename,status='unknown')

! ecriture de l'entete du fichier PostScript
  if(GENERATE_EPS) then
    write(24,20)
  else
    write(24,10)
  endif

  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/MK {mark} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '/GG {0 setgray ST} def'
  write(24,*) '/GC {Colmesh ST} def'
  write(24,*) '/RF {setrgbcolor fill} def'
  write(24,*) '/SF {setgray fill} def'
  write(24,*) '/GS {gsave} def'
  write(24,*) '/GR {grestore} def'
  write(24,*) '% macro dessin fleche'
  write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(24,*) '%'
  write(24,*) '.01 CM setlinewidth'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizex,' CM translate 1. 1. scale'
  write(24,*) '%'
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont 1 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) ' 1. 1. scale'
  write(24,*) '%'

! definition of the time scale
  tscale = LENGTH_TIME_AXIS / dble(NSTEP)

! compute relative coordinates between receiver and source for rotation
 if(ispec2D_xi == NEX_XI/2 + 1 .and. ispec2D_gamma == NEX_GAMMA/2 + 1) then

! receiver is exactly at the source location, therefore angle is undefined
! set to zero to avoid NaN problem or division by zero
   x_diff_source = 0.d0
   z_diff_source = 0.d0
   dist_diff_source = 0.d0
   angle_diff_source_rad = 0.d0

 else
  x_diff_source = dble(ispec2D_xi - ispec2D_xi_source) * dabs(UTM_X_MAX-UTM_X_MIN)/dble(NEX_XI)
  z_diff_source = dble(ispec2D_gamma - ispec2D_gamma_source) * dabs(Z_DEPTH_BLOCK)/dble(NEX_GAMMA)

! normalize relative distance
  dist_diff_source = dsqrt(x_diff_source**2 + z_diff_source**2)
  x_diff_source = x_diff_source / dist_diff_source
  z_diff_source = z_diff_source / dist_diff_source

! compute angle between source and receiver
  angle_diff_source_rad = dacos(x_diff_source)
 endif

  angle_diff_source_deg = angle_diff_source_rad * 180.d0 / PI

! we are in upper-right corner therefore angle must be between 0. and 90. deg
  print *,'angle in deg for this receiver = ',angle_diff_source_deg
  if(angle_diff_source_deg <  -0.00001d0) stop 'angle < 0 deg detected'
  if(angle_diff_source_deg > +90.00001d0) stop 'angle > 90 deg detected'

! check if element has been muted
  if(elem2D_has_been_muted(ispec2D)) then

    write(24,*) '/Times-Roman findfont 1.6 CM scalefont setfont'
    xpos = orig_x + 0.05*LENGTH_TIME_AXIS
    zpos = orig_z
    write(24,510) xpos,zpos
    write(24,*) '(Element near source or edge muted) show'

  else

! normalize displacement locally for this trace
  if(.not. NORMALIZE_GLOBALLY) then
    maxnorm = maxval(sqrt(displ2D_all_time_steps(1,ispec2D,:)**2 + displ2D_all_time_steps(2,ispec2D,:)**2 + displ2D_all_time_steps(3,ispec2D,:)**2))
    displ2D_all_time_steps(:,ispec2D,:) = displ2D_all_time_steps(:,ispec2D,:) / maxnorm
  endif

!! DK DK perform rotation to direction of propagation
!! DK DK check if this is really what we want to do in an anisotropic medium

! draw the polarogram
  do it=1,NSTEP,ISAMP_POLAROGRAMS
    xpos = (it-1)*tscale + orig_x
    zpos = orig_z
    write(24,510) xpos,zpos

! rotate components in the (X,Z) plane, which corresponds to face Ymin
    xval = dcos(angle_diff_source_rad)*displ2D_all_time_steps(1,ispec2D,it) + &
           dsin(angle_diff_source_rad)*displ2D_all_time_steps(3,ispec2D,it)
! we are on vertical face Ymin therefore "vertical" component is along Y
    zval = displ2D_all_time_steps(2,ispec2D,it)
    xend = xval*SIZEVECT
    zend = zval*SIZEVECT
    write(24,520) xend,zend
  enddo

! ligne reliant les extremites de deux vecteurs consecutifs
  do it=1,NSTEP-ISAMP_POLAROGRAMS,ISAMP_POLAROGRAMS

! rotate components in the (X,Z) plane, which corresponds to face Ymin
    xval = dcos(angle_diff_source_rad)*displ2D_all_time_steps(1,ispec2D,it) + &
           dsin(angle_diff_source_rad)*displ2D_all_time_steps(3,ispec2D,it)
! we are on vertical face Ymin therefore "vertical" component is along Y
    zval = displ2D_all_time_steps(2,ispec2D,it)
    xpos = (it-1)*tscale + orig_x + xval*SIZEVECT
    zpos = orig_z + zval*SIZEVECT
    write(24,510) xpos,zpos

! rotate components in the (X,Z) plane, which corresponds to face Ymin
    xval = dcos(angle_diff_source_rad)*displ2D_all_time_steps(1,ispec2D,it+ISAMP_POLAROGRAMS) + &
           dsin(angle_diff_source_rad)*displ2D_all_time_steps(3,ispec2D,it+ISAMP_POLAROGRAMS)
! we are on vertical face Ymin therefore "vertical" component is along Y
    zval = displ2D_all_time_steps(2,ispec2D,it+ISAMP_POLAROGRAMS)
    xpos = (it-1+ISAMP_POLAROGRAMS)*tscale + orig_x + xval*SIZEVECT
    zpos = orig_z + zval*SIZEVECT
    write(24,530) xpos,zpos
  enddo

! draw a line to represent the time axis
  xpos = orig_x
  zpos = orig_z
  write(24,510) xpos,zpos

  xend = LENGTH_TIME_AXIS + additional_shift_arrow
  zend = 0.
  write(24,520) xend,zend

! fleche pour l'axe des temps
  xpos = LENGTH_TIME_AXIS + additional_shift_arrow + orig_x
  zpos = orig_z
  write(24,510) xpos,zpos

  xpos = - length_arrow_cm
  zpos = - length_arrow_cm
  write(24,540) xpos,zpos

  xpos = 0.d0
  zpos = + 2.d0 * length_arrow_cm
  write(24,540) xpos,zpos

  write(24,*) 'CP fill ST'

! end of test to check if element has been muted
  endif

! bounding box
  write(24,*) '1.5 CM 3 CM MV'
  write(24,*) '29.7 1.9 sub CM 3 CM L'
  write(24,*) '29.7 1.9 sub CM 15.3 CM L'
  write(24,*) '1.5 CM 15.3 CM L'
  write(24,*) 'CP ST'

! caption indicating location of receiver
  write(24,*) '/Times-Roman findfont 1 CM scalefont setfont'
  write(24,*) '2 CM 3.4 CM MV'
  if(NORMALIZE_GLOBALLY) then
    write(24,600) ispec2D_xi,ispec2D_gamma,angle_diff_source_deg
  else
    write(24,610) ispec2D_xi,ispec2D_gamma,angle_diff_source_deg
  endif

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  if(GENERATE_EPS) write(24,30)

  close(24)

! add this figure to LaTeX file to assemble all the figures
  ifigure_latex = ifigure_latex + 1
  write(IOUT,170) ispec2D_xi,ispec2D_gamma
  if(ifigure_latex < 4) then
    write(IOUT,*) '\hspace{0.5mm}%'
  else
    ifigure_latex = 0
    write(IOUT,*) '\\[0.5mm]%'
    write(IOUT,*) '%'
    write(IOUT,*) '%'
    write(IOUT,*) '%'
  endif

 170 format('\epsfig{angle=-90,file=polarograms/polarogram_ixi_',i3.3,'_igamma_',i3.3,'.eps,width=4.1cm}%')

! end of test on upper-right corner only
  endif

! end of loops on all the spectral elements on vertical face Ymin
  enddo
  enddo

! close LaTeX file to assemble all the figures
  write(IOUT,*) '\end{center}'
  write(IOUT,*) '\end{document}'
  close(IOUT)

  print *
  print *,'program ended successfully'
  print *

 10   format('%!PS-Adobe-2.0',/,'%%',/,'%% Creator: Dimitri Komatitsch',/,'%%')

 20   format('%!PS-Adobe-2.0 EPSF-1.2',/,'%%Pages: 1',/, &
        '%%DocumentFonts: Times-Roman',/,'%%BoundingBox: 162 42 514 792',/, &
        'save countdictstack mark newpath',/,'/showpage {} def',/,'%%EndProlog',/,'%%Page 1 1')

 30   format('%%Trailer',/,'cleartomark countdictstack exch sub { end } repeat restore',/,'%%EOF')

 510  format(f12.6,' CM ',f12.6,' CM M')
 520  format(f12.6,' CM ',f12.6,' CM LR ST')
 530  format(f12.6,' CM ',f12.6,' CM L ST')
 540  format(f12.6,' CM ',f12.6,' CM LR')

 600  format('(Receiver ixi=',i3,' igamma=',i3,' angle=',f9.2,' deg, norm. global) show')
 610  format('(Receiver ixi=',i3,' igamma=',i3,' angle=',f9.2,' deg, norm. local) show')

  end program plot_polarograms

