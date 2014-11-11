!----------------------------------------------------
program flexwin

use seismo_variables
implicit none

character*120, dimension (1000):: basename
character*240, dimension (1000):: obs_name, syn_name
integer :: n_seis, i
integer :: ier

! read the parameter file
call read_parameter_file()

write(*,*) 'Number of seismograms to measure : '
read(*,*) n_seis
write(*,*) n_seis

write(*,*) 'For each seismogram: '
do i = 1, n_seis
  write(*,*) 'Observed seismogram'
  read(*,'(a)') obs_name(i)
  write(*,'(a)') trim(obs_name(i))
  write(*,*) 'Synthetic seismogram'
  read(*,'(a)') syn_name(i)
  write(*,'(a)') trim(syn_name(i))
  write(*,*) 'Output basename'
  read(*,'(a)') basename(i)
  write(*,'(a)') trim(basename(i))

  if (DEBUG) write(*,*) 'DEBUG : reading sac files'
  call read_sac_files(syn_name(i),obs_name(i),ier)
  if( ier /= 0 ) then
    write(*,*) 'error files read in:',trim(syn_name(i))
    write(*,*) 'skipping files'
    cycle
  endif

  if (DEBUG) write(*,*) 'DEBUG : selecting windows'
  call select_windows_stalta2()


  if(MAKE_SEISMO_PLOTS) then
     if (DEBUG) write(*,*) 'DEBUG : writing output seismos'
     call write_seismos_gmt(basename(i))
  endif

  if(MAKE_WINDOW_FILES) then
     if (DEBUG) write(*,*) 'DEBUG : writing mt input'
     call write_mt_input_2(basename(i),obs_name(i),syn_name(i))
  endif

enddo

end program flexwin
