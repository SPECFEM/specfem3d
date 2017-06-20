program create
  implicit none
  integer nb_tot_proc,nb_dist_by_proc,nb_frequ
  integer nb_colors,color,imin,imax,nbproc
  integer  nb_proc_by_colors,nb_frq_by_colors
  integer i,j,k

  write(*,*) 'uniform frequency distribution for MPI procs'
  write(*,*) 'total number MPI proc ?'
  read(*,*) nb_tot_proc
  write(*,*) 'subset number ?'
  read(*,*) nb_colors
  write(*,*) ' frequency number?'
  read(*,*) nb_frequ



  !write(*,*) nb_colors

  imin=0;k=0;

  nb_proc_by_colors = nb_tot_proc/nb_colors
  nb_frq_by_colors = nb_frequ/nb_colors

  open(10,file='dblepara.txt')
  write(10,*) nb_colors,nb_tot_proc
  do i=1,nb_colors
     imax=imin+nb_frq_by_colors-1
     write(10,*) i-1,imin,imax,nb_proc_by_colors
     write(*,*) i-1,imin,imax,nb_proc_by_colors
     imin=imax+1
     do  j=1,nb_proc_by_colors
       write(10,*) k
       k=k+1
     enddo
  enddo
  close(10)
end program create
