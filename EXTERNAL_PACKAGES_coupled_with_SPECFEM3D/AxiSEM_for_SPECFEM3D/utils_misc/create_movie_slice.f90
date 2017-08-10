program create_movie_slice

  implicit none
  include '../../../setup/constants.h'

  character(len=256) indir,outdir,name_file,local_data_file
  character(len=256) local_file,filename,prname,prname_lp
  character(len=256) vtk_output,iteration
  character(len=1) slice_code
  integer imin,imax,istep
  real(kind=CUSTOM_REAL)  slice_coord,ycte,xcte
  integer, allocatable :: ibool(:,:,:,:),index_el_to_store(:,:),nb_index_to_run(:)
  real(kind=CUSTOM_REAL), allocatable :: xstore(:),ystore(:),zstore(:),datap(:,:,:,:),data_slice(:,:)
  real(kind=CUSTOM_REAL), allocatable :: xgrid(:),ygrid(:),zgrid(:)
  integer, allocatable :: connectiv(:,:,:)
  integer, allocatable :: NSPEC_AB_iproc(:),index_to_run(:,:)
  integer iel_to_store,ier,iglob,iglob0,iglob1,iproc
  integer nglob_ab,nproc,nspec_ab,ispec,nspec
  integer kindex,igll,jgll,kgll,nbpoints
  integer ipoint, ielm,nelem,i,it,itsn
  !integer, parameter :: IOVTK=111


  !! input files
  open(10,file='slice_movie.par')
  read(10,'(a)') indir                 !! input dir where are the spcefem solutions
  read(10,'(a)') outdir                !! output dir where will be vtk file
  read(10,'(a)') filename              !! name of the field to be plotted
  read(10,*) nproc                     !! number of proc in simu specfem3D
  read(10,*) istep,imax                !! step and end of the snapshot
  read(10,*) slice_code, slice_coord   !! slice_code=X or Y or Z and coord of the plane
                                       !! eg X 1000. will plot slice X=1000.
  close(10)
  !slice_code ! X Y or Z
  !slice_coord ! cut value

  !!
  !! 1  READ SEPCFEM3D mesh files =============================

  !! 1.1 counting the number of elements in the whole mesh
  if (slice_code == 'Y')  ycte = slice_coord
  if (slice_code == 'X')  xcte = slice_coord
  imin = istep
  NSPEC = 0

  allocate(NSPEC_AB_iproc(nproc))
  allocate(data_slice(NGLLX,NGLLZ))
  allocate(nb_index_to_run(nproc))

  do iproc=1,nproc
     write(*,'(a,i6.6,a)') trim(indir)//'proc',iproc-1,'_'
     write(name_file,'(a,i6.6,a)') trim(indir)//'proc',iproc-1,'_'
     open(unit=27,file=trim(name_file)//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
     read(27) NSPEC_AB
     read(27) NGLOB_AB
     close(27)
     NSPEC = max(NSPEC_AB,NSPEC)
     NSPEC_AB_iproc(iproc)=NSPEC_AB

  enddo

  allocate(index_el_to_store(NSPEC,nproc),index_to_run(NSPEC,nproc))
  index_el_to_store=0

  iel_to_store=0
  nbpoints=0
  nelem=0

  !! 1.2 read mesh files and select the slice ==================================
  do iproc=1,nproc

     write(name_file,'(a,i6.6,a)') trim(indir)//'proc',iproc-1,'_'
     open(unit=27,file=trim(name_file)//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
     read(27) NSPEC_AB
     read(27) NGLOB_AB

     ! ibool and global point arrays file
     allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)

     if (ier /= 0) stop 'error allocating array ibool'

     allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)

     if (ier /= 0) stop 'error allocating array xstore etc.'
     read(27) ibool
     read(27) xstore
     read(27) ystore
     read(27) zstore
     close(27)




     !! 2 detect slice to plot
     !! tableau (ipoint_to_save, iproc)  -> ispec,
     !!
     iel_to_store=0

     do ispec=1,NSPEC_AB


        select case (slice_code)

         case('Y')
            !! CAS Y=cste
            !
            ! on regarde si ycte est dans l'element

            iglob0=ibool(3,1,3,ispec)
            iglob1=ibool(3,5,3,ispec)

            if ( ycte > ystore(iglob0) .and. ycte < ystore(iglob1) ) then

              ! liste des elements a garder
              iel_to_store=iel_to_store+1
              index_el_to_store(ispec,iproc)=1 !!
              index_to_run(iel_to_store,iproc)=ispec
              nbpoints=nbpoints+NGLLX*NGLLZ
              nelem=nelem+1
            endif

           case('X')


            iglob0=ibool(1,3,3,ispec)
            iglob1=ibool(5,3,3,ispec)

            if ( xcte > xstore(iglob0) .and. xcte < xstore(iglob1) ) then

               ! liste des elements a garder
               iel_to_store=iel_to_store+1
               index_el_to_store(ispec,iproc)=1 !!
               index_to_run(iel_to_store,iproc)=ispec
               nbpoints=nbpoints+NGLLY*NGLLZ
               nelem=nelem+1
            endif

         end select

     enddo
     deallocate(xstore,ystore,zstore)
     deallocate(ibool)

     nb_index_to_run(iproc)=iel_to_store

  enddo

  allocate(xgrid(nbpoints),ygrid(nbpoints),zgrid(nbpoints))
  allocate(connectiv(NGLLX,NGLLY,nelem))

  !! 1.3 store grid point slice and connectivity ==================================
  ipoint = 0
  ielm = 0
  do iproc=1,nproc

     write(name_file,'(a,i6.6,a)') trim(indir)//'proc',iproc-1,'_'
     open(unit=27,file=trim(name_file)//'external_mesh.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
     read(27) NSPEC_AB
     read(27) NGLOB_AB

     ! ibool and global point arrays file
     allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)

     if (ier /= 0) stop 'error allocating array ibool'

     allocate(xstore(NGLOB_AB),ystore(NGLOB_AB),zstore(NGLOB_AB),stat=ier)

     if (ier /= 0) stop 'error allocating array xstore etc.'
     read(27) ibool
     read(27) xstore
     read(27) ystore
     read(27) zstore
     close(27)

    select case (slice_code)

     case('Y')
       do i=1, nb_index_to_run(iproc)
          ispec=index_to_run(i,iproc)
          ielm=ielm+1
          do kgll=1,NGLLZ
             do jgll=3,3
                do igll=1,NGLLX
                   iglob=ibool(igll,jgll,kgll,ispec)
                   ipoint=ipoint+1

                   xgrid(ipoint) = xstore(iglob)
                   ygrid(ipoint) = ystore(iglob)
                   zgrid(ipoint) = zstore(iglob)

                   connectiv(igll, kgll,ielm) = ipoint
                enddo
             enddo
          enddo
       enddo

      case('X')

        do i=1, nb_index_to_run(iproc)
          ispec=index_to_run(i,iproc)
          ielm=ielm+1
          do kgll=1,NGLLZ
             do jgll=1,NGLLY
                do igll=3,3
                   iglob=ibool(igll,jgll,kgll,ispec)
                   ipoint=ipoint+1

                   xgrid(ipoint) = xstore(iglob)
                   ygrid(ipoint) = ystore(iglob)
                   zgrid(ipoint) = zstore(iglob)

                   connectiv(jgll, kgll,ielm) = ipoint
                enddo
             enddo
          enddo
       enddo

     end select

     deallocate(xstore,ystore,zstore)
     deallocate(ibool)

  enddo

  !! 2 read field file and write vtk slice file  ========================

  !! snapshot time step loop
  do it = imin,imax,istep

     !! opening vtk file for the snapshot
     write(iteration,'(a3,i6.6)')  '_it',it
     write(vtk_output,'(a,a4)') trim(outdir)//'/'//trim(filename)//trim(iteration),'.vtk'
     write(*,*) 'open vtk',trim(vtk_output)
     open(IOVTK,file=trim(vtk_output))

     !! new file
     write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
     write(IOVTK,'(a)') 'material model VTK file'
     write(IOVTK,'(a)') 'ASCII'
     write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'

     write(IOVTK, '(a,i12,a)') 'POINTS ', nbpoints, ' float'
     !! write coord
     do i=1,nbpoints
        write(IOVTK,'(3e18.6)') xgrid(i),ygrid(i),zgrid(i)
     enddo

     write(IOVTK,*) ""
     write(IOVTK,'(a,i12,i12)') "CELLS ",nelem*(NGLLZ-1)*(NGLLX-1),nelem*(NGLLZ-1)*(NGLLX-1)*5
     !! write connectivity
     do ielm = 1, nelem
        do kgll=1,NGLLZ-1
           do igll=1,NGLLX-1
              write(IOVTK,'(5i12)') 4,connectiv(igll,kgll,ielm)-1,connectiv(igll+1,kgll,ielm)-1, &
                   connectiv(igll+1,kgll+1,ielm)-1,connectiv(igll,kgll+1,ielm)-1
           enddo
        enddo
     enddo
     write(IOVTK,*) ""
     !! read field to be stored in vtk file
     write(IOVTK,'(a,i12)') "CELL_TYPES ",nelem*(NGLLZ-1)*(NGLLX-1)
     write(IOVTK,'(6i12)') (9,itsn=1,nelem*(NGLLZ-1)*(NGLLX-1))
     write(IOVTK,*) ""
     write(IOVTK,'(a,i12)') "POINT_DATA ",nbpoints
     write(IOVTK,'(a)') "SCALARS "//trim(filename)//" float"
     write(IOVTK,'(a)') "LOOKUP_TABLE default"

     kindex=0
     do iproc=1,nproc  !!
        if (nb_index_to_run(iproc) > 0) then
          !! open specfem field
          allocate(datap(NGLLX,NGLLY,NGLLZ,NSPEC_AB_iproc(iproc)),stat=ier)
          write(prname,'(a,i6.6,a)') trim(indir)//'/'//'proc',iproc-1,'_'
          local_data_file = trim(prname)//trim(filename)//trim(iteration)//'.bin'
          write(*,*) 'reading data ',  trim(local_data_file)
          open(unit = 28,file = trim(local_data_file),status='old', &
              action='read',form ='unformatted',iostat=ier)
          read(28) datap
          close(28)


          ! read field in specfem3D output
          do i=1, nb_index_to_run(iproc) !! loop over elements meeting  the slice
             ispec=index_to_run(i,iproc) !! index of the element


            select case (slice_code)

            case('Y')
            !! we plot a 2D slice in the (XZ) plane
            do kgll=1,NGLLZ
               do jgll=3,3  !! we choose the middle of the element (Y=cte)
                  do igll=1,NGLLX

                     kindex=kindex+1
                     data_slice(igll,kgll)=datap(igll,jgll,kgll,ispec)  !! store data in the slice

                  enddo
               enddo
             enddo

             !! write field in vtk file
             do kgll=1,NGLLZ
                do igll=1,NGLLX
                   write(IOVTK,*) data_slice(igll,kgll)
                enddo
             enddo

            case('X')
              do kgll=1,NGLLZ
               do jgll=1,NGLLY  !! we choose the middle of the element (Y=cte)
                  do igll=3,3

                     kindex=kindex+1
                     data_slice(jgll,kgll)=datap(igll,jgll,kgll,ispec)  !! store
                                                                        !! data in the slice

                  enddo
               enddo
              enddo

              !! write field in vtk file
              do kgll=1,NGLLZ
                 do jgll=1,NGLLY
                    write(IOVTK,*) data_slice(jgll,kgll)
                 enddo
               enddo





           end select

          enddo

          deallocate(datap)
         endif
     enddo
     write(IOVTK,*) ""
     close(IOVTK)
  enddo

end program create_movie_slice
