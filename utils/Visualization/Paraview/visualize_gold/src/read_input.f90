! this subroutine reads the input information from a structured ASCII text file
! REVISION:
!   April 09,2010, Hom Nath Gharti
! FEEDBACK:
!   homnath_AT_norsar_DOT_no 
subroutine read_input (inp_fname)
use visualize_par
use string_process
implicit none

character(len=*),intent(in) :: inp_fname
character(len=256) :: line,tag
character(len=80) :: strval,token
character(len=1) :: tmp_char
character(len=80),dimension(50) :: args
integer :: ind,ios,narg,slen

integer :: i_proc,proc_count,proc_mode,slice_count1,slice_count2,slice_npmax
integer :: input_stat,output_stat,procinfo_stat,proclist_stat,server_stat
integer,dimension(3) :: proc_ind ! indicial parameters for processor list

slice_count1=0
slice_count2=0

input_stat=-1
output_stat=-1
procinfo_stat=-1
proclist_stat=-1
server_stat=-1

! default values
out_path='../output'

! open file to read
open(unit=11,file=trim(inp_fname),status='old', action='read',iostat=ios)
if (ios /= 0)then
  write(*,'(/,a)')'ERROR: input file "'//trim(inp_fname)//'" cannot be opened!'
  stop
endif

do
  read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
  if (ios/=0)exit
  
  ! check for blank and comment line
  if (isblank(line) .or. iscomment(line,'#'))cycle 
   
  ! look for line continuation
  tag=trim(line)
  call last_char(line,tmp_char,ind)
  if (tmp_char=='&')then
        slen=len(line)
        tag=trim(line(1:ind-1))
        read(11,'(a)',iostat=ios)line ! This will read a line and proceed to next line
        tag=trim(tag)//trim(line)         
  endif
  
  call first_token(tag,token)
  ! read input information
  if (trim(token)=='input:')then
    call split_string(tag,',',args,narg)
    inp_path=get_string('path',args,narg)
    dt=get_real('dt',args,narg)
    t_nstep=get_integer('nstep',args,narg)
    t_start=get_integer('start',args,narg)
    t_inc=get_integer('step',args,narg)
    t_width=get_integer('width',args,narg)
    inp_ncomp=get_integer('ncomp',args,narg); allocate(inp_head(inp_ncomp))
    if (inp_ncomp==1)then
      inp_head=get_string('head',args,narg)
    elseif (inp_ncomp==3)then ! vector
      inp_head=get_string_vect('head',inp_ncomp,args,narg)
    elseif (inp_ncomp==6)then ! tensor
      inp_head=get_string_vect('head',inp_ncomp,args,narg)
    else
      write(*,'(/,a)')'ERROR: wrong ncomp value in input: line!'
      stop
    endif
    inp_ext=get_string('ext',args,narg); inp_ext='.'//trim(inp_ext)
    dat_topo=get_integer('topo',args,narg)
    !print*,narg,arg
    input_stat=0      
    cycle      
  endif
  
  ! read output information
  if (trim(token)=='output:')then
    call split_string(tag,',',args,narg)
    !out_path=get_string('path',args,narg)     
    call seek_string('path',strval,args,narg)
    if (.not. isblank(strval))out_path=trim(strval)
    out_ncomp=get_integer('ncomp',args,narg)
    !stop
    out_head=get_string('head',args,narg)
    out_vname=get_string('vname',args,narg)       
    out_res=get_integer('res',args,narg)
    out_format=get_integer('form',args,narg)
    !print*,narg,arg
    output_stat=0      
    cycle      
  endif
  
  ! read processor information
  if (trim(token)=='procinfo:')then
    call split_string(tag,',',args,narg)      
    proc_head=get_string('head',args,narg)
    proc_width=get_integer('width',args,narg)            
    out_nslice=get_integer('nslice',args,narg)
    slice_npmax=get_integer('npmax',args,narg)
      allocate(slice_nproc(out_nslice))
      allocate(slice_proc_list(out_nslice,slice_npmax)) 
      if (out_format==1 .and. out_nslice>1)then
        ! allocate memory for server_name and server_exec
        allocate(server_name(out_nslice))
        allocate(server_exec(out_nslice))
      endif    
    procinfo_stat=0      
    cycle      
  endif
  
  ! read processor list
  if (trim(token)=='proclist:')then
    proclist_stat=-1
    call split_string(tag,',',args,narg)      
    slice_count1=slice_count1+1
    if (slice_count1>out_nslice)then
      write(*,'(/,a)')'ERROR: number of slices exceeds the actual number!'
      stop
    endif
    slice_nproc(slice_count1)=get_integer('np',args,narg)
    proc_mode=get_integer('mode',args,narg)
    if (proc_mode==0)then
      slice_proc_list(slice_count1,1:slice_nproc(slice_count1)) &
      =get_integer_vect('list',slice_nproc(slice_count1),args,narg)        
    elseif (proc_mode==1)then ! Indicial
      proc_ind=get_integer_vect('list',3,args,narg) ! start, end, step
                
      proc_count=0
      do i_proc=proc_ind(1),proc_ind(2),proc_ind(3)
        proc_count=proc_count+1
        if (proc_count>slice_npmax)then
          write(*,'(/,a)')'ERROR: number of processors per slice exceeds the maximum number!'
          stop
        endif
        slice_proc_list(slice_count1,proc_count)=i_proc
      enddo
      slice_nproc(slice_count1)=proc_count
      
    else
      write(*,'(/,a)')'ERROR: wrong proc_mode value!'
      stop
    endif            
                
    proclist_stat=0      
    cycle      
  endif
  
  ! read server information
  if (output_stat==0 .and. out_format==1 .and. out_nslice>1)then    
    if (trim(token)=='server:')then      
      server_stat=-1        
      call split_string(tag,',',args,narg)
      slice_count2=slice_count2+1
      if (slice_count2>out_nslice)then
        write(*,'(/,a)')'ERROR: number of slices exceeds the actual number!'
        stop
      endif          
      server_name(slice_count2)=get_string('name',args,narg)
      server_exec(slice_count2)=get_string('exec',args,narg)      
      server_stat=0      
      cycle      
    endif
  endif  


enddo ! do

! check for proclist: line number
if (slice_count1<out_nslice)then
  write(*,'(/,a)')'ERROR: number of lines with "proclist:" are less than the number of output slices ',out_nslice
  stop
endif

! check for server: line number
if (out_format==1 .and. out_nslice>1 .and. slice_count2>1 .and. slice_count2<out_nslice)then
  write(*,'(/,a)')'ERROR: number of lines with "server:" are less than the number of output slices ',out_nslice
  stop
endif

! set same server information if necessary
if (out_format==1 .and. slice_count2==1 .and. out_nslice>1)then
  ! set server name and executable for all other
  server_name(2:out_nslice)=server_name(1)
  server_exec(2:out_nslice)=server_exec(1)    
endif

! check input status
if (input_stat /= 0)then
  write(*,'(/,a)')'ERROR: error reading input information! make sure the line with "input:" token is correct.'
  stop
endif

! check output status
if (output_stat /= 0)then
  write(*,'(/,a)')'ERROR: error reading output information! make sure the line with "output:" token is correct.'
  stop
endif

! check procinfo status
if (procinfo_stat /= 0)then
  write(*,'(/,a)')'ERROR: error reading processor information! make sure the line with "procinfo:" token is correct.'
  stop
endif

! check proclist status
if (proclist_stat /= 0)then
  write(*,'(/,a)')'ERROR: error reading processor list! make sure the line/s with "proclist:" token is/are correct.'
  stop
endif

! check server status
if (out_format==1 .and. out_nslice>1 .and. server_stat /= 0)then
  write(*,'(/,a)')'ERROR: error reading server information! make sure the line/s with "server:" token is/are correct.'
  stop
endif

end subroutine read_input
