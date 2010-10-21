! this module contains various routines to process string
! REVISION:
!   April 09,2010, Hom Nath Gharti
! FEEDBACK:
!   homnath_AT_norsar_DOT_no
module string_process

!=====================================================

contains

! check if line is blank
logical function isblank(str)
character(len=*) :: str

isblank=.false.

if (len(trim(str))==0)isblank=.true.
return
end function isblank
!=====================================================

! check if line is blank
logical function iscomment(str,rch)
character(len=*),intent(in) :: str
character(len=1),intent(in) :: rch ! reference character for commenting
character(len=1) :: ch
integer :: ind

call first_char(str,ch,ind)

iscomment=.false.

if (ch==rch)iscomment=.true.

return
end function iscomment
!=====================================================

! get first token of a string
subroutine first_token(str,token)
implicit none
character(len=*),intent(inout) :: str
character(len=len_trim(adjustl(str))) :: tmp_str
character(len=*),intent(out) :: token
integer :: i,slen

tmp_str=trim(adjustl(str))
slen=len(tmp_str)

! set default values
token=tmp_str
str=' '

! first token is a word before first space
do i=1,slen
  if (tmp_str(i:i)==' ')then
    token=tmp_str(1:i-1)
    str=tmp_str(i+1:slen)
    exit
  endif
enddo

return

end subroutine first_token
!=====================================================

! get first non-space character of a string if any
subroutine first_char(str,ch,ind)
implicit none
character(len=*),intent(in) :: str
character(len=1),intent(out) :: ch ! first non-space character
integer,intent(out) :: ind ! index of first non-space character
integer :: i,slen

slen=len(str)

! set default values
ch=str(1:1)
ind=1

! find first character
do i=1,slen
  if (str(i:i)/=' ')then
    ch=str(i:i)
    ind=i
    exit
  endif
enddo
return
end subroutine first_char
!=====================================================

! get last non-space character of a string if any
subroutine last_char(str,ch,ind)
implicit none
character(len=*),intent(in):: str
character(len=1),intent(out) :: ch ! first non-space character
integer,intent(out) :: ind ! index of first non-space character
integer :: i,slen

slen=len(str)

! set default values
ch=str(slen:slen)
ind=slen

! find last character
do i=slen,1,-1
  if (str(i:i)/=' ')then
    ch=str(i:i)
    ind=i
    exit
  endif
enddo
return
end subroutine last_char
!=====================================================

! split string with a delimter (delm) into several parts
subroutine split_string(str,delm,args,narg)
implicit none
character(len=*),intent(in) :: str
character(len=1),intent(in) :: delm ! delimeter to split string
character(len=*),dimension(*),intent(out) :: args ! list of slitted strings
integer,intent(out) :: narg

character(len=len_trim(str)) :: tmp_str
integer :: i,i1,slen
integer,dimension(100) :: ind

slen=len_trim(str)
tmp_str=trim(str)

! find and count indices of all delimeters 
narg=0
do i=1,slen
  if(tmp_str(i:i)==delm)then
    narg=narg+1
    ind(narg)=i
  endif
enddo
narg=narg+1
ind(narg)=slen+1

! split string and set to args
i1=1
do i=1,narg
  args(i)=tmp_str(i1:ind(i)-1)
  i1=ind(i)+1
enddo

return
end subroutine split_string
!=====================================================

! get string value from string list which contain a character '=' that separates
! variable name and variable vlue 
character(len=80) function get_string(vname,slist,nvar)
character(len=*),intent(in) :: vname
character(len=*),dimension(*) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)get_string     
     return
  endif
enddo

write(*,'(/,a)')'ERROR: cannot read the string variable "'//vname//'"!'
stop
end function get_string
!=====================================================

! seek string value from string list which contain a character '=' that separates
! variable name and variable vlue
subroutine seek_string(vname,strval,slist,nvar)
character(len=*),intent(in) :: vname
character(len=*),intent(out) :: strval
character(len=*),dimension(*) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

strval=''
do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)strval
     return
  endif
enddo

return
!write(*,'(/,a)')'ERROR: cannot read the string variable "'//vname//'"!'
!stop
end subroutine seek_string
!=====================================================

! get string vector from string list which contain a character '=' that separates
! variable name and variable vlue
function get_string_vect(vname,n,slist,nvar)
implicit none
character(len=80),dimension(n) :: get_string_vect
character(len=*),intent(in) :: vname
integer,intent(in) :: n
character(len=*),dimension(*),intent(in) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)get_string_vect
     return
  endif
enddo

write(*,'(/,a)')'ERROR: cannot read the integer variable "'//vname//'"!'
stop
end function get_string_vect
!=====================================================

! get integer value from string list which contain a character '=' that separates
! variable name and variable vlue
integer function get_integer(vname,slist,nvar)
character(len=*),intent(in) :: vname
character(len=*),dimension(*),intent(in) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)get_integer
     return
  endif
enddo

write(*,'(/,a)')'ERROR: cannot read the integer variable "'//vname//'"!'
stop
end function get_integer
!=====================================================

! get integer vector from string list which contain a character '=' that separates
! variable name and variable vlue
function get_integer_vect(vname,n,slist,nvar)
implicit none
integer :: n
integer,dimension(n) :: get_integer_vect
character(len=*),intent(in) :: vname
character(len=*),dimension(*),intent(in) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)get_integer_vect(1:n)
     return
  endif
enddo

write(*,'(/,a)')'ERROR: cannot read the integer variable "'//vname//'"!'
stop
end function get_integer_vect
!=====================================================

! get real value from string list which contain a character '=' that separates
! variable name and variable vlue
real function get_real(vname,slist,nvar)
character(len=*),intent(in) :: vname
character(len=*),dimension(*),intent(in) :: slist
integer,intent(in) :: nvar
character(len=80),dimension(2) :: args 
integer :: i,narg

do i=1,nvar
  call split_string(slist(i),'=',args,narg)
  if (narg/=2)cycle
  if (vname==trim(adjustl(args(1))))then
     read(args(2),*)get_real
     return
  endif
enddo

write(*,'(/,a)')'ERROR: cannot read the real variable "'//vname//'"!'
stop
end function get_real
!=====================================================

end module string_process

