!
! This utility reorders atoms in the [Atoms] section (imod=0) or in the
! [GTO] section (imod=1).
! imod=0 has been tested only for MOLPRO's MOLDEN file, whereas
! imod=1 has been tested only for CFour's MOLDEN file
! 
! How to use this program?
!
! 1) $F90 -O3 ReOrdAtm.f90 -o roa.exe
! 2) roa.exe -m {imod} < old_molden > new_molden
!    where {imod} can be 0 (default) or 1.
!
program ReOrdAtm

implicit none
integer(kind=4) :: iinp=5, iout=6, igto=71, i, NAtm=0, imod=0
integer(kind=4),allocatable :: iord(:)
character*100 :: ctmp
character*8 :: starint

! read argument "-m {imod}"
call ReadInp(ctmp,imod)

! #atoms
call CountAtm(iinp,iout,NAtm,ctmp)

! save new MOLDEN file
write(iout,"('[Molden Format]')")
if(imod .eq. 1) write(iout,"('[Program] CFour')")

call searchar(iinp,iout,7,"[ATOMS]",ctmp)
! NOTE: some programs do not recognize ANGS; Angs must be used instead
if(index(ctmp,'ANGS').ne.0)then
  write(iout,"('[Atoms] Angs')")
else
  write(iout,"('[Atoms] AU')")
end if

if(imod .eq. 0)then
  allocate(iord(NAtm))
  ! reference atomic ordering in [GTO]
  call RefOrd(iinp,iout,NAtm,iord,ctmp)
  do i=1,NAtm
    call RdAtmI(iinp,iout,NAtm,iord(i),i,ctmp)
  end do
else if(imod .eq. 1)then
  call RdAtm(iinp,iout,NAtm,ctmp)
end if

write(iout,"('[GTO]')")
call searchar(iinp,iout,5,"[GTO]",ctmp)
if(imod .eq. 0)then
  i=0
  do while(.true.)
    read(iinp,"(100a)")ctmp
    if(len_trim(ctmp) .eq. 0) cycle
    i = i +1
    write(iout,"(i4,'  0')")i
  
!   search the next empty line
    do while(.true.)
      read(iinp,"(100a)")ctmp
      write(iout,"(a)")trim(ctmp)
      if(len_trim(ctmp) .eq. 0) exit
    end do
    if(i .eq. NAtm) exit
  end do
else if(imod .eq. 1)then
  open(igto,file='gto123456789.tmp')
  call gtocpy(NAtm,iinp,igto)
  do i=1,NAtm
  	write(starint,"('***',i5)") i
  	call searchar(igto,iout,8,starint,ctmp)
    write(iout,"(i4,'  0')")i
    do while(.true.)
      read(igto,"(100a)")ctmp
      if(index(ctmp,"***") .eq. 1) exit
      write(iout,"(a)")trim(ctmp)
    end do
    write(iout,*)
  end do
  close(igto,status='delete')
end if

write(iout,"('[MO]')")
call searchar(iinp,iout,4,"[MO]",ctmp)
do while(.true.)
  read(iinp,"(100a)",err=100,end=100)ctmp
  if(index(ctmp,'[') .ne. 0 .and. index(ctmp,']') .ne. 0) exit
  write(iout,"(a)")trim(ctmp)
end do

100   stop
end

!
! read & write coordinates
!
subroutine RdAtm(iinp,iout,NAtm,ctmp)

implicit none
integer(kind=4) :: iinp, iout, NAtm, i, ia, iz
real(kind=8) :: xyz(3)
character*100 :: ctmp

call searchar(iinp,iout,7,"[ATOMS]",ctmp)

do i=1,NAtm
  read(iinp,*)ctmp, ia, iz, xyz(1), xyz(2), xyz(3)
  if(ia .ne. i)goto 100
  write(iout,"(a4,2i5,3f20.10)")trim(ctmp), i, iz, xyz
end do
return

100   write(iout,"(' Error! This MOLDEN file is not supported!')")
stop

return
end

!
! read coordinates of the iatm-th atom
! note that the atoms may be not ordered ascendingly.
!
subroutine RdAtmI(iinp,iout,NAtm,iold,inew,ctmp)

implicit none
integer(kind=4) :: iinp, iout, NAtm, iold, inew, i, ia, iz
real(kind=8) :: xyz(3)
character*100 :: ctmp

call searchar(iinp,iout,7,"[ATOMS]",ctmp)

do i=1,NAtm
  read(iinp,*)ctmp, ia, iz, xyz(1), xyz(2), xyz(3)
  if(ia .eq. iold)then
    write(iout,"(a4,2i5,3f20.10)")trim(ctmp), inew, iz, xyz
    return
  end if
end do

return
end

! make a copy of basis functions (imod = 1)
Subroutine gtocpy(NAtm,iinp,igto)

implicit none
integer(kind=4) :: i,ia,NAtm,iinp,igto
character*100 :: ctmp

i=0
do while(.true.)
  read(iinp,"(100a)")ctmp
  if(len_trim(ctmp) .eq. 0) cycle
  i = i +1
  read(ctmp,*) ia
  write(igto,"('***',i5)") ia

! search the next empty line
  do while(.true.)
    read(iinp,"(100a)")ctmp
    if(len_trim(ctmp) .eq. 0) exit
    write(igto,"(a)")trim(ctmp)
  end do
  if(i .eq. NAtm) exit
end do
write(igto,"('***',i5)") 0

Return
End

!
! atomic ordering in [GTO] (imod = 0)
!
subroutine RefOrd(iinp,iout,NAtm,iord,ctmp)

implicit none
integer(kind=4) :: iinp, iout, NAtm, iord(*), iatm=0
character*100 :: ctmp

call searchar(iinp,iout,5,"[GTO]",ctmp)

do while(.true.)
  read(iinp,"(100a)",err=100,end=100)ctmp
  if(len_trim(ctmp) .eq. 0) cycle
  iatm = iatm +1
  read(ctmp,*) iord(iatm)
  if(iatm .eq. NAtm) return

! search the next empty line
  do while(.true.)
    read(iinp,"(100a)",err=100,end=100)ctmp
    if(len_trim(ctmp) .eq. 0) exit
  end do
end do

100   write(iout,"(' Error when reading [GTO]!')")
stop

return
end

! read input argument "-m {imod}"
subroutine ReadInp(ctmp,imod)

implicit none
integer(kind=4) :: i, istr, iend, imod, nonspace
character*100 :: ctmp

i = 0
do while(.true.)
  i = i + 1
  call GetCmdUp(i,ctmp)
  istr=nonspace(ctmp)
  iend=len_trim(ctmp)

  if(iend .eq. 0) then
    Return
  else if(ctmp(istr:iend) .eq. '-M') then
    i = i + 1
    call GetCmdUp(i,ctmp)
    istr=nonspace(ctmp)
    iend=len_trim(ctmp)
    if(ctmp(istr:iend) .eq. '1') imod = 1
  else
    write(*,"(/,' Unknown argument: ',a,/)")ctmp(istr:iend)
    stop
  end if
  cycle
end do

Return
end

! read an argument in upper case
Subroutine GetCmdUp(i,ctmp)

implicit none
integer(kind=4) :: i
character*100 :: ctmp

call get_command_argument(i,ctmp)
! call getarg(i,ctmp)
call charl2u(ctmp,len_trim(ctmp))

Return
End

! position of the first non-space character in a string.
function nonspace(string)

implicit none
integer(kind=4) :: length, i, nonspace
character*(*) :: string

length=LEN_TRIM(string)
if(length .le. 1) then
  i=length
else
  do i=1,length
    if(string(i:i) .ne. ' ') exit
  end do
endif

nonspace=i

return
end

!
! search [Atoms] and count the number of atoms
!
subroutine CountAtm(iinp,iout,NAtm,ctmp)

implicit none
integer(kind=4) :: iinp, iout, NAtm
character*100 :: ctmp

call searchar(iinp,iout,7,"[ATOMS]",ctmp)

do while(.true.)
  read(iinp,"(100a)",err=100,end=100)ctmp
  if(index(ctmp,'[').ne.0 .and. index(ctmp,']').ne.0) exit
  if(len_trim(ctmp) .eq. 0) cycle
  NAtm = NAtm +1
end do

100   if(NAtm .lt. 1)then
  write(iout,"(' Error! NAtm < 1.')")
  stop
end if

return
end

!
! search LOGO
!
subroutine searchar(iinp,iout,lenth,LOGO,ctmp)

implicit none
integer(kind=4) :: iinp, iout, lenth
character*(*) :: LOGO
character*100 :: ctmp

rewind(iinp)

do while(.true.)
  read(iinp,"(100a)",err=100,end=100)ctmp
  call charl2u(ctmp,len_trim(ctmp))
  if(index(ctmp,LOGO(1:lenth)).ne.0) return
end do

100   write(iout,"(' Error! ',a,/,' was not found!')")trim(LOGO(1:lenth))
stop

return
end

!
! tmp --> TMP
!
subroutine charl2u(tmp,nc)

implicit none
integer(kind=4) :: i, nc
character*(*) :: tmp
character*1,external :: L2U

do i=1,nc
  tmp(i:i)=L2U(tmp(i:i))
end do

return
end

!
! l-->L
!
function L2U(letter)

implicit none
character*1 letter,L2U

if((ichar(letter).ge.97).and.(ichar(letter).le.122))then
  L2U=char(ichar(letter)-32)
else
  L2U=letter
endif

return
end

