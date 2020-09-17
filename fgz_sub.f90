!
! Self contained example of calling zlib
!   David Duffy David.Duffy@qimr.edu.au 
!   version 0.1 2009-11-11
!   version 0.2 2009-11-20
!
! Definition of a port
!   slots: associated file name
!          1=uncompressed 2=gzipped 3=unzipped copy
!          Fortran style logical unit number
!          gzip C-style file handle
!
! zlib actually will read uncompressed files transparently, 
! so "1" is not actually necessary
!
module ioports
  use, intrinsic :: iso_c_binding
  type, public :: ioport
    character (len=256) :: filnam
    integer :: filtyp
    integer :: fstream
    type (c_ptr) :: gzfile = c_null_ptr
  end type ioport
end module ioports
!
! Fortran interface to zlib 
!   based on looking at fgzlib, fgsl and Janus Weil's example
!   on comp.lang.fortran May 2009
!   currently enough functionality to read gzipped text files
!
module f95zlib
  use, intrinsic :: iso_c_binding
  use ioports
  integer, parameter, public :: fgz_pathmax = 2048
  integer, parameter, public :: fgz_strmax = 128
! buffer for gzread
  integer, parameter :: ZBUFLEN = 65536
  character (len=ZBUFLEN), target :: zbuffer
! current character and end of zbuffer
  integer :: zbufpos=0, zbufend=ZBUFLEN
! gzopen  
  interface
    function gzopen(path, mode) bind(C, name='gzopen')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: path, mode
      type (c_ptr) :: gzopen
    end function
  end interface
! gzread  
  interface
    function gzread(filehandle, buf, len) bind(C, name='gzread')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzread
      type (c_ptr), value :: filehandle
      type (c_ptr), value :: buf
      integer(c_int), value :: len
    end function
  end interface
! gzgetc  
  interface
    function gzgetc(filehandle) bind(C, name='gzgetc')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzgetc
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzrewind 
  interface
    function gzrewind(filehandle) bind(C, name='gzrewind')
      use, intrinsic :: iso_c_binding
      integer(c_int) :: gzrewind
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzclose  
  interface
    function gzclose(filehandle) bind(C, name='gzclose')
    use, intrinsic :: iso_c_binding
    integer(c_int) :: gzclose
    type (c_ptr), value :: filehandle
    end function
  end interface
contains
!
! Wrapper for gzopen
!   also reinitializes gzread's buffer
!
  subroutine fgz_open(path, mode, fd, ios)
    use iso_c_binding
    character(kind=c_char, len=*), intent(in) :: path, mode
    type (ioport) :: fd
    integer :: ios

    character(kind=c_char,len=fgz_pathmax) :: lpath
    character(kind=c_char,len=fgz_strmax) :: lmode

    ios=0
    fd%filnam=path
    fd%filtyp=2
    fd%fstream=-1
    if (len(trim(path)) < fgz_pathmax .and. len(trim(mode)) < fgz_strmax) then
      lpath = trim(path) // c_null_char
      lmode = trim(mode) // c_null_char
      fd%gzfile = gzopen(lpath, lmode)
    else
      fd%gzfile = c_null_ptr
    end if
    if (.not.c_associated(fd%gzfile)) ios=-1
    zbufpos=0
  end subroutine fgz_open
!
! Wrapper for gzrewind
!
  subroutine fgz_rewind(fd, ios)
    use iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ir
    ios = 0
    ir = gzrewind(fd%gzfile)
    if (ir /= 0) ios=ir
  end subroutine fgz_rewind
!
! Wrapper for gzread
!   read one line of text from buffer
!
  subroutine fgz_read(fd, lin, advance, ios)
    use iso_c_binding
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    integer :: eos, i, j, linlen, nchar, newzpos, pos
    integer(c_int) :: blen, rlen
!
!  eol morez more
!   F    T    T    read buffer, copy to output
!   F    T    F    read buffer, output full
!   T    F    F    found <NL>
!  advancing
!   no             after output full, exit with buffer pos at end of text
!   yes            after output full, exit with buffer pos at next <NL>
!
    logical :: advancing, eol, more, morez

    type (c_ptr) :: buf = c_null_ptr

    advancing=.true.
    if (present(advance)) advancing=(advance == 'yes') 
    linlen=len(lin)
    ios=0
    lin=' '
    sta=1
    nchar=-1
    pos=0
    j=0
    eol=.false.
    more=.true.
    morez=.true.
    do while (morez)
      j=j+1
! refill buffer if necessary
      if (zbufpos == 0) then
        blen=ZBUFLEN
        buf=c_loc(zbuffer(1:1))
        rlen=gzread(fd%gzfile, buf, blen)
        if (rlen <= 0) then
          ios=-1
          return
        end if
        zbufpos=1
        zbufend=rlen
      end if
! place buffer index at <NL> or buffer end
! if <NL> will exit after updating output
      newzpos=zbufend+1
      nchar=zbufend-zbufpos+1
      do i=zbufpos, zbufend
        if (zbuffer(i:i) == achar(10)) then
          eol=.true.
          morez=.false.
          newzpos=i+1
          nchar=i-zbufpos
          exit
        end if
      end do
! read in min(buffer, remaining output)
! if not advancing move buffer idx back to last character read and exit
      if (more) then
        if (linlen < pos+nchar) then
          more=.false.
          nchar=linlen-pos
          if (.not.advancing) then
            newzpos=zbufpos+nchar
            morez=.false.
          end if
        end if
        lin((pos+1):(pos+nchar))=zbuffer(zbufpos:(zbufpos+nchar-1))
        pos=pos+nchar
      end if
      zbufpos=newzpos
      if (zbufpos > zbufend) then
        zbufpos=0
      end if
    end do
    if (.not.advancing .and. eol) ios=-2
  end subroutine fgz_read
!
! Wrapper for gzclose
!
  subroutine fgz_close(fd, ios)
    use iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ic 
    ios = 0
    ic = gzclose(fd%gzfile)
    if (ic /= 0) ios = ic
  end subroutine fgz_close
end module f95zlib
!
! Readline subroutine for either plain or gzipped files -- 
!
module fileio
  use ioports
  use f95zlib
  public :: close_infile, open_infile, readline, rewind_infile
contains
!
! Find a free Fortran style unit
!
  subroutine newlun(strm)
    integer, intent(out) :: strm
    integer :: iport
    logical :: ios
    do iport=8, 99
      inquire(iport, opened=ios)
      if (.not.ios) then
        strm=iport
        return
      end if
    end do
    write(*,'(a)') 'ERROR: No available i/o streams!'
    stop
  end subroutine newlun 
!
! Open a (plain or gzipped) file for reading
!
  subroutine open_infile(filnam, port, ios)
    character(len=*), intent(in) :: filnam
    type (ioport) :: port
    integer, intent(out) :: ios

    integer :: eon, strm
    logical :: gzipped
    character(len=len(filnam)) :: wrkfil

    ios=0
    gzipped=.false.
    eon=len_trim(filnam)
    if (eon == 0) then
      write(*,'(a)') 'ERROR: No file name given.'
      ios=1
      return
    end if
    gzipped=(filnam(max(1,(eon-2)):eon) == '.gz')
    if (gzipped) then
      call fgz_open(filnam, 'rb', port, ios)
    else
      call newlun(strm)
      open(strm, file=filnam)
      port%filnam=filnam
      port%filtyp=1
      port%fstream=strm
    end if
  end subroutine open_infile
!
! Reopen a file for reading
!
  subroutine rewind_infile(port, ios)
    type(ioport), intent(inout) :: port
    integer, intent(out) :: ios

    ios=0
    if (port%filtyp == 2) then
      call fgz_rewind(port, ios)
    else
      rewind(port%fstream)
    end if
  end subroutine rewind_infile
!
! Read one record from file
!
  subroutine readline(port, lin, advance, ios)
    type (ioport), intent(in) :: port
    character(len=*) :: lin
    character(len=*), optional :: advance
    integer, intent(out) :: ios  
    character (len=3) :: mode
    ios=0
    mode='yes'
    if (present(advance)) then
      mode=advance
    end if
    if (port%filtyp == 2) then
      call fgz_read(port, lin, advance=mode, ios=ios)
    else
      read(port%fstream,'(a)', advance=mode, iostat=ios) lin
    end if
  end subroutine readline
!
! Close file for reading -  
!
  subroutine close_infile(port, ios)
    type (ioport), intent(in) :: port
    integer, intent(out) :: ios
    ios=0
    if (port%filtyp == 1) then
      close(port%fstream)
    else if (port%filtyp == 2) then
      call fgz_close(port, ios)
    end if
  end subroutine close_infile
end module fileio


