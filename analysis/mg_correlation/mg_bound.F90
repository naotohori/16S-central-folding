program dcd_bound_mg

   implicit none
   integer, parameter :: PREC=8
   integer, parameter :: MAXP = 1000

   integer, parameter :: FDCD = 10
   integer, parameter :: FPSF = 11
   integer, parameter :: FOUT = 20

   integer :: nmp_dcd   ! Read from DCD
   integer :: nmp_RNA   ! Read from PSF
   integer :: nmp_P     ! Read from PSF
   integer :: nmp_Mg    ! Read from PSF
   real(PREC) :: box_size, box_size_half  ! Given as arguments
   real(PREC) :: dist_cut, dist2_cut

   integer :: ix, imp, iframe, nframe
   integer :: iP, imp_P, imp_Mg
   integer :: n
   integer :: istatus
   integer :: idummy
   integer :: iarg
   integer :: iargc

   integer :: mp_P(MAXP)
   real(4), allocatable :: xyz_dcd(:,:)

   real(PREC) :: v(3)
   character(256) :: cfile_dcd, cfile_psf, cfile_out, cinp

   iarg = iargc()
   if (iarg /= 5) then
      write(*,*) 'Usage: PROGRAM [DCD file] [PSF file] [box size] [P-Mg dist to be bound] [output]'
      stop
   endif

   call getarg(1, cfile_dcd)
   call getarg(2, cfile_psf)
   call getarg(3, cinp)
   read(cinp, *) box_size
   box_size_half = 0.5 * box_size
   call getarg(4, cinp)
   read(cinp, *) dist_cut
   dist2_cut = dist_cut * dist_cut
   call getarg(5, cfile_out)

   open(FOUT, file=cfile_out, status='UNKNOWN', action='WRITE', iostat=istatus)
   if (istatus > 0) then
      write(*,*) 'Error in opening output file'
      stop
   endif

   open(FDCD, file=cfile_dcd, status='OLD', action="READ", iostat=istatus, &
              form='UNFORMATTED', access='STREAM') 
   if (istatus > 0) then
      write(*,*) 'Error in opening DCD file'
      stop
   endif

   open(FPSF, file=cfile_psf, status='OLD', action='READ', iostat=istatus)
   if (istatus > 0) then
      write(*,*) 'Error in opening PSF file'
      stop
   endif

   call read_psf(FPSF, nmp_RNA, nmp_P, nmp_Mg, mp_P)
   !write(*,*) 'nmp_RNA', nmp_RNA
   !write(*,*) 'nmp_Mg',  nmp_Mg
 
        
   call dcd_count_frame(FDCD, nframe, nmp_dcd, istatus)
   if (istatus > 0) then
      write(*,   *) 'Warning: the DCD file has an abnormal end'
      write(FOUT,*) '## Warning: the DCD file has an abnormal end'
   endif 
   !write(FOUT,*) '## Number of frames: ',nframe
   !write(FOUT,*) '## The frame number shown below starts from 0.'
   !write(FOUT,*) '## Number of MP in DCD: ',nmp_dcd

   call dcd_skip_header(FDCD)

   allocate(xyz_dcd(3,nmp_dcd))
 
   !!! iframe counting starts from 0
   do iframe = 0, nframe-1
 
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(1,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(2,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(3,imp),imp=1,nmp_dcd)
      read (FDCD) idummy

      do iP = 1, nmp_P

         imp_P = mp_P(iP)
         n = 0

         do imp_Mg = nmp_RNA + 1, nmp_RNA + nmp_Mg
         
            v(:) = xyz_dcd(:,imp_P) - xyz_dcd(:,imp_Mg)

            ! Periodic boundary
            do ix = 1, 3
               if (v(ix) > box_size_half) then
                  v(ix) = v(ix) - box_size
               else if (v(ix) < -box_size_half) then
                  v(ix) = v(ix) + box_size
               endif
            enddo

            if (dot_product(v,v) <= dist2_cut) then
               n = n + 1
            endif

         enddo

         if (n > 9) then
            write(FOUT,*) 'ERROR: n > 9 at (frame, iP)=', iframe, iP
            flush(FOUT)
            stop
         endif
         write(FOUT, '(i1)', advance='no') n
      enddo

      write(FOUT, *) ''
 
   enddo

   deallocate(xyz_dcd)

   stop

contains
   subroutine read_psf(f, nmp_RNA, nmp_P, nmp_Mg, mp_P)
      integer, intent(in) :: f
      integer, intent(out) :: nmp_RNA, nmp_P, nmp_Mg
      integer, intent(out) :: mp_P(MAXP)
      integer :: istatus

      integer :: idummy
      real(PREC) :: rdummy
      character(8) :: cdummy
      character(2) :: c
      character(1) :: t

      nmp_RNA = 0
      nmp_Mg = 0
      nmp_P = 0
      mp_P(:) = 0

      ! First line
      read (f, *, iostat=istatus) cdummy

      do
         read(f, *, iostat=istatus) cdummy, idummy, c, cdummy, t, rdummy, rdummy, idummy
         if (istatus /= 0) exit

         if (c(1:1) .eq. 'R') then
            nmp_RNA = nmp_RNA + 1
            if (t == 'P') then
               nmp_P = nmp_P + 1
               mp_P(nmp_P) = nmp_RNA
            endif
         else if (c == 'Mg') then
            nmp_Mg = nmp_Mg + 1
         else if (c == 'Cl') then
            exit
         endif

      enddo
   endsubroutine read_psf

   subroutine dcd_skip_header(f)
      integer, intent(in) :: f
      integer :: i, nblock_size, idummy

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size
   endsubroutine dcd_skip_header

   subroutine dcd_count_frame(f, n, nmp, istatus)
      integer, intent(in) :: f
      integer, intent(out) :: n, nmp, istatus
      real(4) :: xdummy
      logical :: flg_first = .True.

      call dcd_skip_header(f)

      n = 0
      istatus = 0
      L1: do 
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         if (flg_first) then
            nmp = idummy / 4
            flg_first = .False.
         endif
         Lx: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lx
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Ly: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Ly
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Lz: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lz
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         n = n + 1
      enddo L1

      rewind(f)

   endsubroutine dcd_count_frame

end program dcd_bound_mg
