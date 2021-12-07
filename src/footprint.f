      PROGRAM reformat
      
c     ***********************************************************************
c     * Change format of a trajectory file                                  *
c     * Michael Sprenger / Spring, summer 2010                              *
c     ***********************************************************************

      implicit none
      
c     ----------------------------------------------------------------------
c     Declaration of variables
c     ----------------------------------------------------------------------

c     Mode & Special parameters depending on mode
      character*80                           mode
      real                                   clon,clat,radius

c     Input and output files
      character*80                           inpfile     ! Input filename
      character*80                           outfile     ! Output filename

c     Trajectories
      integer                                ntra        ! Number of trajectories
      integer                                ntim        ! Number of times
      integer                                ncol        ! Number of columns
      real,allocatable, dimension (:,:,:) :: tra         ! Trajectories (ntra,ntim,ncol)
      real,allocatable, dimension (:,:)   :: proxy       ! Footprint
      character*80                           vars(100)   ! Variable names
      integer                                refdate(6)  ! Reference date

c     Auxiliary variables
      integer                                inpmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      real                                   dist
      real                                   lon0,lat0,lon1,lat1
      real                                   boost
      real                                   hhmm,tfrac
      integer                                nout

c     Externals
      real                                   sdis
      external                               sdis

c     ----------------------------------------------------------------------
c     Preparations
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='footprint.param')
       read(10,*) inpfile
       read(10,*) outfile
       read(10,*) ntra,ntim,ncol
       read(10,*) mode
       if ( mode.eq.'proxy' ) then
          read(10,*) clon
          read(10,*) clat
          read(10,*) radius
       endif
      close(10)
      
c     Check that a valid mode is selected
      if ( mode.eq.'proxy' ) goto 10
      print*,' Unknown mode ',trim(mode)
      stop
 10   continue

c     Determine the formats
      call mode_tra(inpmode,inpfile)
      if (inpmode.eq.-1) inpmode=1

c     Allocate memory
      allocate(tra(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tra      ***' 
      if ( mode.eq. 'proxy' ) then
           allocate(proxy(ntra,ncol+1),stat=stat)
           if (stat.ne.0) print*,'*** error allocating array proxy  ***'  
      endif

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)
      
c     ----------------------------------------------------------------------
c     Mode = 'proxy': get nearest distance to a lat/lon point
c     ----------------------------------------------------------------------

      if ( mode.ne.'proxy') goto 200

c     Transform into fractional time
      do i=1,ntra
         do j=1,ntim
            hhmm = tra(i,j,1)
            call hhmm2frac(hhmm,tfrac)
            tra(i,j,1) = tfrac
         enddo
      enddo

c     Get proxy and get number of output steps
      call get_proxy (proxy,clon,clat,tra,ntra,ntim,ncol,radius)
      
c     Transform into hhmm time       
      do i=1,ntra
         tfrac = proxy(i,1)
         call frac2hhmm(tfrac,hhmm)
         proxy(i,1) = hhmm
      enddo
      
c     Write output
      vars(ncol+1) = 'dist'
      open(10,file=outfile)
         call write_hea(10,refdate,vars,1,1,ncol+1,1)
         do i=1,ntra
            if ( proxy(i,ncol+1).lt.radius ) then
                write(10,'(1f7.2,f9.2,f8.2,i6,100f10.3)') 
     >               (proxy(i,j),j=1,3),             ! time, lon, lat
     >               nint(proxy(i,4)),               ! p
     >               (proxy(i,j),j=5,ncol+1)         ! fields
           endif
        enddo
      close(10)
   
 200  continue
 
      end


c     --------------------------------------------------------------------
c     Subroutines to write the netcdf output file
c     --------------------------------------------------------------------

      subroutine writecdf2D(cdfname,
     >                      varname,arr,time,
     >                      dx,dy,xmin,ymin,nx,ny,
     >                      crefile,crevar)

c     Create and write to the netcdf file <cdfname>. The variable
c     with name <varname> and with time <time> is written. The data
c     are in the two-dimensional array <arr>. The list <dx,dy,xmin,
c     ymin,nx,ny> specifies the output grid. The flags <crefile> and
c     <crevar> determine whether the file and/or the variable should
c     be created

      IMPLICIT NONE

c     Declaration of input parameters
      character*80 cdfname,varname
      integer nx,ny
      real arr(nx,ny)
      real dx,dy,xmin,ymin
      real time
      integer crefile,crevar

c     Further variables
      real varmin(4),varmax(4),stag(4)
      integer ierr,cdfid,ndim,vardim(4)
      character*80 cstname
      real mdv
      integer datar(14),stdate(5)
      integer i

C     Definitions
      varmin(1)=xmin
      varmin(2)=ymin
      varmin(3)=1050.
      varmax(1)=xmin+real(nx-1)*dx
      varmax(2)=ymin+real(ny-1)*dy
      varmax(3)=1050.
      cstname=trim(cdfname)//'_cst'
      ndim=4
      vardim(1)=nx
      vardim(2)=ny
      vardim(3)=1
      stag(1)=0.
      stag(2)=0.
      stag(3)=0.
      mdv=-999.98999

C     Create the file
      if (crefile.eq.0) then
         call cdfwopn(cdfname,cdfid,ierr)
         if (ierr.ne.0) goto 906
      else if (crefile.eq.1) then
         call crecdf(cdfname,cdfid,
     >        varmin,varmax,ndim,cstname,ierr)
         if (ierr.ne.0) goto 903

C        Write the constants file
         datar(1)=vardim(1)
         datar(2)=vardim(2)
         datar(3)=int(1000.*varmax(2))
         datar(4)=int(1000.*varmin(1))
         datar(5)=int(1000.*varmin(2))
         datar(6)=int(1000.*varmax(1))
         datar(7)=int(1000.*dx)
         datar(8)=int(1000.*dy)
         datar(9)=1
         datar(10)=1
         datar(11)=0            ! data version
         datar(12)=0            ! cstfile version
         datar(13)=0            ! longitude of pole
         datar(14)=90000        ! latitude of pole
         do i=1,5
            stdate(i)=0
         enddo
c         call wricst(cstname,datar,0.,0.,0.,0.,stdate)
      endif

c     Write the data to the netcdf file, and close the file
      if (crevar.eq.1) then
         call putdef(cdfid,varname,ndim,mdv,
     >        vardim,varmin,varmax,stag,ierr)
         if (ierr.ne.0) goto 904
      endif
      call putdat(cdfid,varname,time,0,arr,ierr)
      if (ierr.ne.0) goto 905
      call clscdf(cdfid,ierr)

      return

c     Error handling
 903  print*,'*** Problem to create netcdf file ***'
      stop
 904  print*,'*** Problem to write definition ***'
      stop
 905  print*,'*** Problem to write data ***'
      stop
 906  print*,'*** Problem to open netcdf file ***'
      stop

      END



