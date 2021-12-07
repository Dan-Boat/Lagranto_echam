      PROGRAM caltra

C     ********************************************************************
C     *                                                                  *
C     * Calculates trajectories                                          *
C     *                                                                  *
C     *	Heini Wernli	   first version:	April 1993                   *
C     * Michael Sprenger   major upgrade:       2008-2009                *
c     * Michael Sprenger   on-the-fly version   2016 / 2017              *
C     *                                                                  *
C     ********************************************************************

      implicit none
      
c     --------------------------------------------------------------------
c     Declaration of parameters
c     --------------------------------------------------------------------

c     Maximum number of levels for input files
      integer   nlevmax
      parameter	(nlevmax=100)

c     Maximum number of input files (dates, length of trajectories)
      integer   ndatmax
      parameter	(ndatmax=5000)

c     Numerical epsilon (for float comparison)
      real      eps
      parameter (eps=0.001)

c     Distance in m between 2 lat circles 
      real	deltay
      parameter	(deltay=1.112E5)

c     Number of iterations for iterative Euler scheme
      integer   numit
      parameter (numit=3)

c     Jump flag (=1: trajectory moving into topography will be lifted)
      integer   jflag
      parameter (jflag = 1)

c     Timecheck ('yes' enforces an explicit check of the time on P files)
      character*80 timecheck
      parameter    (timecheck='no')

c     Filename prefix (typically 'P')
      character*1 prefix
      parameter   (prefix='P')
      
c     Scaling factor for vertical wind speed (must be in [Pa/s])
      real        wfactor
      parameter   (wfactor=1.)
      
c     Name of the tracevar file
      character*80 tracevar
      parameter    (tracevar = 'tracevars' )

c     --------------------------------------------------------------------
c     Declaration of variables
c     --------------------------------------------------------------------

c     Input parameters
      integer                                fbflag          ! Flag for forward/backward mode
      integer                                aglflag         ! Flag for starting positions (1 = AGL, 0 = absolute)
      integer                                nearestflag     ! Nearest-neighbor interpolation
      integer                                numdat          ! Number of input files
      character*13                           dat(ndatmax)    ! Dates of input files
      integer                                itime(ndatmax)  ! Time corresponding to dat (in min)
      integer                                timeinc         ! Time increment between input files (in minutes)
      integer                                ntra            ! Number of trajectories
      character*80                           cdfname         ! Name of output files
      integer                                ts              ! Time step
      integer                                deltout         ! Output time interval (in minutes)
      character*80                           strname         ! File with start positions
      character*80                           caltra_dir      ! base directory for caltra input netCDF files
      character*80                           caltra_struct   ! structure of input netCDFs

c     Trajectories
      integer                                ncol            ! Number of columns for insput trajectories
      real,allocatable, dimension (:,:,:) :: traout          ! Output trajectories (ntra,ntim,4+ntrace)
      integer                                reftime(6)      ! Reference date
      character*80                           vars(200)       ! Field names
      real,allocatable, dimension (:)     :: xx0,yy0,pp0     ! Position of air parcels
      integer,allocatable, dimension (:)  :: tt0             ! Time of air parcel
      integer,allocatable, dimension(:)   :: active          ! Active (1) or inactive (0) trajectory
      integer,allocatable, dimension (:)  :: leftflag        ! Flag for domain-leaving
      real                                   xx1,yy1,pp1     ! Updated position of air parcel
      integer                                leftcount       ! Number of domain leaving trajectories
      integer                                ntim            ! Number of output time steps
      integer                                timerange       ! Timerange (in min) of trajectories

c     Input and output format for trajectories (see iotra.f)
      integer   inpmode
      integer   outmode

c     Meteorological fields
      real,allocatable, dimension (:)     :: spt0,spt1       ! Surface pressure
      real,allocatable, dimension (:)     :: uut0,uut1       ! Zonal wind
      real,allocatable, dimension (:)     :: vvt0,vvt1       ! Meridional wind
      real,allocatable, dimension (:)     :: wwt0,wwt1       ! Vertical wind
      real,allocatable, dimension (:)     :: p3t0,p3t1       ! 3d-pressure 
      real,allocatable, dimension (:)     :: fld0,fld1       ! 3d tracing field 

      real                                   pollon,pollat   ! Longitude/latitude of pole
      real                                   polgam          ! Rotation of grid
      real                                   xmin,xmax       ! Zonal grid extension
      real                                   ymin,ymax       ! Meridional grid extension
      integer                                nx,ny,nz        ! Grid dimensions
      real                                   dx,dy           ! Horizontal grid resolution
      real                                   mdv             ! Missing data value
      real                                   per             ! Periodic domain
      integer                                hem             ! Hemispheric flag
      real                                   ak(nlevmax)     ! Vertical layers and levels
      real                                   bk(nlevmax) 

c     Variable tracing
      integer                                ntrace
      character*80                           tvar(25)        ! Name of tracing variable
      real                                   tfac(25)        ! Scaling factor for output
      character*80                           tfil(25)        ! File name prefix

c     Auxiliary variables                 
      real                                   rd
      integer	                             itm,iloop,i,j,k,filo,lalo
      integer                                ierr,stat
      integer                                cdfid,fid
      integer	                             time0,time1,time
      real                                   reltpos0,reltpos1
      real                                   xind,yind,pind,pp,sp,stagz
      character*80                           filename,varname
      integer                                reftmp(6)
      character                              ch
      real                                   frac,tload
      integer                                itim
      real                                   lon,lat
      integer                                timemin,timemax
      character*200                          arg
      integer                                first
      real                                   thhmm,tfrac
      integer                                date1(5),date2(5)
      logical                                file_exists
      integer                                nactive,delta
      integer                                indt,indt0,indt1
      integer                                file_load
      real                                   x0,y0,p0,f0,zind,x1,y1
      character*80                           datestr,datformat
      real                                   rtimerange
      character*80                           yyyy,mm

c     Externals 
      real                                   int_index4
      external                               int_index4
      integer                                iargc

c     --------------------------------------------------------------------
c     Start of program, Read parameters
c     --------------------------------------------------------------------

c     Write start message
      print*,'========================================================='
      print*,'              *** START OF PROGRAM CALTRA ***'
      print*

c     Usage
      if ( iargc().eq.0 ) then
        print*,'caltra {startf} {timerange in [h]} {outfile} [optional]'
        print*
        print*,
     >    ' -i value       : time increents (in min) between files'
        print*,
     >    ' -ts value      : time step (in min) [default = timeinc/12 ]'
        print*,
     >    ' -o value       : output interval (in min) [default = ts ]'
        print*,
     >    ' -ref date      : reference date {year}{month}{day}_{hour}'
        print*,
     >    ' -cdf dir struct: base dir and structure (flat,yyyy.mm)'
        print*,
     >    ' -agl           : starting pressures are AGL'
         print*,
     >    ' -nearest       : nearest-neighbor interpolation in trace'
        stop
      endif

c     Get all mandatory arguments
      call getarg(1,arg)
      strname = trim(arg)
      call getarg(2,arg)
      read(arg,*) rtimerange
      call getarg(3,arg)
      cdfname = trim(arg)

c     Get optional arguments
      i             = 4
      timeinc       = 60
      ts            = 5.
      deltout       = 5
      first         = 1
      aglflag       = 0
      nearestflag   = 0
      reftime(1)    = 0
      reftime(2)    = 0
      reftime(3)    = 0
      reftime(4)    = 0
      reftime(5)    = 0
      reftime(6)    = 0
      caltra_dir    = ''
      caltra_struct = 'flat'
      do while ( i.le.iargc() )
        
        call getarg(i,arg)
        
        if ( trim(arg).eq.'-i' ) then
           if ( first.eq.1 ) then
            print*,
     >      '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           call getarg(i+1,arg)
           i = i + 2
           read(arg,*) timeinc
           print*,'  optional argument -i set to  : ',timeinc, ' min'
        endif
        
        if ( trim(arg).eq.'-ts' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           call getarg(i+1,arg)
           i = i + 2
           read(arg,*) ts
           print*,'  optional argument -ts set to : ',ts, ' min'
        endif
        
        if ( trim(arg).eq.'-agl' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           i = i + 1
           aglflag = 1
           print*,'  optional argument -agl set to : ',aglflag
        endif
               
        if ( trim(arg).eq.'-nearest' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           i = i + 1
           nearestflag = 1
           print*,'  optional argument -nearest set to : ',nearestflag
        endif
        
        if ( trim(arg).eq.'-o' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           call getarg(i+1,arg)
           i = i + 2
           read(arg,*) deltout
           print*,'  optional argument -o set to  : ',deltout, ' min'
        endif
        
        if ( trim(arg).eq.'-ref' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           call getarg(i+1,arg)
           i = i + 2
           read(arg(1:4),*)   reftime(1)  ! year
           read(arg(5:6),*)   reftime(2)  ! month
           read(arg(7:8),*)   reftime(3)  ! day
           read(arg(10:11),*) reftime(4)  ! hour
           if ( len_trim(arg).eq.13 ) then
              read(arg(12:13),*) reftime(5)  ! min
              datformat = 'yyyymmdd_hhmm'
           else 
              reftime(5) = 0
              datformat = 'yyyymmdd_hh'
           endif
           print*,'  optional argument -ref set to  : ',trim(arg)
        endif
              
        if ( trim(arg).eq.'-cdf' ) then
           if ( first.eq.1 ) then
            print*,
     >       '---- OPTIONAL PARAMETERS --------------------------------'
              print*
              first = 0
           endif
           call getarg(i+1,caltra_dir)
           call getarg(i+2,caltra_struct)
           i = i + 3
           print*,'  optional argument -cdf set to  : ',
     >               trim(caltra_dir),' + ',trim(caltra_struct)
        endif
      enddo
      if ( first.eq.0 ) print*
      
c     Convert time range from h to min
      call hhmm2frac(rtimerange,tfrac)
      timerange = nint(60. * tfrac)
      

c     Check that output interval is consistent with time step
      if ( mod(deltout,ts).ne.0 ) then
        print*,' ERROR: output interval must be multiple of timestep'
        stop
      endif
      
c     Check that timerange is consistent with output interval
      if ( mod(timerange,deltout).ne.0 ) then
        print*,' ERROR: timerange must be multiple of output interval'
        stop
      endif

c     Set flag for forward/backward mode
      if ( timerange.gt.0 ) then
         fbflag = 1
      else
         fbflag = -1
      endif
      
c     Set the formats of the input and output files
      call mode_tra(outmode,cdfname)
      if (outmode.eq.-1) outmode=1

c     Write some status information
      print*,'---- INPUT PARAMETERS -----------------------------------'
      print* 
      print*,'  Forward/Backward       : ',fbflag
      print*,'  time increment         : ',timeinc
      print*,'  Output file            : ',trim(cdfname)
      print*,'  Time step (min)        : ',ts
      print*,'  Output time interval   : ',deltout
      print*,'  Jump flag              : ',jflag
      print*,'  Vertical wind (scale)  : ',wfactor
      print*,'  Trajectory pos file    : ',trim(strname)
      print*,'  Input format           : (time,lon,lat,p)-list'
      print*,'  Output format          : ',outmode
      print*,'  Time check             : ',trim(timecheck)
      print*,'  Date format            : ',trim(datformat)
      print*
      
c     Check whether a tracing file exists; if yes, read it
      INQUIRE(FILE=tracevar, EXIST=file_exists)
      if ( file_exists.eqv..false. ) then
         ntrace = 0
      else
         ntrace = 0
         open(10,file=tracevar)
  10         continue
             ntrace = ntrace + 1
             read(10,*,end=20) tvar(ntrace), tfac(ntrace), tfil(ntrace)
             goto 10
  20         continue
             ntrace = ntrace -1            
         close(10)
      endif
      
c     Write status info on screen
      if ( ntrace.gt.0 ) then
         print*,
     >       '---- TRACING --------------------------------------------'
         print* 
         print*,' tracing file      : ',trim(tracevar)
         print*
         do i=1,ntrace
            write(*,'(2x,a10,f10.2,2x,a5)') tvar(i),tfac(i),tfil(i)
         enddo
         print*
      endif
     
c     --------------------------------------------------------------------
c     Get starting positions, allocate memory for trajectories
c     --------------------------------------------------------------------

c     Get the number of starting positions
      open(10,file=strname)
         ntra = 0
 30      ntra = ntra + 1
         read(10,*,end=35) thhmm,x0,y0,p0
         goto 30
 35   close(10)
      ntra = ntra - 1
      ncol = 4 + ntrace
      ntim = abs(timerange/deltout) + 1  
         
C     Get memory for trajectory arrays
      allocate(traout(ntra,ntim,ncol),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array traout   ***' ! Output trajectories
      allocate(xx0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array xx0      ***' ! X position (longitude)
      allocate(yy0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array yy0      ***' ! Y position (latitude)
      allocate(pp0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array pp0      ***' ! Pressure
      allocate(tt0(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tt0      ***' ! Time
      allocate(leftflag(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array leftflag ***' ! Leaving-domain flag
      allocate(active(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array active ***'   ! Active/inactive trajectories

c     Read start coordinates from file - Format (time,lon,lat,lev)
      open(10,file=strname)
          do i=1,ntra
             read(10,*) thhmm,xx0(i),yy0(i),pp0(i)
             call hhmm2frac(thhmm,tfrac)
             tt0(i) = nint( 60. * tfrac )
          enddo
      close(10)
      
c     Check that a reference date is set
      if ( reftime(1).eq.0 ) then
          print*
          print*,' ERROR: reference date is missing / cannot be defined'
          stop
      endif

c     Set value and sign of time range
      reftime(6) = timerange
         
c     Write some status information      
      print*,'---- REFERENCE DATE---------- ---------------------------'
      print*
      print*,' Reference time (year)  :',reftime(1)
      print*,'                (month) :',reftime(2)
      print*,'                (day)   :',reftime(3)
      print*,'                (hour)  :',reftime(4)
      print*,'                (min)   :',reftime(5)
      print*,' Time range             :',reftime(6),' min'
      print*
      print*,
     >       '---- STARTING POSITIONS + TRAJECTORIES ------------------'
      print*
      print*,' ntra                    : ',ntra
      print*,' ntim                    : ',ntim
      print*,' ncol                    : ',ncol
      print*,' ntrace                  : ',ntrace
      print*
     
c     -----------------------------------------------------------------------
c     Set the list of needed input files
c     -----------------------------------------------------------------------

c     Now get a list of trajectory times that have to be calculated
      timemin = tt0(1) 
      timemax = tt0(1) 
      do i=2,ntra
         if ( tt0(i).lt.timemin ) timemin = tt0(i)
         if ( tt0(i).gt.timemax ) timemax = tt0(i)
      enddo
      
c     Add timerange to time window (times in min)
      if ( fbflag.eq.1 ) then
         timemax = timemax + timerange
      else 
         timemin = timemin + timerange
      endif
      
c     Find nearest multiples of timeinc
      if ( mod(timemax,timeinc).ne.0 ) then
         if ( timemax.gt.0 ) then
            timemax = (timemax/timeinc+1) * timeinc
         else
            timemax = (timemax/timeinc  ) * timeinc
         endif
      endif
      if ( mod(timemin,timeinc).ne.0 ) then
         if ( timemin.gt.0 ) then
            timemin = (timemin/timeinc  ) * timeinc
         else
            timemin = (timemin/timeinc-1) * timeinc
         endif
      endif
  
c     Set the number of files needed
      numdat = ( timemax - timemin ) / timeinc + 1

c     Wriet total time range to screen
      print*,'---- TIME RANGE OF ALL TRAJECTORIES ---------------------'
      print*
      print*,' Min(Time) [min]        :',timemin
      print*,' Max(Time) [min]        :',timemax
      print*

c     Now get the list of dates needed 
      if ( fbflag.eq.1 ) then
         do i=1,numdat
         
           itime(i) = timemin + (i-1) * timeinc
           date1(1) = reftime(1)
           date1(2) = reftime(2)
           date1(3) = reftime(3)
           date1(4) = reftime(4)
           date1(5) = reftime(5)
           call newdate(date1,real(itime(i)),date2)
           call datestring ( dat(i), date2(1),date2(2),date2(3),
     >                               date2(4),date2(5))
         enddo

      elseif (fbflag.eq.-1) then
         do i=1,numdat

           itime(i) = timemin + (numdat-i) * timeinc
           date1(1) = reftime(1)
           date1(2) = reftime(2)
           date1(3) = reftime(3)
           date1(4) = reftime(4)
           date1(5) = reftime(5)
           call newdate(date1,real(itime(i)),date2)
           call datestring ( dat(i), date2(1),date2(2),date2(3),
     >                               date2(4),date2(5))
         enddo

      endif
      
c     Check whether the list of files needed is of correct format
      do i=1,numdat
        datestr = dat(i)
        if ( (datformat.eq.'yyyymmdd_hh' ).and.
     >       (datestr(12:13).ne.'00') ) 
     >  then
           print*,' ERROR: date format of input must be yyyymmdd_hh'
           print*,'       ',trim(datestr)
           stop 
        endif
      enddo      

c     Write list of dates
      do i=1,numdat
           if ( (mod(i,10).eq.0).or.(i.eq.1).or.
     >          (i.eq.numdat).and.(mod(i,10).ne.0).or.
     >          (numdat.lt.20) )
     >     then
              write(*,'(i4,i9,f10.2,a20)')
     >         i,itime(i),real(itime(i))/60.,prefix//trim(dat(i))
           endif
      enddo
      print*

c     Now check that all files are available
      do i=1,numdat
         if ( datformat.eq.'yyyymmdd_hhmm' ) then
            datestr = dat(i)
            if ( caltra_dir.eq.'nil' ) then
               filename = ''
            else
               filename = caltra_dir
            endif
            if ( caltra_struct.eq.'yyyy.mm' ) then
               yyyy = datestr(1:4)
               mm   = datestr(5:6)
               filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
            endif         
            filename = trim(filename)//'/'//prefix//trim(datestr)
            INQUIRE(FILE=filename, EXIST=file_exists)
            if ( file_exists.eqv..true. ) goto 40
         endif
         if ( datformat.eq.'yyyymmdd_hh' ) then
            datestr = dat(i)
            datestr = datestr(1:11)
            if ( caltra_dir.eq.'nil' ) then
               filename = ''
            else
               filename = caltra_dir
            endif
            if ( caltra_struct.eq.'yyyy.mm' ) then
               yyyy = datestr(1:4)
               mm   = datestr(5:6)
               filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
            endif    
            filename = trim(filename)//prefix//trim(datestr)
            INQUIRE(FILE=filename, EXIST=file_exists)
            if ( file_exists.eqv..true. ) goto 40
         endif
         print*,' ERROR: file is missing : ',trim(filename)
         stop
  40     continue       
      enddo

c     Write some status information
      print* 
      print*,'  Forward/Backward       : ',fbflag
      print*,'  #input files           : ',numdat
      print*,'  First/last input file  : ',trim(dat(1)),' ... ',
     >                                     trim(dat(numdat))
      print*,'  time increment         : ',timeinc
      print*,'  Output file            : ',trim(cdfname)
      print*,'  Time step (min)        : ',ts
      print*,'  Output time interval   : ',deltout
      print*,'  Input pressure AGL     : ',aglflag
      print*

c     Check whether starting times are compatible with input data
      do i=1,ntra
         delta = mod( tt0(i) - timemin, ts )
         if ( delta.ne.0 ) then
            print*,' ERROR: starting time',tt0(i),
     >             'not compatible with input files:',ts
            stop
         endif
      enddo

c     -----------------------------------------------------------------------
c     Get the constant grid parameters (allocate memory for met arrays)
c     -----------------------------------------------------------------------

c     Read the constant grid parameters (nx,ny,nz,xmin,xmax,ymin,ymax,pollon,pollat)
c     The negative <-fid> of the file identifier is used as a flag for parameter retrieval  
      if ( datformat.eq.'yyyymmdd_hh' ) then
        datestr = dat(1)
        datestr = datestr(1:11)
      elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
        datestr = dat(1)
      endif
      if ( caltra_dir.eq.'nil' ) then
         filename = ''
      else
         filename = caltra_dir
      endif
      if ( caltra_struct.eq.'yyyy.mm' ) then
         yyyy = datestr(1:4)
         mm   = datestr(5:6)
         filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
      endif  
      filename = trim(filename)//prefix//trim(datestr)
      varname  = 'U'
      nx       = 1
      ny       = 1
      nz       = 1
      tload   = 0.
      call input_open (fid,filename)
      call input_grid (-fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >                 tload,pollon,pollat,rd,rd,nz,rd,rd,rd,timecheck)
      call input_close(fid)

C     Check if the number of levels is too large
      if (nz.gt.nlevmax) goto 993

C     Set logical flag for periodic data set (hemispheric or not)
      hem = 0
      delta=xmax-xmin-360.
      if (abs(delta+dx).lt.eps) then               ! Program aborts: arrays must be closed
         goto 992
      else if (abs(delta).lt.eps) then              ! Periodic and hemispheric
         hem=1
         per=360.
      endif

C     Allocate memory for some meteorological arrays
      allocate(spt0(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt0 ***'   ! Surface pressure
      allocate(spt1(nx*ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array spt1 ***'
      allocate(uut0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array uut0 ***'   ! Zonal wind
      allocate(uut1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array uut1 ***'
      allocate(vvt0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array vvt0 ***'   ! Meridional wind
      allocate(vvt1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array vvt1 ***'
      allocate(wwt0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array wwt0 ***'   ! Vertical wind
      allocate(wwt1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array wwt1 ***'
      allocate(p3t0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t0 ***'   ! Pressure
      allocate(p3t1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array p3t1 ***'

C     Allocate memory for tracing arrays
      allocate(fld0(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array fld0 ***'   ! Tracing field
      allocate(fld1(nx*ny*nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array fld1 ***'
      
c     Write some status information
      print*
      print*,'---- CONSTANT GRID PARAMETERS ---------------------------'
      print*
      print*,'  xmin,xmax     : ',xmin,xmax
      print*,'  ymin,ymax     : ',ymin,ymax
      print*,'  dx,dy         : ',dx,dy
      print*,'  pollon,pollat : ',pollon,pollat
      print*,'  nx,ny,nz      : ',nx,ny,nz
      print*,'  per, hem      : ',per,hem
      print*
      
c     -----------------------------------------------------------------------
c     Loop to calculate trajectories
c     -----------------------------------------------------------------------

c     Write some status information
      print*,'---- TRAJECTORIES ----------- ---------------------------'
      print*   

C     Save starting positions, and set all other times to missing data
c      itim = 1
c      do i=1,ntra
c         traout(i,itim,1) = real(tt0(i))
c         traout(i,itim,2) = xx0(i)
c         traout(i,itim,3) = yy0(i)
c         traout(i,itim,4) = pp0(i)    
c      enddo
      
c     Set all other fields to mdv (sufficeint to mask p column)
      do i=1,ntra
        do j=2,ntim
          do k=1,ncol
             traout(i,j,k) = -999.
          enddo
        enddo
      enddo
      
c     Init the flag and the counter for trajectories leaving the domain
      leftcount=0
      do i=1,ntra
         leftflag(i)=0
      enddo
      
c     If input data is AGL then determine absolute pressure
c     In any case check if below topography
    
      do i=2,numdat
    
C       read first file
        if ( datformat.eq.'yyyymmdd_hh' ) then
            datestr = dat(i-1)
            datestr = datestr(1:11)
        elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
            datestr = dat(i-1)
        endif
        if ( caltra_dir.eq.'nil' ) then
           filename = ''
        else
           filename = caltra_dir
        endif
        if ( caltra_struct.eq.'yyyy.mm' ) then
           yyyy = datestr(1:4)
           mm   = datestr(5:6)
           filename = trim(filename)//'/'//trim(yyyy)
     >                //'/'//trim(mm)//'/'
        endif  
        filename = trim(filename)//trim(prefix)//trim(datestr)
                         call input_open (fid,filename)
        varname='P'                               
        call input_grid
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t0,spt0,nz,ak,bk,stagz,timecheck)
        call input_close(fid)
        
C       read second file
        if ( datformat.eq.'yyyymmdd_hh' ) then
            datestr = dat(i)
            datestr = datestr(1:11)
        elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
            datestr = dat(i)
        endif
        if ( caltra_dir.eq.'nil' ) then
          filename = ''
        else
           filename = caltra_dir
        endif
        if ( caltra_struct.eq.'yyyy.mm' ) then
           yyyy = datestr(1:4)
           mm   = datestr(5:6)
           filename = trim(filename)//'/'//trim(yyyy)
     >                //'/'//trim(mm)//'/'
        endif  
        filename = trim(filename)//trim(prefix)//trim(datestr)
                         call input_open (fid,filename)
        varname='P'                               
        call input_grid
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >       tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
        call input_close(fid)
      
         do j=1,ntra
        
            if ((fbflag .eq. 1 .and. (tt0(j) .ge. itime(i-1) .and. 
     >            tt0(j)  .lt. itime(i))) .or.
     >           (fbflag .eq. -1 .and. (tt0(j) .le. itime(i-1) .and. 
     >            tt0(j) .gt. itime(i)))) then
    
                reltpos0 = fbflag*real(tt0(j)-itime(i-1))/real(timeinc)
                
C                   Interpolate surface pressure to actual position (from first input file)
                    x1 = xx0(j)
                    y1 = yy0(j)
                    call get_index4 (xind,yind,pind,x1,y1,1050.,
     >                       reltpos0,
     >                       p3t0,p3t1,spt0,spt1,3,
     >                       nx,ny,nz,xmin,ymin,dx,dy,mdv)
                    sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,
     >                                   1.,0.,mdv) 

                if (aglflag .eq. 1) then

C                   get actual pressure
                    pp0(j) = sp - pp0(j)
                
                end if
                
                if ( pp0(j).gt.sp ) then
                    write(*,'(a30,4f10.2)')
     >                  'WARNING: starting point below topography ',
     >                   xx0(j),yy0(j),pp0(j),sp
                    leftflag(j) = 1
                end if

            end if
        
        end do
      
      end do
      
c     Get 3D and surface pressure from first data file 
      if ( datformat.eq.'yyyymmdd_hh' ) then
        datestr = dat(1)
        datestr = datestr(1:11)
      elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
        datestr = dat(1)
      endif
      if ( caltra_dir.eq.'nil' ) then
         filename = ''
      else
         filename = caltra_dir
      endif
      if ( caltra_struct.eq.'yyyy.mm' ) then
          yyyy = datestr(1:4)
          mm   = datestr(5:6)
          filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
      endif  
      filename = trim(filename)//prefix//trim(datestr)
      call input_open (fid,filename)
      varname = 'P'
      call input_grid
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)
      call input_close(fid)
     
c     Read wind fields and vertical grid from first file
      if ( datformat.eq.'yyyymmdd_hh' ) then
        datestr = dat(1)
        datestr = datestr(1:11)
      elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
        datestr = dat(1)
      endif
      if ( caltra_dir.eq.'nil' ) then
         filename = ''
      else
         filename = caltra_dir
      endif
      if ( caltra_struct.eq.'yyyy.mm' ) then
         yyyy = datestr(1:4)
         mm   = datestr(5:6)
         filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
      endif  
      filename = trim(filename)//prefix//trim(datestr)
      tload    = real(itime(1))
      write(*,'(a16,a20,f10.2)') '  (file,time) : ',
     >                      prefix//trim(datestr),tload/60.

      call input_open (fid,filename)

      varname='U'                                      ! U
      call input_wind 
     >    (fid,varname,uut1,tload,stagz,mdv,
     >     xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

      varname='V'                                      ! V
      call input_wind 
     >    (fid,varname,vvt1,tload,stagz,mdv,
     >     xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

      varname='OMEGA'                                  ! OMEGA
      call input_wind
     >       (fid,varname,wwt1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

      varname='P'                                      ! GRID - P,PS
      call input_grid                                  !
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)

      call input_close(fid)

c     Mark all trajectories as inactive
      do i=1,ntra
        active(i) = 0
      enddo

c     Loop over all input files (time step is <timeinc>)
      do itm=1,numdat-1

c       Calculate actual and next time
        time0 = itime(itm  ) 
        time1 = itime(itm+1) 
        
c       Copy old velocities to new ones
        do i=1,nx*ny*nz
           uut0(i)=uut1(i)
           vvt0(i)=vvt1(i)
           wwt0(i)=wwt1(i)
           p3t0(i)=p3t1(i)
        enddo
        do i=1,nx*ny
           spt0(i)=spt1(i)
        enddo

c       Read wind fields and surface pressure at next time
        if ( datformat.eq.'yyyymmdd_hh' ) then
           datestr = dat(itm+1)
           datestr = datestr(1:11)
        elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
           datestr = dat(itm+1)
        endif
        if ( caltra_dir.eq.'nil' ) then
          filename = ''
        else
          filename = caltra_dir
        endif
        if ( caltra_struct.eq.'yyyy.mm' ) then
            yyyy = datestr(1:4)
            mm   = datestr(5:6)
            filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
        endif  
        filename = trim(filename)//prefix//trim(datestr)
        tload    = real( itime(itm+1) )

        call frac2hhmm(real(time1),tload)
        write(*,'(a16,a20,f10.2)') '  (file,time) : ',
     >                          prefix//trim(datestr),tload/60.
     
        call input_open (fid,filename)

        varname='U'                                     ! U
        call input_wind 
     >       (fid,varname,uut1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

        varname='V'                                     ! V
        call input_wind 
     >       (fid,varname,vvt1,tload,stagz,mdv,
     >        xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

        varname='OMEGA'                              ! OMEGA
        call input_wind
     >          (fid,varname,wwt1,tload,stagz,mdv,
     >           xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

        varname='P'                                      ! P,PS
        call input_grid                                  !
     >       (fid,varname,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        tload,pollon,pollat,p3t1,spt1,nz,ak,bk,stagz,timecheck)

        call input_close(fid)

c       Loop over all trajectories - 
        nactive   = 0
        file_load = 0
        do iloop=0,(timeinc/ts-1)+1

C         Calculate relative time position in the interval timeinc (0=beginning, 1=end)
          reltpos0 = real( iloop    *ts ) / real(timeinc)
          reltpos1 = real( (iloop+1)*ts ) / real(timeinc)

c         Set the time of the trajectory step
          time = time0 + iloop*ts*fbflag

C         Initialize counter for domain leaving trajectories
          leftcount=0

c         Timestep for all trajectories
          do i=1,ntra

c           Set the flag whether the trajectory is active or not
            active(i) = 1
            if ( fbflag.eq.1 ) then
               if ( time.lt.tt0(i)           ) active(i)=0
               if ( time.gt.tt0(i)+timerange ) active(i)=0
            else
               if ( time.lt.tt0(i)+timerange ) active(i)=0
               if ( time.gt.tt0(i)           ) active(i)=0
            endif
            if ( active(i).eq.1 ) nactive = nactive + 1

          enddo
          
c         Decide whether trajectory has to be saved (active(i) = 2)
          do i=1,ntra
            delta = mod(time-tt0(i),deltout)
            if ( (delta.eq.0).and.(active(i).eq.1) ) then
                active(i)  = 2
                file_load  = 1
            endif
          enddo

C         Save positions only every deltout minutes
          do i=1,ntra
            if ( active(i).eq.2 ) then
              indt = fbflag * (time - tt0(i)) / deltout + 1
              call frac2hhmm(real(time),tload)
              traout(i,indt,1) = tload
              if ( abs( pp0(i) + 999. ).lt.eps ) then
                   traout(i,indt,2) = -999.
                   traout(i,indt,3) = -999.
                   traout(i,indt,4) = -999.
              else
                   traout(i,indt,2) = xx0(i)
                   traout(i,indt,3) = yy0(i)
                   traout(i,indt,4) = pp0(i)
              endif
           endif
          enddo

c         Avoid double-updating at time instances with input files
          if ( iloop.eq.(timeinc/ts) ) goto 100

c         Update the trajectory position
          do i=1,ntra

c           If the trajectory is not active, nothing to do 
            if ( active(i).eq.0 ) goto 100

C           Check if trajectory has already left the data domain
            if (leftflag(i).ne.1) then	

c             Iterative Euler timestep (x0,y0,p0 -> x1,y1,p1)
              call euler_3d(
     >               xx1,yy1,pp1,leftflag(i),
     >               xx0(i),yy0(i),pp0(i),reltpos0,reltpos1,
     >               real(ts)*60.,numit,jflag,mdv,wfactor,fbflag,
     >               spt0,spt1,p3t0,p3t1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >               xmin,ymin,dx,dy,per,hem,nx,ny,nz)
         
c             Update trajectory position, or increase number of trajectories leaving domain
              if (leftflag(i).eq.1) then
                leftcount=leftcount+1
                print*,'     -> Trajectory ',i,' leaves domain'
                xx0(i)=-999.
                yy0(i)=-999.
                pp0(i)=-999.
              else
                xx0(i)=xx1
                yy0(i)=yy1
                pp0(i)=pp1
              endif

c          Trajectory has already left data domain (mark as <mdv>)
           else
              xx0(i)=-999.
              yy0(i)=-999.
              pp0(i)=-999.
           endif
           
c          Exit point for inactive trajectories
 100       continue

          enddo

c       End loop over times steps between two input files
        enddo

c       Count the number of active trajectories
        print*,'     -------------------------------------------------',
     >         nint(100. * real(nactive)/real((timeinc/ts)*ntra)),' %'

c       Tracing - but only if any of the positions was updated
        if ( file_load.eq.1 ) then

c         Loop over all tracing variables
          do j=1,ntrace

             print*,'   -> Now tracing ',trim(tvar(j)),' ',trim(tfil(j))

c            Decide whether tracing from constant file
             INQUIRE(FILE=tfil(j), EXIST=file_exists)
             if ( file_exists.eqv..true. ) then
                  varname  = tvar(j)
                  filename = tfil(j)
                  call input_open (fid,filename)
                  call input_wind
     >                 (fid,varname,fld0,tload,stagz,mdv,
     >                  xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
                  call input_close(fid)
                  print*,'   -> R ',trim(tvar(j)),' from ', 
     >                                     trim(tfil(j))
                  do k=1,nx*ny*nz
                    fld1(k) = fld0(k)
                  enddo   
                  goto 150
             endif     

c            Load field at first time
             if ( datformat.eq.'yyyymmdd_hh' ) then
                datestr = dat(itm)
                datestr = datestr(1:11)
             elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
                datestr = dat(itm)
             endif
             if ( caltra_dir.eq.'nil' ) then
               filename = ''
             else
               filename = caltra_dir
             endif
             if ( caltra_struct.eq.'yyyy.mm' ) then
               yyyy = datestr(1:4)
               mm   = datestr(5:6)
               filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
             endif  
             filename = trim(filename)//trim(tfil(j))//trim(datestr)
             varname  = tvar(j)
             tload    = real( itime(itm) )
             call input_open (fid,filename)
             call input_wind
     >                 (fid,varname,fld0,tload,stagz,mdv,
     >                  xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
             call input_close(fid)
             print*,'   -> R ',trim(tvar(j)),' from ', 
     >                                     trim(tfil(j))//trim(datestr)

c            Load field at second time
             if ( datformat.eq.'yyyymmdd_hh' ) then
                datestr = dat(itm+1)
                datestr = datestr(1:11)
             elseif ( datformat.eq.'yyyymmdd_hhmm' ) then
                datestr = dat(itm+1)
             endif
             if ( caltra_dir.eq.'nil' ) then
               filename = ''
             else
               filename = caltra_dir
             endif
             if ( caltra_struct.eq.'yyyy.mm' ) then
               yyyy = datestr(1:4)
               mm   = datestr(5:6)
               filename = trim(filename)//'/'//trim(yyyy)
     >                                  //'/'//trim(mm)//'/'
             endif  
             filename = trim(filename)//trim(tfil(j))//trim(datestr)
             varname  = tvar(j)
             tload    = real( itime(itm+1) )
             call input_open (fid,filename)
             call input_wind
     >                 (fid,varname,fld1,tload,stagz,mdv,
     >                  xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)
             call input_close(fid)
             print*,'   -> R ',trim(tvar(j)),' from ',
     >                                   trim(tfil(j))//trim(datestr)

c            Loop over trajectories and all times
  150        continue           
             do i=1,ntra
             
c                Set which time indices are covered by dat(itm)...dat(itm+1)
                 indt0 = fbflag* (time0 - tt0(i)) / deltout + 1          
                 indt1 = fbflag* (time1 - tt0(i)) / deltout + 1  
             
                 do indt=indt0,indt1
                 
c                    Skip if outside trajectory range
                     if ( (indt.lt.1).or.(indt.gt.ntim) ) goto 200 
                 
c                    Set the time of the the trajectory step
                     time = nint( traout(i,indt,1) ) 

C                    Calculate relative time position in the interval timeinc (0=beginning, 1=end)
                     reltpos0 = fbflag * 
     >                    real(time - time0) / real(timeinc)
                     if ( (reltpos0.lt.-eps).or.
     >                    (reltpos0.gt.(1.+eps) ) ) goto 200

c                    Extract position of trajectory                    
                     x0   = traout(i,indt,2)
                     y0   = traout(i,indt,3)
                     p0   = traout(i,indt,4)

c                    Do the interpolation
                     if ( (abs(x0-mdv).gt.eps).and.
     >                    (abs(x0-mdv).gt.eps) )
     >               then
                         call get_index4 (xind,yind,pind,x0,y0,p0,
     >                                   reltpos0,p3t0,p3t1,spt0,spt1,3,
     >                                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
                     else
                          xind = mdv
                          yind = mdv
                          pind = mdv
                     endif

                     if ( pind.lt.1. ) pind = 1.

c                    Bilinear interpolation
                     if ( nearestflag.eq.0 ) then
                                         
                       if ( (xind.ge.1.).and.(xind.le.real(nx)).and.
     >                      (yind.ge.1.).and.(yind.le.real(ny)).and.
     >                      (pind.ge.1.).and.(pind.le.real(nz)) )
     >                 then
                            f0 = int_index4(fld0,fld1,nx,ny,nz,
     >                                      xind,yind,pind,reltpos0,mdv)
                       else
                            f0 = mdv
                       endif
                     
c                    Nearest-neighbor interpolation
                     else
                       
                       xind = real ( nint(xind) )
                       yind = real ( nint(yind) )
                       pind = real ( nint(pind) )
                       
                       if ( xind.lt.       1 ) xind = 1.
                       if ( xind.gt.real(nx) ) xind = real(nx)
                       if ( yind.lt.       1 ) yind = 1.
                       if ( yind.gt.real(ny) ) yind = real(ny)
                       if ( pind.lt.       1 ) pind = 1.
                       if ( pind.gt.real(nz) ) pind = real(nz)

                       if ( reltpos0.lt.0.5 ) then
                            f0 = int_index4(fld0,fld0,nx,ny,nz,
     >                                      xind,yind,pind,0.,mdv)
                       else
                            f0 = int_index4(fld1,fld1,nx,ny,nz,
     >                                      xind,yind,pind,0.,mdv)
                       endif
                                           
                     endif
                     
c                    Save the interpolated  field                     
                     traout(i,indt,4+j) = f0

c                    Exit point for time loop
200                  continue 
                  
c                 End loop over times                    
                  enddo
                  
c              End loop over trajectories                  
               enddo
               
c            End loop over tracing variables               
             enddo

          endif
      
c     End loop over input files
      enddo
     
c     -----------------------------------------------------------------------
c     Write output
c     -----------------------------------------------------------------------

c     Change time format to hh.mm
      do i=1,ntra
         do j=1,ntim
         
            tload = traout(i,j,1)/60.
            call frac2hhmm(tload,traout(i,j,1))
         
         enddo
      enddo
      
c     Apply scaling to tracing fields
      if ( ntrace.gt.0 ) then
         do i=1,ntra  
            do j=1,ntim
               do k=1,ntrace
                  traout(i,j,4+k) = tfac(k) * traout(i,j,4+k)
               enddo
            enddo
         enddo
      endif

c     Write trajectory file
      vars(1)  ='time'
      vars(2)  ='lon'
      vars(3)  ='lat'
      vars(4)  ='p'
      do i=1,ntrace
         vars(4+i) = tvar(i)
      enddo
      call wopen_tra
     >      (cdfid,cdfname,ntra,ntim,4+ntrace,reftime,vars,outmode)
      call write_tra(cdfid,traout,ntra,ntim,4+ntrace,outmode)
      call close_tra(cdfid,outmode)   

c     Write some status information, and end of program message
      print*  
      print*,'---- STATUS INFORMATION --------------------------------'
      print*
      print*,'  #leaving domain    ', leftcount
      print*,'  #staying in domain ', ntra-leftcount
      print*
      print*,'              *** END OF PROGRAM CALTRA ***'
      print*,'========================================================='

      stop

c     ------------------------------------------------------------------
c     Exception handling
c     ------------------------------------------------------------------

 991  write(*,*) '*** ERROR: all start points outside the data domain'
      call exit(1)
      
 992  write(*,*) '*** ERROR: close arrays on files (prog. closear)'
      call exit(1)

 993  write(*,*) '*** ERROR: problems with array size'
      call exit(1)

      end 


c     *******************************************************************
c     * Time step 
c     *******************************************************************

C     Time-step from (x0,y0,p0) to (x1,y1,p1)
C
C     (x0,y0,p0) input	coordinates (long,lat,p) for starting point
C     (x1,y1,p1) output	coordinates (long,lat,p) for end point
C     deltat	 input	timestep in seconds
C     numit	 input	number of iterations
C     jump	 input  flag (=1 trajectories don't enter the ground)
C     left	 output	flag (=1 if trajectory leaves data domain)

c     -------------------------------------------------------------------
c     Iterative Euler time step (KINEMATIC 3D TRAJECTORIES)
c     -------------------------------------------------------------------

      subroutine euler_3d(x1,y1,p1,left,x0,y0,p0,reltpos0,reltpos1,
     >                deltat,numit,jump,mdv,wfactor,fbflag,
     >                spt0,spt1,p3d0,p3d1,uut0,uut1,vvt0,vvt1,wwt0,wwt1,
     >                xmin,ymin,dx,dy,per,hem,nx,ny,nz)

      implicit none

c     Declaration of subroutine parameters
      integer      nx,ny,nz
      real         x1,y1,p1
      integer      left
      real         x0,y0,p0
      real         reltpos0,reltpos1
      real         deltat
      integer      numit
      integer      jump
      real         wfactor
      integer      fbflag
      real         spt0(nx*ny)   ,spt1(nx*ny)
      real         uut0(nx*ny*nz),uut1(nx*ny*nz)
      real         vvt0(nx*ny*nz),vvt1(nx*ny*nz)
      real         wwt0(nx*ny*nz),wwt1(nx*ny*nz)
      real         p3d0(nx*ny*nz),p3d1(nx*ny*nz)
      real         xmin,ymin,dx,dy
      real         per
      integer      hem
      real         mdv

c     Numerical and physical constants
      real         deltay
      parameter    (deltay=1.112E5)  ! Distance in m between 2 lat circles
      real         pi                       
      parameter    (pi=3.1415927)    ! Pi

c     Auxiliary variables
      real         xmax,ymax
      real         xind,yind,pind
      real         u0,v0,w0,u1,v1,w1,u,v,w,sp
      integer      icount
      character    ch

c     Externals    
      real         int_index4
      external     int_index4
    
c     Reset the flag for domain-leaving
      left=0

c     Set the esat-north bounray of the domain
      xmax = xmin+real(nx-1)*dx
      ymax = ymin+real(ny-1)*dy

C     Interpolate wind fields to starting position (x0,y0,p0)
      call get_index4 (xind,yind,pind,x0,y0,p0,reltpos0,
     >                 p3d0,p3d1,spt0,spt1,3,
     >                 nx,ny,nz,xmin,ymin,dx,dy,mdv)
      u0 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      v0 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)
      w0 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos0,mdv)

c     Force the near-surface wind to zero
      if (pind.lt.1.) w0=w0*pind

C     For first iteration take ending position equal to starting position
      x1=x0
      y1=y0
      p1=p0

C     Iterative calculation of new position
      do icount=1,numit

C        Calculate new winds for advection
         call get_index4 (xind,yind,pind,x1,y1,p1,reltpos1,
     >                    p3d0,p3d1,spt0,spt1,3,
     >                    nx,ny,nz,xmin,ymin,dx,dy,mdv)
         u1 = int_index4(uut0,uut1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         v1 = int_index4(vvt0,vvt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)
         w1 = int_index4(wwt0,wwt1,nx,ny,nz,xind,yind,pind,reltpos1,mdv)

c        Force the near-surface wind to zero
         if (pind.lt.1.) w1=w1*pind
 
c        Get the new velocity in between
         u=(u0+u1)/2.
         v=(v0+v1)/2.
         w=(w0+w1)/2.

C        Calculate new positions
         x1 = x0 + fbflag*u*deltat/(deltay*cos(y0*pi/180.))
         y1 = y0 + fbflag*v*deltat/deltay
         p1 = p0 + fbflag*wfactor*w*deltat/100.
         
c       Handle pole problems (crossing and near pole trajectory)
        if ((hem.eq.1).and.(y1.gt.90.)) then
          y1=180.-y1
          x1=x1+per/2.
        endif
        if ((hem.eq.1).and.(y1.lt.-90.)) then
          y1=-180.-y1
          x1=x1+per/2.
        endif
        if (y1.gt.89.99) then
           y1=89.99
        endif

c       Handle crossings of the dateline
        if ((hem.eq.1).and.(x1.gt.xmin+per-dx)) then
           x1=xmin+amod(x1-xmin,per)
        endif
        if ((hem.eq.1).and.(x1.lt.xmin)) then
           x1=xmin+per+amod(x1-xmin,per)
        endif

C       Interpolate surface pressure to actual position
        call get_index4 (xind,yind,pind,x1,y1,1050.,reltpos1,
     >                   p3d0,p3d1,spt0,spt1,3,
     >                   nx,ny,nz,xmin,ymin,dx,dy,mdv)
        sp = int_index4 (spt0,spt1,nx,ny,1,xind,yind,1.,reltpos1,mdv)       
 
c       Handle trajectories which cross the lower boundary (jump flag)
        if ((jump.eq.1).and.(p1.gt.sp)) p1=sp-10.
 
C       Check if trajectory leaves data domain
        if ( ( (hem.eq.0).and.(x1.lt.xmin)    ).or.
     >       ( (hem.eq.0).and.(x1.gt.xmax-dx) ).or.
     >         (y1.lt.ymin).or.(y1.gt.ymax).or.(p1.gt.sp) )
     >  then
          left=1
          goto 100
        endif

      enddo

c     Exit point for subroutine
 100  continue

      return

      end

