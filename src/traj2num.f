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
      real,allocatable, dimension (:)     :: num         ! Output number

c     Auxiliary variables
      integer                                inpmode
      integer                                stat
      integer                                fid
      integer                                i,j,k
      real                                   dist
      real                                   lon0,lat0,lon1,lat1

c     Externals
      real                                   sdis
      external                               sdis

c     ----------------------------------------------------------------------
c     Preparations
c     ----------------------------------------------------------------------

c     Read parameters
      open(10,file='traj2num.param')
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
      if ( mode.eq.'boost' ) goto 10
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
      allocate(num(ntra),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array num      ***'
      if ( mode.eq. 'proxy' ) then
           allocate(proxy(ntra,ncol+1),stat=stat)
           if (stat.ne.0) print*,'*** error allocating array proxy  ***'  
      endif

c     Read inpufile
      call ropen_tra(fid,inpfile,ntra,ntim,ncol,refdate,vars,inpmode)
      call read_tra (fid,tra,ntra,ntim,ncol,inpmode)
      call close_tra(fid,inpmode)
    
c     ----------------------------------------------------------------------
c     Mode = 'boost': get the maximum distance traveled in one time step
c     ----------------------------------------------------------------------

      if ( mode.ne.'boost') goto 100

      do i=1,ntra

        num(i) = 0.

        do j=2,ntim

c          Get spherical distance between data points
           lon0 = tra(i,j-1,2)
           lat0 = tra(i,j-1,3)
           lon1 = tra(i,j  ,2)
           lat1 = tra(i,j  ,3)
           dist = sdis( lon1,lat1,lon0,lat0 )

           if ( dist.gt.num(i) ) num(i) = dist

        enddo

      enddo

 100  continue
 
c     ----------------------------------------------------------------------
c     Mode = 'proxy': get nearest distance to a lat/lon point
c     ----------------------------------------------------------------------

      if ( mode.ne.'proxy') goto 200

          call get_proxy (proxy,clon,clat,tra,ntra,ntim,ncol,radius)
          do i=1,ntra
            num(i) = proxy(i,ncol+1)
          enddo
   
 200  continue
 
c     ----------------------------------------------------------------------
c     Write output file
c     ----------------------------------------------------------------------

      open(10,file=outfile)
        do i=1,ntra
            write(10,*) num(i)
        enddo
      close(10)
      
      end

c     ***********************************************************************
c     * SUBROUTINES                                                         *
c     ***********************************************************************

c     -----------------------------------------------------------------------
c     Spherical distance between lat/lon points
c     -----------------------------------------------------------------------

      real function sdis(xp,yp,xq,yq)
c
c     calculates spherical distance (in km) between two points given
c     by their spherical coordinates (xp,yp) and (xq,yq), respectively.
c
      real      re
      parameter (re=6370.)
      real      pi180
      parameter (pi180=3.14159/180.)
      real      xp,yp,xq,yq,arg

      arg=sin(pi180*yp)*sin(pi180*yq)+
     >    cos(pi180*yp)*cos(pi180*yq)*cos(pi180*(xp-xq))
      if (arg.lt.-1.) arg=-1.
      if (arg.gt.1.) arg=1.

      sdis=re*acos(arg)

      end
      
c     -----------------------------------------------------------------------
c     Nearest distance to a lat/lon point
c     -----------------------------------------------------------------------

      subroutine get_proxy (proxy,clon,clat,tra,ntra,ntim,ncol,radius)
     
c     Given a set of trajectories tra(ntra,ntim,ncol), find the closest point
c     to clon/clat and return the footprint of this trajectory at this
c     point: proxy(ntra,ncol+1); the parameter <ncol+1> of the trajectory will
c     become the minimum distance.

      implicit none
                
c     Declaration of parameters
      integer    ntra,ntim,ncol
      real       tra(ntra,ntim,ncol)
      real       proxy(ntra,ncol+1)
      real       radius
      real       clon,clat
      
c     Flag for tests
      integer      test
      parameter    (test=0)
      character*80 testfile
      parameter    (testfile='TEST.nc')
      
c     Number of grid points for the radius mask
      integer      nmax 
      parameter    (nmax = 400)     
      real         degkm
      parameter    (degkm = 111.)
      
c     Number of binary search steps
      integer      nbin
      parameter    (nbin=10)

c     Auxiliry variables
      integer      i,j,k
      real         lon,lat,rlon,rlat
      real         dist
      character*80 varname
      integer      rnx,rny
      real         rmask (nmax,nmax)
      real         rcount(nmax,nmax)
      real         rxmin,rymin,rdx,rdy
      integer      indx,indy
      integer      isnear
      real         near,near0,near1,nearm
      integer      j0,j1
      real         r0,r1,rm,frac
      real         rlon0,rlon1,rlonm,rlat0,rlat1,rlatm
      
c     Externals
      real       sdis
      external   sdis

c     Transform into clon/clat centered system
      do i=1,ntra
        do j=1,ntim
            lon = tra(i,j,2)
            lat = tra(i,j,3)
            call getenvir_f (clon,clat,0.,rlon,rlat,lon,lat,1)
            tra(i,j,2) = rlon
            tra(i,j,3) = rlat
        enddo
      enddo
      
c     Init the radius mask - define the grid resolution
      rdx   = 2.5 * radius / degkm / real(nmax)
      rdy   = 2.5 * radius / degkm / real(nmax)
      rnx   = nmax
      rny   = nmax
      rxmin = -real(rnx/2)*rdx
      rymin = -real(rny/2)*rdy
      do i=1,rnx
         do j=1,rny
              rlon = rxmin + real(i-1) * rdx
              rlat = rymin + real(j-1) * rdy
              dist = sdis(rlon,rlat,0.,0.)
              if ( dist.lt.radius ) then
                 rmask(i,j) = dist
              else
                rmask(i,j) = radius
              endif
         enddo
      enddo

c     Init counter for test
      if ( test.eq.1 ) then
        do i=1,rnx
          do j=1,rny
            rcount(i,j) = 0.
          enddo
        enddo
      endif


c     Loop over all trajectories
      do i=1,ntra
      
c        Decide whether trajectory comes close to point   
         isnear = 0
         near   = radius
         do j=1,ntim
            indx = nint( ( tra(i,j,2) - rxmin ) / rdx + 1. )
            indy = nint( ( tra(i,j,3) - rymin ) / rdy + 1. )
            if ( (indx.ge.1).and.(indx.le.rnx).and.
     >           (indy.ge.1).and.(indy.le.rny) )
     >      then
               if (rmask(indx,indy).lt.near) then
                   near   = rmask(indx,indy)
                   isnear = j
               endif
            endif
         enddo

c        No close point was found - go to next trajectory
         do k=1,ncol
            proxy(i,k) = tra(i,1,k)
         enddo
         proxy(i,ncol+1) = radius
         if ( isnear.eq.0 ) goto 310

c        Get the exact position and time with binary search
         j0 = isnear - 1
         if ( j0.eq.0 ) j0 = 1
         rlon0 = tra(i,j0,2)
         rlat0 = tra(i,j0,3)
         near0 = sdis(rlon0,rlat0,0.,0.)
         r0    = real(j0)

         j1 = isnear + 1
         if ( j1.gt.ntim ) j1 = ntim
         rlon1 = tra(i,j1,2)
         rlat1 = tra(i,j1,3)
         near1 = sdis(tra(i,j1,2),tra(i,j1,3),0.,0.)
         r1    = real(j1)

         do k=1,nbin
            rm    = 0.5 * ( r0 + r1 )
            rlatm = 0.5 * ( rlat0 + rlat1 )
            rlonm = 0.5 * ( rlon0 + rlon1 )
            nearm = sdis(rlonm,rlatm,0.,0.)
            if ( near0.lt.near1 ) then
               r1    = rm
               rlon1 = rlonm
               rlat1 = rlatm
               near1 = nearm
            else
               r0    = rm
               rlon0 = rlonm
               rlat0 = rlatm
               near0 = nearm
            endif
         enddo

c        Now get the final distance and position
         proxy(i,ncol+1) = nearm
         if ( test.eq.1 ) then
            indx = nint( ( rlonm - rxmin ) / rdx + 1. )
            indy = nint( ( rlatm - rymin ) / rdy + 1. )
            if ( (indx.ge.1).and.(indx.le.rnx).and.
     >           (indy.ge.1).and.(indy.le.rny) )
     >      then
               rcount(indx,indy) = rcount(indx,indy) + 1.
            endif
         endif

c        Get all values at exactly this position
         j0   = int(rm)
         frac = rm - real(j0)
         j1 = j0 + 1
         if (j1.gt.ntim ) j1 = ntim
         do k=1,ncol
            proxy(i,k) = (1.-frac)*tra(i,j0,k)+frac*tra(i,j1,k)
         enddo

c        Next trajectory
310      continue

      enddo
   
c     Write radius mask to netCDF file for tests
      if ( test.eq.1 ) then
        varname = 'MASK'
        call writecdf2D(testfile,varname,rmask,0.,
     >                    rdx,rdy,rxmin,rymin,rnx,rny,1,1)
        varname = 'COUNT'
        call writecdf2D(testfile,varname,rcount,0.,
     >                    rdx,rdy,rxmin,rymin,rnx,rny,0,1)
      endif

      end



c     --------------------------------------------------------------------------------
c     Forward coordinate transformation (True lon/lat -> Rotated lon/lat)
c     --------------------------------------------------------------------------------

      SUBROUTINE getenvir_f (clon,clat,rotation,
     >                       rlon,rlat,lon,lat,n)

      implicit none

c     Declaration of input and output parameters
      integer     n
      real        clon,clat,rotation
      real        lon(n), lat(n)
      real        rlon(n),rlat(n)

c     Auxiliary variables 
      real         pollon,pollat
      integer      i
      real         rlon1(n),rlat1(n)

c     Externals
      real     lmtolms,phtophs
      external lmtolms,phtophs

c     First rotation
      pollon=clon-180.
      if (pollon.lt.-180.) pollon=pollon+360.
      pollat=90.-clat
      do i=1,n

c        First rotation
         pollon=clon-180.
         if (pollon.lt.-180.) pollon=pollon+360.
         pollat=90.-clat
         rlon1(i)=lmtolms(lat(i),lon(i),pollat,pollon)
         rlat1(i)=phtophs(lat(i),lon(i),pollat,pollon)            

c        Second coordinate transformation
         pollon=-180.
         pollat=90.+rotation
         rlon(i)=90.+lmtolms(rlat1(i),rlon1(i)-90.,pollat,pollon)
         rlat(i)=phtophs(rlat1(i),rlon1(i)-90.,pollat,pollon)   
         
      enddo

      END   

c     --------------------------------------------------------------------------------
c     Backward coordinate transformation (Rotated lon/lat -> True lon/lat)
c     --------------------------------------------------------------------------------

      SUBROUTINE getenvir_b (clon,clat,rotation,
     >                       lon,lat,rlon,rlat,n)

      implicit none

c     Declaration of input and output parameters
      integer     n
      real        clon,clat,rotation
      real        lon(n), lat(n)
      real        rlon(n),rlat(n)

c     Auxiliary variables 
      real         pollon,pollat
      integer      i
      real         rlon1(n),rlat1(n)

c     Externals
      real         lmstolm,phstoph
      external     lmstolm,phstoph

c     First coordinate transformation (make the local coordinate system parallel to equator)
      pollon=-180.
      pollat=90.+rotation
      do i=1,n
         rlon1(i)=90.+lmstolm(rlat(i),rlon(i)-90.,pollat,pollon)
         rlat1(i)=phstoph(rlat(i),rlon(i)-90.,pollat,pollon)            
      enddo

c     Second coordinate transformation (make the local coordinate system parallel to equator)
      pollon=clon-180.
      if (pollon.lt.-180.) pollon=pollon+360.
      pollat=90.-clat
      do i=1,n
         lon(i)=lmstolm(rlat1(i),rlon1(i),pollat,pollon)
         lat(i)=phstoph(rlat1(i),rlon1(i),pollat,pollon)            
      enddo

      END
      
c     --------------------------------------------------------------------------------
c     Transformation routine: LMSTOLM and PHSTOPH from library gm2em
c     --------------------------------------------------------------------------------

      REAL FUNCTION LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C
C**** LMSTOLM  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C**                IM ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHIS   = ZPIR18*PHIS
      ZLAMS   = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS   = ZPIR18*ZLAMS
 
      ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) -
     2          COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) +
     2          SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMSTOLM =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMSTOLAM =  90.0
            ELSE
              LMSTOLAM = -90.0
            ENDIF
      ELSE
        LMSTOLM = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END

      REAL FUNCTION PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C
C**** PHSTOPH  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C****                 ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE BREITE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      SINPOL = SIN(ZPIR18*POLPHI)
      COSPOL = COS(ZPIR18*POLPHI)
      ZPHIS  = ZPIR18*PHIS
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = ZPIR18*ZLAMS
      ARG     = COSPOL*COS(ZPHIS)*COS(ZLAMS) + SINPOL*SIN(ZPHIS)
 
      PHSTOPH = ZRPI18*ASIN(ARG)
 
      RETURN
      END

      REAL FUNCTION LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** LMTOLMS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL =     ZPIR18*POLLAM
      ZPHI    =     ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
 
      ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
      ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMTOLMS =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMTOLMS =  90.0
            ELSE
              LMTOLMS = -90.0
            ENDIF
      ELSE
        LMTOLMS = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END

      REAL FUNCTION PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** PHTOPHS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   ROTIERTE BREITE PHIS ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHI    = ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
      ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
 
      PHTOPHS = ZRPI18*ASIN(ZARG)
 
      RETURN
      END

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



