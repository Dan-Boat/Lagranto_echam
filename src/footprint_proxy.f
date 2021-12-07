   
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
      
c     Transform into geo system
      do i=1,ntra
         rlon = proxy(i,2)
         rlat = proxy(i,3)
         call getenvir_f (clon,clat,0.,lon,lat,rlon,rlat,1)
         proxy(i,2) = lon
         proxy(i,3) = lat
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
      
