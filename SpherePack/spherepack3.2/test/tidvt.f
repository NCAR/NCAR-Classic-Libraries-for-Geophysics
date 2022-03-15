c
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     *                                                               *
c     *                  copyright (c) 1998 by UCAR                   *
c     *                                                               *
c     *       University Corporation for Atmospheric Research         *
c     *                                                               *
c     *                      all rights reserved                      *
c     *                                                               *
c     *                      SPHEREPACK version 3.2                   *
c     *                                                               *
c     *       A Package of Fortran77 Subroutines and Programs         *
c     *                                                               *
c     *              for Modeling Geophysical Processes               *
c     *                                                               *
c     *                             by                                *
c     *                                                               *
c     *                  John Adams and Paul Swarztrauber             *
c     *                                                               *
c     *                             of                                *
c     *                                                               *
c     *         the National Center for Atmospheric Research          *
c     *                                                               *
c     *                Boulder, Colorado  (80307)  U.S.A.             *
c     *                                                               *
c     *                   which is sponsored by                       *
c     *                                                               *
c     *              the National Science Foundation                  *
c     *                                                               *
c     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c
c
c     3/97
c
c     a program for testing all divergence,vorticity and idvt(ec,es,gc,gs) routines
c
c     (1) first set a valid vector field by setting a stream function sf and velocity
c         potential function sv as polys in x,y,z restricted to the sphere.  Then
c         derive (v,w) and dv,vt from sf and sv analytically by differentiation.
c         (see tvha.f)
c
c     (2) compute the coefficients br,bi,cr,ci of (v,w) using vector analysis
c
c     (3) compute the divergence and vorticity of (v,w) using div,vrt (es,ec,gc,gs)
c
c     (4) compare with divergence and vorticity from (1)
c
c     (5) invert dv,vt with idvt(ec,es,gc,gs) and compare with vector field from (1)
c
      program tidvt
c
c     set dimensions with parameter statements
c
      parameter(nnlat= 25,nnlon= 16, nnt = 3)
      parameter(mmdab = (nnlon+2)/2, mmdb = (nnlon+1)/2)
      parameter (lleng= 5*nnlat*nnlat*nnlon,llsav=15*nnlat*nnlat*nnlon)
      dimension work(lleng),wsave(llsav)
      parameter (lldwork = 4*nnlat*nnlat )
      double precision dwork(lldwork)
      dimension br(mmdb,nnlat,nnt),bi(mmdb,nnlat,nnt)
      dimension cr(mmdb,nnlat,nnt),ci(mmdb,nnlat,nnt)
      dimension ad(mmdab,nnlat,nnt),bd(mmdab,nnlat,nnt)
      dimension av(mmdab,nnlat,nnt),bv(mmdab,nnlat,nnt)
      dimension dv(nnlat,nnlon,nnt),vt(nnlat,nnlon,nnt)
      dimension dvsav(nnlat,nnlon,nnt),vtsav(nnlat,nnlon,nnt)
      dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
      dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
      dimension vsav(nnlat,nnlon,nnt),wsav(nnlat,nnlon,nnt)
      dimension ptrbd(nnt),ptrbv(nnt)
      double precision dtheta, dwts
c
c     set dimension variables
c
      nlat = nnlat
      nlon = nnlon
      mdab = mmdab
      mdb = mmdb
      nmax = max0(nlat,nlon)

      lwork = lleng
      lsave = llsav
      nt = nnt
      call iout(nlat,"nlat")
      call iout(nlon,"nlon")
      call iout(nt,"  nt")
      isym = 0
      ityp = 0
c
c     set equally spaced colatitude and longitude increments
c
      pi = 4.0*atan(1.0)
      dphi = (pi+pi)/nlon
      dlat = pi/(nlat-1)
c
c     compute nlat gaussian points in thetag
c
      ldwork = lldwork
      call gaqd(nlat,dtheta,dwts,dwork,ldwork,ier)
      do  i=1,nlat
	thetag(i) = dtheta(i)
      end do
      call name("gaqd")
      call iout(ier," ier")
      call vecout(thetag,"thtg",nlat)
c
c     test all divergence and inverse divergence subroutines
c
      do icase=1,4
c
c     icase=1 corresponds to "ec"
c     icase=2 corresponds to "es"
c     icase=3 corresponds to "gc"
c     icase=4 corresponds to "gs"

      call name("****")
      call name("****")
      call iout(icase,"icas")
c
c     set vector field v,w
c
      do k=1,nt
	do j=1,nlon
	  phi = (j-1)*dphi
	  sinp = sin(phi)
	  cosp = cos(phi)
	  do i=1,nlat
	    theta = (i-1)*dlat
	    if (icase.gt.2) theta=thetag(i)
	    cost = cos(theta)
	    sint = sin(theta)
	    x = sint*cosp
	    y = sint*sinp
	    z = cost
	    if (k .eq.1) then
c              sf = x
c              sv = y
c
c          v = -1/sint*dstdp + dsvdt
c
c          w =  1/sint*dsvdp + dstdt
c
c          dv = 1/sint*[d(sint*v)/dt + dwdp]  = dvdt + ct/st*v + 1/st*dwdp
c
c          vt = 1/sint*[-dv/dp + d(sint*w)/dt) = dwdt + ct/st*w - 1/st*dvdp

	       v(i,j,k) = sinp + cost*sinp
	       w(i,j,k) = cosp + cost*cosp
	       dv(i,j,k) = -2.0*sint*sinp
	       vt(i,j,k) = -2.0*sint*cosp
	    else if (k.eq.2) then
C              sf = y
c              sv = z
	       v(i,j,k) = -cosp-sint
	       w(i,j,k) = cost*sinp
	       vt(i,j,k) = -2.*sint*sinp
	       dv(i,j,k) = -2.*cost
	    else if (k.eq.3) then
c              st = x
c              sv = z
	       v(i,j,k) = sinp - sint
	       w(i,j,k) = cost*cosp
	       vt(i,j,k) = -2.*sint*cosp
	       dv(i,j,k) = -2.*cost
	    end if
c
c      save derived vector field,vorticity,divergence for latter comparison
c
	    vtsav(i,j,k) = vt(i,j,k)
	    dvsav(i,j,k) = dv(i,j,k)
	    vsav(i,j,k) = v(i,j,k)
	    wsav(i,j,k) = w(i,j,k)
	  end do
	end do
      end do

      if (icase.eq.1) then

      call name("**ec")
c
c     analyze vector field
c
      call vhaeci(nlat,nlon,wsave,lsave,work,lwork,dwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhaec(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdb,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")
c
c     compute divergence,vorticity of (v,w) in dv,vt
c

      call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divec(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")

      call vrtec(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("vrt ")
      call iout(ierror,"ierr")


      else if (icase.eq.2) then

      call name("**es")
c
c     analyze vector field
c
      call vhaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")
      call vhaes(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdb,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call dives(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")
      call vrtes(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("vrt ")
      call iout(ierror,"ierr")

      else if (icase .eq. 3) then

      call name("**gc")
c
c     analyze vector field
c
      call vhagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")
      call vhagc(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdb,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divgc(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")
      call vrtgc(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("vrt ")
      call iout(ierror,"ierr")

      else if (icase .eq. 4) then

      call name("**gs")
c
c     analyze vector field
c
      call vhagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")
      call vhags(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdb,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divgs(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")
      call vrtgs(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("vrt ")
      call iout(ierror,"ierr")

      end if
c
c     compute "error" in dv,vt
c
      err2d = 0.0
      err2v = 0.0
      do k=1,nt
	do j=1,nlon
	  do i=1,nlat
	      err2d = err2d + (dv(i,j,k)-dvsav(i,j,k))**2
	      err2v = err2v + (vt(i,j,k)-vtsav(i,j,k))**2
	  end do
	end do
      end do
c
c     set and print least squares error in dv
c
      err2d = sqrt(err2d/(nt*nlat*nlon))
      call vout(err2d,"errd")
      err2v = sqrt(err2v/(nt*nlat*nlon))
      call vout(err2v,"errv")
c
c     now recompute (v,w) inverting dv,vt using idvt(ec,es,gc,gs)
c
      do kk=1,nt
      do j=1,nlon
      do i=1,nlat
      v(i,j,kk) = 0.0
      w(i,j,kk) = 0.0
      end do
      end do
      end do


      if (icase.eq.1) then

      call name("**ec")


      call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaec(nlat,nlon,isym,nt,dv,nlat,nlon,ad,bd,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call shaec(nlat,nlon,isym,nt,vt,nlat,nlon,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idvtec(nlat,nlon,isym,nt,v,w,nlat,nlon,ad,bd,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ptrbd,ptrbv,ierror)
      call name("idvt")
      call iout(ierror,"ierr")
      call vecout(prtbd,"prtd",nt)
      call vecout(prtbv,"prtv",nt)

      else if (icase.eq.2) then

      call name("**es")

      call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaes(nlat,nlon,isym,nt,dv,nlat,nlon,ad,bd,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call shaes(nlat,nlon,isym,nt,vt,nlat,nlon,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)

      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idvtes(nlat,nlon,isym,nt,v,w,nlat,nlon,ad,bd,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ptrbd,ptrbv,ierror)
      call name("idvt")
      call iout(ierror,"ierr")
      call vecout(prtbd,"prtd",nt)
      call vecout(prtbv,"prtv",nt)

      else if (icase.eq.3) then

      call name("**gc")

      call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shagc(nlat,nlon,isym,nt,dv,nlat,nlon,ad,bd,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call shagc(nlat,nlon,isym,nt,vt,nlat,nlon,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)

      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idvtgc(nlat,nlon,isym,nt,v,w,nlat,nlon,ad,bd,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ptrbd,ptrbv,ierror)
      call name("idvt")
      call iout(ierror,"ierr")
      call vecout(prtbd,"prtd",nt)
      call vecout(prtbv,"prtv",nt)

      else if (icase.eq.4) then

      call name("**gs")

      call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shags(nlat,nlon,isym,nt,dv,nlat,nlon,ad,bd,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call shags(nlat,nlon,isym,nt,vt,nlat,nlon,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)

      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idvtgs(nlat,nlon,isym,nt,v,w,nlat,nlon,ad,bd,av,bv,
     +mdab,nlat,wsave,lsave,work,lwork,ptrbd,ptrbv,ierror)
      call name("idvt")
      call iout(ierror,"ierr")
      call vecout(prtbd,"prtd",nt)
      call vecout(prtbv,"prtv",nt)

      end if

c
c     compare this v,w with original derived from sf,sv
c
      err2v = 0.0
      err2w = 0.0
      do k=1,nt
	do j=1,nlon
	  do i=1,nlat
	    err2v = err2v + (v(i,j,k)-vsav(i,j,k))**2
	    err2w = err2w + (w(i,j,k)-wsav(i,j,k))**2
	  end do
	end do
      end do
      err2v = sqrt(err2v/(nlat*nlon*nt))
      err2w = sqrt(err2w/(nlat*nlon*nt))
      call vout(err2v,"errv")
      call vout(err2w,"errw")
c
c     end of icase loop
c
      end do
      end
c
      subroutine iout(ivar,nam)
      real nam
      write(6,10) nam , ivar
   10 format(1h a4, 3h = ,i8)
      return
      end
c
      subroutine vout(var,nam)
      real nam
      write(6,10) nam , var
   10 format(1h a4,3h = ,e12.5)
      return
      end
c
      subroutine name(nam)
      real nam
      write(6,100) nam
  100 format(1h a8)
      return
      end
c
      subroutine vecout(vec,nam,len)
      dimension vec(len)
      real nam
      write(6,109) nam, (vec(l),l=1,len)
  109 format(1h a4,/(1h 8e11.4))
      return
      end
