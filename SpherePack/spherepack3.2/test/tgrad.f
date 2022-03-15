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
c     12/96
c
c     a program for testing all gradient and inverse gradient routines
c
c     (1) first a scalar field is set in st by restricting a poly in x,y,z
c         to the sphere surface
c
c     (2) a scalar analysis is used to compute the coefs a,b of st
c
c     (3) a,b are input to the various gradient routines to compute a vector field
c         (v,w)
c
c     (4) the vector field (v,w) is compared with the gradient of st obtained
c         analytically
c
c     (5) the inverse gradient of (v,w) is computed and compared with (1)
c
      program tgrad
c
c     set dimensions with parameter statements
c
      parameter(nnlat= 33,nnlon= 18, nnt = 4)
      parameter  (mmdb=(nnlon+1)/2, mmdab=(nnlon+2)/2)
      parameter (lleng= 5*nnlat*nnlat*nnlon,llsav= 5*nnlat*nnlat*nnlon)
      parameter (lldwork = 4*nnlat*nnlat)
      double precision dwork(lldwork)
      dimension work(lleng),wsave(llsav)
      dimension br(mmdb,nnlat,nnt),bi(mmdb,nnlat,nnt)
      dimension cr(mmdb,nnlat,nnt),ci(mmdb,nnlat,nnt)
      dimension a(mmdab,nnlat,nnt),b(mmdab,nnlat,nnt)
      dimension sf(nnlat,nnlon,nnt)
      dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
      dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
      double precision dtheta, dwts
c
c     set dimension variables
c
      nlat = nnlat
      nlon = nnlon
      lwork = lleng
      lsave = llsav
      mdb = mmdb
      mdab = mmdab
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
c     test all analysis and synthesis subroutines
c
      do icase=1,4
c
c     icase=1 test gradec,igradec
c     icase=2 test grades,igrades
c     icase=3 test gradgc,igradgc
c     icase=4 test gradgs,igradgs
c
      call name("****")
      call name("****")
      call iout(icase,"icas")
c
c
c     set scalar field as polys in x,y,z and set vector field (v,w) using
c     v = dsfdt, w = 1/sint*dsfdp
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
	    dxdt = cost*cosp
	    dxdp = -sint*sinp
	    dydt = cost*sinp
	    dydp = sint*cosp
	    dzdt = -sint
	    dzdp = 0.0
	    if (k.eq.1) then
	       sf(i,j,k) = x*y
	       dsfdt = x*dydt+y*dxdt
	       dsfdp = x*dydp+y*dxdp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = (cosp*dydp+sinp*dxdp)
	    else if (k.eq.2) then
	       sf(i,j,k) = x*z
	       dsfdp = x*dzdp+z*dxdp
	       dsfdt = x*dzdt+z*dxdt
	       v(i,j,k) = dsfdt
	       w(i,j,k) = cosp*dzdp-z*sinp
	    else if (k.eq.3) then
	       sf(i,j,k) = y*z
	       dsfdt = y*dzdt + z*dydt
	       dsfdp = y*dzdp+ z*dydp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = sinp*dzdp + z*cosp
	    else if (k.eq.4) then
	       sf(i,j,k) = x*y*z
	       dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
	       dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
	    end if
	  end do
	end do
      end do

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(sf(1,1,kk),"  sf",nlat,nlon)
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     end do

      if (icase.eq.1) then

      call name("**ec")
c
c     analyze scalar field st
c
      call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaec(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(a(1,1,kk),"   a",nlat,nlat)
c     call aout(b(1,1,kk),"   b",nlat,nlat)
c     end do

c
c     compute gradient of st in v,w
c

      call vhseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhci")
      call iout(ierror,"ierr")

      call gradec(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("grad")
      call iout(ierror,"ierr")

      else if (icase.eq.2) then

      call name("**es")
c
c     analyze scalar field st
c
      call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaes(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(a(1,1,kk),"   a",nlat,nlat)
c     call aout(b(1,1,kk),"   b",nlat,nlat)
c     end do

c
c     compute gradient of st in v,w
c

      call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhsi")
      call iout(ierror,"ierr")

      call grades(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("grad")
      call iout(ierror,"ierr")

      else if (icase .eq. 3) then


      call name("**gc")
      call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shagc(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(a(1,1,kk),"   a",nlat,nlat)
c     call aout(b(1,1,kk),"   b",nlat,nlat)
c     end do

c
c     compute gradient of st in v,w
c

      call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhgc")
      call iout(ierror,"ierr")

      call gradgc(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("grad")
      call iout(ierror,"ierr")


      else if (icase .eq. 4) then

      call name("**gs")

      call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shags(nlat,nlon,isym,nt,sf,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(a(1,1,kk),"   a",nlat,nlat)
c     call aout(b(1,1,kk),"   b",nlat,nlat)
c     end do

c
c     compute gradient of st in v,w
c

      call vhsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhgs")
      call iout(ierror,"ierr")

      call gradgs(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,mdab,nlat,wsave,
     +lsave,work,lwork,ierror)
      call name("grad")
      call iout(ierror,"ierr")
      end if


c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     end do

c
c     compute "error" in v,w
c
      err2v = 0.0
      err2w = 0.0
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
	    dxdt = cost*cosp
	    dxdp = -sint*sinp
	    dydt = cost*sinp
	    dydp = sint*cosp
	    dzdt = -sint
	    dzdp = 0.0
	    if (k.eq.1) then
	       dsfdt = x*dydt+y*dxdt
	       dsfdp = x*dydp+y*dxdp
	       ve = dsfdt
	       we = (cosp*dydp+sinp*dxdp)
	    else if (k.eq.2) then
	       dsfdp = x*dzdp+z*dxdp
	       dsfdt = x*dzdt+z*dxdt
	       ve = dsfdt
	       we = cosp*dzdp-z*sinp
	    else if (k.eq.3) then
	       dsfdt = y*dzdt + z*dydt
	       dsfdp = y*dzdp+ z*dydp
	       ve = dsfdt
	       we = sinp*dzdp + z*cosp
	    else if (k.eq.4) then
	       dsfdt = x*y*dzdt + x*z*dydt + y*z*dxdt
	       dsfdp = x*y*dzdp + x*z*dydp + y*z*dxdp
	       ve = dsfdt
	       we = cosp*y*dzdp+cosp*z*dydp+sinp*z*dxdp
	    end if
	       err2v = err2v + (v(i,j,k)-ve)**2
	       err2w = err2w + (w(i,j,k)-we)**2
	  end do
	end do
      end do
c
c     set and print least squares error in v,w
c
      err2v = sqrt(err2v/(nt*nlat*nlon))
      err2w = sqrt(err2w/(nt*nlat*nlon))
      call vout(err2v,"errv")
      call vout(err2w,"errw")
c
c     now recompute sf by inverting (v,w) using igrad(ec,es,gc,gs)
c
      if (icase.eq.1) then

      call vhaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhaec(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shec")
      call iout(ierror,"ierr")

      call igradec(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("igra")
      call iout(ierror,"ierr")

      else if (icase.eq.2) then
c
c     analyze vector field (v,w)
c
      call vhaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhaes(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(br(1,1,kk),"  br",nlat,nlat)
c     call aout(bi(1,1,kk),"  bi",nlat,nlat)
c     call aout(cr(1,1,kk),"  cr",nlat,nlat)
c     call aout(ci(1,1,kk),"  ci",nlat,nlat)
c     end do

      call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shes")
      call iout(ierror,"ierr")

      call igrades(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("igra")
      call iout(ierror,"ierr")

      else if (icase .eq. 3) then

      call vhagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhagc(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shgc")
      call iout(ierror,"ierr")

      call igradgc(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("igra")
      call iout(ierror,"ierr")

      else if (icase .eq. 4) then
      call vhagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhags(nlat,nlon,ityp,nt,v,w,nlat,nlon,br,bi,cr,ci,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

      call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shgs")
      call iout(ierror,"ierr")

      call igradgs(nlat,nlon,isym,nt,sf,nlat,nlon,br,bi,
     +mdb,nlat,wsave,lsave,work,lwork,ierror)
      call name("igra")
      call iout(ierror,"ierr")

      end if

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(sf(1,1,kk),"  sf",nlat,nlon)
c     end do

c
c     compare this sf with original
c
      err2s = 0.0
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
	    if (k.eq.1) then
	       se = x*y
	    else if (k.eq.2) then
	       se = x*z
	    else if (k.eq.3) then
	       se = y*z
	    else if (k.eq.4) then
	       se = x*y*z
	    end if
	       err2s = err2s + (sf(i,j,k)-se)**2
	  end do
	end do
      end do
      err2s = sqrt(err2s/(nlat*nlon*nt))
      call vout(err2s,"errs")
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
