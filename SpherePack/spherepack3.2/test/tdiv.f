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
c     a program for testing all divergence and inverse divergence routines
c
c     (1) first set an irrotational vector field v,w
c
c     (2) compute the coefficients br,bi,cr,ci of (v,w) using vector analysis
c
c     (3) compute the divergence of (v,w) using divec,dives,divgc,divgs
c
c     (4) analystically compute the divergence and compare with (3)
c
c     (5) invert the divergence and compare with the irrotational (v,w)
c
      program tdiv
c
c     set dimensions with parameter statements
c
      parameter(nnlat= 15,nnlon=  9, nnt = 2)
c *** be sure to set
c     mmdb = min0(nlat,(nlon+1)/2), mmdab = min0(nlat,(nnlon+2)/2)
      parameter (mmdab = (nnlon+2)/2, mmdb = (nnlon+1)/2)
      parameter (lleng= 5*nnlat*nnlat*nnlon,llsav=15*nnlat*nnlat*nnlon)
      parameter (lldwork = 4*nnlat*nnlat)
      dimension work(lleng),wsave(llsav)
      dimension br(mmdb,nnlat,nnt),bi(mmdb,nnlat,nnt)
      dimension dwork(lldwork)
      dimension cr(mmdb,nnlat,nnt),ci(mmdb,nnlat,nnt)
      dimension a(mmdab,nnlat,nnt),b(mmdab,nnlat,nnt)
      dimension dv(nnlat,nnlon,nnt)
      dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
      dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
      dimension pertrb(nnt)
      double precision dtheta, dwts
c
c     set dimension variables
c
      nlat = nnlat
      nlon = nnlon
      nmax = max0(nlat,nlon)
      mdab = mmdab
      mdb = mmdb

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
	    dxdt = cost*cosp
	    d2xdt2 = -sint*cosp
	    dxdp = -sint*sinp
	    d2xdp2 = -sint*cosp
	    dydt = cost*sinp
	    d2ydt2 = -sint*sinp
	    dydp = sint*cosp
	    d2ydp2 = -sint*sinp
	    dzdt = -sint
	    d2zdt2 = -cost
	    dzdp = 0.0
	    d2zdp2 = 0.0
	    if (k.eq.1) then
	       sf = x*y
	       dsfdt = x*dydt+y*dxdt
	       dsfdp = x*dydp+y*dxdp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = cosp*dydp+sinp*dxdp
c              dv = 1/sint*(d(sint*v)/dt + dw/dp)
	       dvdt = x*d2ydt2 + 2.*dxdt*dydt + y*d2xdt2
c              1/sint*dwdp = 1/sint*(cosp*d2ydp2-sinp*dydp+sinp*d2xdp2+cosp*dxdp)
c                          = -4.*sinp*cosp
	       dv(i,j,k) = dvdt + cost*(cosp*dydt+sinp*dxdt) -4.*cosp*sinp
	    else if (k.eq.2) then
	       sf = x*z
	       dsfdt = x*dzdt+z*dxdt
	       dsfdp = x*dzdp+z*dxdp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = -cost*sinp
	       dvdt = x*d2zdt2+2.*dzdt*dxdt + z*dx2dt2
c              dv = 1/sint*(d(sint*v)/dt + dw/dp)
	       dv(i,j,k) = -6.*cost*sint*cosp
	    else if (k.eq.3) then
	       sf = y*z
	       dsfdt = y*dzdt+z*dydt
	       dsfdp = y*dzdp+z*dydp
	       v(i,j,k) = dsfdt
	       w(i,j,k) = z*cosp
	       dv(i,j,k) = -6.*cost*sint*sinp
	    end if
	  end do
	end do
      end do

c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     call aout(dv(1,1,kk),"  dv",nlat,nlon)
c     end do
c     end if

      if (icase.eq.1) then

      call name("**ec")
c
c     analyze vector field
c
      call vhaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shs ")
      call iout(ierror,"ierr")

      call vhaec(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdb,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

cc    if (nmax.lt.10) then
cc    do kk=1,nt
cc    call iout(kk,"**kk")
cc    call aout(br(1,1,kk),"  br",mdb,nlat)
cc    call aout(bi(1,1,kk),"  bi",mdb,nlat)
cc    call aout(cr(1,1,kk),"  cr",mdb,nlat)
cc    call aout(ci(1,1,kk),"  ci",mdb,nlat)
cc    end do
cc    end if
c
c     compute divergence of (v,w) in dv
c

      call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divec(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")
      call iout(nlat,"nlat")
      call iout(nlon,"nlon")


      else if (icase.eq.2) then

      call name("**es")
      call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call dives(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)

      call name("div ")
      call iout(ierror,"ierr")

      else if (icase .eq. 3) then

      call name("**gc")

      call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divgc(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")

      else if (icase .eq. 4) then

      call name("**gs")

      call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shsi")
      call iout(ierror,"ierr")

      call divgs(nlat,nlon,isym,nt,dv,nlat,nlon,br,bi,mdb,nlat,
     +wsave,lsave,work,lwork,ierror)
      call name("div ")
      call iout(ierror,"ierr")

      end if

c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(dv(1,1,kk),"  dv",nlat,nlon)
c     end do
c     end if
c
c     compute "error" in dv
c
      err2 = 0.0
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
	    d2xdt2 = -sint*cosp
	    dxdp = -sint*sinp
	    d2xdp2 = -sint*cosp
	    dydt = cost*sinp
	    d2ydt2 = -sint*sinp
	    dydp = sint*cosp
	    d2ydp2 = -sint*sinp
	    dzdt = -sint
	    d2zdt2 = -cost
	    dzdp = 0.0
	    d2zdp2 = 0.0
	    if (k.eq.1) then
	       sf = x*y
	       dvdt = x*d2ydt2 + 2.*dxdt*dydt + y*d2xdt2
	       dve = dvdt + cost*(cosp*dydt+sinp*dxdt) - 4.*cosp*sinp
	    else if (k.eq.2) then
c              sf = x*z
	       dve = -6.*sint*cost*cosp
	    else if (k.eq.3) then
c              sf = y*z
	       dve = -6.*cost*sint*sinp
	    end if
	       err2 = err2 + (dv(i,j,k)-dve)**2
	  end do
	end do
      end do
c
c     set and print least squares error in dv
c
      err2 = sqrt(err2/(nt*nlat*nlon))
      call vout(err2,"err2")
c
c     now recompute (v,w) inverting dv using idiv(ec,es,gc,gs)
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

      call shaec(nlat,nlon,isym,nt,dv,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(a(1,1,kk),"   a",nlat,nlat)
c     call aout(b(1,1,kk),"   b",nlat,nlat)
c     end do
c     end if

      call vhseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idivec(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("idiv")
      call iout(ierror,"ierr")
      call vecout(pertrb,"prtb",nt)

      else if (icase.eq.2) then

      call name("**es")


      call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaes(nlat,nlon,isym,nt,dv,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idives(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("idiv")
      call iout(ierror,"ierr")
      call vecout(pertrb,"prtb",nt)

      else if (icase.eq.3) then

      call name("**gc")


      call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shagc(nlat,nlon,isym,nt,dv,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idivgc(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("idiv")
      call iout(ierror,"ierr")
      call vecout(pertrb,"prtb",nt)

      else if (icase.eq.4) then

      call name("**gs")


      call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shags(nlat,nlon,isym,nt,dv,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("idvi")
      call iout(ierror,"ierr")

      call idivgs(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("idiv")
      call iout(ierror,"ierr")
      call vecout(pertrb,"prtb",nt)

      end if


c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     end do
c     end if

c
c     compare this v,w with original
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
	       sf = x*y
	       dsfdt = x*dydt+y*dxdt
	       dsfdp = x*dydp+y*dxdp
	       ve = dsfdt
	       we = cosp*dydp+sinp*dxdp
	    else if (k.eq.2) then
	       sf = x*z
	       dsfdt = x*dzdt+z*dxdt
	       dsfdp = x*dzdp+z*dxdp
	       ve = dsfdt
	       we = -cost*sinp
	    else if (k.eq.3) then
	       sf = y*z
	       dsfdt = y*dzdt+z*dydt
	       dsfdp = y*dzdp+z*dydp
	       ve = dsfdt
	       we = z*cosp
	    else if (k.eq.4) then
	    end if
	    err2v = err2v + (v(i,j,k)-ve)**2
	    err2w = err2w + (w(i,j,k)-we)**2
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
   10 format(1h a4, " = ",i8)
      return
      end
c
      subroutine vout(var,nam)
      real nam
      write(6,10) nam , var
   10 format(1h a4," = ",e12.5)
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
