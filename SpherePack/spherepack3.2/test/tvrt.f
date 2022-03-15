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
c     6/98
c
c     a program for testing all vorticity and ivnerse vorticity routines
c
c
c     (1) first set a stream function and velocity potential scalar fields as
c         polys in x,y,z restricted to the sphere
c
c     (2) derive a vector field (v,w) from (1)
c
c     (3) compute the vorticity vt of (2) and compare with the vorticity
c         computed analytically
c
c     (4) compute vector field (ve,we) using br,bi,cr,ci from (v,w) with
c         br=bi=0.0
c
c     (5) invert the vorticity in (3) and compare with (4)
c

      program tvrt
c
c     set dimensions with parameter statements
c
      parameter(nnlat= 24,nnlon= 14, nnt = 3)
      parameter (mmdab = (nnlon+2)/2, mmdc = (nnlon+1)/2)
      parameter (lleng= 5*nnlat*nnlat*nnlon,llsav=5*nnlat*nnlat*nnlon)
      dimension work(lleng),wsave(llsav)
      parameter (lldwork = 4*nnlat*nnlat)
      dimension dwork(lldwork)
      dimension br(mmdc,nnlat,nnt),bi(mmdc,nnlat,nnt)
      dimension cr(mmdc,nnlat,nnt),ci(mmdc,nnlat,nnt)
      dimension a(mmdab,nnlat,nnt),b(mmdab,nnlat,nnt)
      dimension vt(nnlat,nnlon,nnt)
      dimension thetag(nnlat),dtheta(nnlat),dwts(nnlat)
      dimension v(nnlat,nnlon,nnt),w(nnlat,nnlon,nnt)
      dimension ve(nnlat,nnlon,nnt),we(nnlat,nnlon,nnt)
      dimension pertrb(nnt)
      double precision dtheta, dwts
c
c     set dimension variables
c
      nlat = nnlat
      nlon = nnlon
      nmax = max0(nlat,nlon)
      mdab = mmdab
      mdc = mmdc
      lwork = lleng
      lsave = llsav
      ldwork = lldwork
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
      call gaqd(nlat,dtheta,dwts,dwork,ldwork,ier)
      do  i=1,nlat
	thetag(i) = dtheta(i)
      end do
      call name("gaqd")
      call iout(ier," ier")
      call vecout(thetag,"thtg",nlat)
c
c     test all vorticity subroutines
c
      do icase=1,4
      call name("****")
      call name("****")
      call iout(icase,"icas")
c
c
c     set scalar stream and velocity potential fields as polys in x,y,z
c     and then set v,w from st,sv scalar fields
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
c              st(i,j,k) = x
c              sv(i,j,k) = y
c
c          v = -1/sin(theta)*dstdp + dsvdt
c
c          w =  1/sin(theta)*dsvdp + dstdt
c
	       v(i,j,k) = sinp + cost*sinp
	       w(i,j,k) = cosp + cost*cosp
	       vt(i,j,k) = -2.0*sint*cosp
	    else if (k.eq.2) then
C              st = y
c              sv = z
	       v(i,j,k) = -cosp-sint
	       w(i,j,k) = cost*sinp
c         sint*vt = -dvdp + sint*dwdt + cost*w
c                 = sinp + sint*(-sint*sinp)+cost*cost*sinp
c                 = sinp + (cost**2-sint**2)*sinp
	       vt(i,j,k) = -2.*sint*sinp
	    else if (k.eq.3) then
c           st = x
c           sv = z
	    v(i,j,k) = sinp - sint
	    w(i,j,k) = cost*cosp
c     sint*vt = -cosp-sint*sint*sinp+cost*cost*cosp
c             = -cosp + (1-2.*sint**2)*cosp =
	    vt(i,j,k) = -2.*sint*cosp
	    end if
	  end do
	end do
      end do

c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     call aout(vt(1,1,kk),"  vt",nlat,nlon)
c     end do

      if (icase.eq.1) then

      call name("**ec")
c
c     analyze vector field
c
      call vhaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhai")
      call iout(ierror,"ierr")

      call vhaec(nlat,nlon,isym,nt,v,w,nlat,nlon,br,bi,cr,ci,mdc,
     +nlat,wsave,lsave,work,lwork,ierror)
      call name("vha ")
      call iout(ierror,"ierr")

c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(br(1,1,kk),"  br",mdc,nlat)
c     call aout(bi(1,1,kk),"  bi",mdc,nlat)
c     call aout(cr(1,1,kk),"  cr",mdc,nlat)
c     call aout(ci(1,1,kk),"  ci",mdc,nlat)
c     end do
c     end if
c
c     compute vorticity of (v,w) in vt
c

      call shseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

      call name("vrti")
      call iout(ierror,"ierr")

      call vrtec(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdc,nlat,
     +wsave,lsave,work,lwork,ierror)

      call name("vrt ")
      call iout(ierror,"ierr")
      call iout(nlat,"nlat")
      call iout(nlon,"nlon")

      else if (icase.eq.2) then

      call name("**es")
      call shsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)

      call name("vrti")
      call iout(ierror,"ierr")

      call vrtes(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdc,nlat,
     +wsave,lsave,work,lwork,ierror)

      call name("vrt ")
      call iout(ierror,"ierr")

      else if (icase .eq. 3) then

      call name("**gc")

      call shsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)

      call name("vrti")
      call iout(ierror,"ierr")

      call vrtgc(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdc,nlat,
     +wsave,lsave,work,lwork,ierror)

      call name("vrt ")
      call iout(ierror,"ierr")

      else if (icase .eq. 4) then

      call name("**gs")

      call shsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)

      call name("vrti")
      call iout(ierror,"ierr")

      call vrtgs(nlat,nlon,isym,nt,vt,nlat,nlon,cr,ci,mdc,nlat,
     +wsave,lsave,work,lwork,ierror)

      call name("vrt ")
      call iout(ierror,"ierr")
      end if

c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(vt(1,1,kk),"  vt",nlat,nlon)
c     end do
c     end if
c
c     compute "error" in vt
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
	       vte = -2.0*sint*cosp
	    else if (k.eq.2) then
	       vte = -2.*sint*sinp
	    else if (k.eq.3) then
	    vte = -2.*sint*cosp
	    end if
	       err2 = err2 + (vt(i,j,k)-vte)**2
	  end do
	end do
      end do
      err2 = sqrt(err2/(nt*nlat*nlon))
      call vout(err2,"err2")
c
c     now recompute (v,w) inverting vt using ivrt(ec,es,gc,gs)
c     and compare with (ve,we) generated by synthesizing br,bi,cr,ci
c     with br=bi=0.0
c

      do k=1,nt
	do i=1,mdc
	  do j=1,nlat
	    br(i,j,k) = 0.0
	    bi(i,j,k) = 0.0
	  end do
	end do
      end do

      if (icase.eq.1) then

      call name("**ec")

c
c     set vector field (ve,we) with br=bi=0.0 for comparison with inverted vt
c
      call vhseci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhsi")
      call iout(ierror,"ierr")

      call vhsec(nlat,nlon,ityp,nt,ve,we,nlat,nlon,br,bi,cr,ci,
     +           mdc,nlat,wsave,lsave,work,lwork,ierror)

      call name("vhs ")
      call iout(ierror,"ierr")

      call shaeci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaec(nlat,nlon,isym,nt,vt,nlat,nlon,a,b,
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
      call name("vhsi")
      call iout(ierror,"ierr")

      call ivrtec(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("ivrt")
      call iout(ierror,"ierr")
      call vout(pertrb,"prtb")

      else if (icase.eq.2) then

      call name("**es")
c
c     set vector field (ve,we) with br=bi=0.0 for comparison with inverted vt
c
      call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhsi")
      call iout(ierror,"ierr")

      call vhses(nlat,nlon,ityp,nt,ve,we,nlat,nlon,br,bi,cr,ci,
     +           mdc,nlat,wsave,lsave,work,lwork,ierror)
      call name("vhs ")
      call iout(ierror,"ierr")


      call shaesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shaes(nlat,nlon,isym,nt,vt,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsesi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("ivti")
      call iout(ierror,"ierr")

      call ivrtes(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("ivrt")
      call iout(ierror,"ierr")
      call vout(pertrb,"prtb")

      else if (icase.eq.3) then

      call name("**gc")
c
c     set vector field (ve,we) with br=bi=0.0 for comparison with inverted vt
c
      call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("vhsi")
      call iout(ierror,"ierr")

      call vhsgc(nlat,nlon,ityp,nt,ve,we,nlat,nlon,br,bi,cr,ci,
     +           mdc,nlat,wsave,lsave,work,lwork,ierror)

      call name("vhs ")
      call iout(ierror,"ierr")

      call shagci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shagc(nlat,nlon,isym,nt,vt,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgci(nlat,nlon,wsave,lsave,dwork,ldwork,ierror)
      call name("ivti")
      call iout(ierror,"ierr")

      call ivrtgc(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)

      call name("ivrt")
      call iout(ierror,"ierr")
      call vout(pertrb,"prtb")

      else if (icase.eq.4) then

      call name("**gs")
c
c     set vector field (ve,we) with br=bi=0.0 for comparison with inverted vt
c
      call vhsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("vhsi")
      call iout(ierror,"ierr")

      call vhsgs(nlat,nlon,ityp,nt,ve,we,nlat,nlon,br,bi,cr,ci,
     +           mdc,nlat,wsave,lsave,work,lwork,ierror)
      call name("vhs ")
      call iout(ierror,"ierr")

      call shagsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("shai")
      call iout(ierror,"ierr")

      call shags(nlat,nlon,isym,nt,vt,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,ierror)
      call name("sha ")
      call iout(ierror,"ierr")
      call iout(lsave,"lsav")
      call iout(lwork,"lwrk")

      call vhsgsi(nlat,nlon,wsave,lsave,work,lwork,dwork,ldwork,ierror)
      call name("ivti")
      call iout(ierror,"ierr")

      call ivrtgs(nlat,nlon,isym,nt,v,w,nlat,nlon,a,b,
     +mdab,nlat,wsave,lsave,work,lwork,pertrb,ierror)
      call name("ivrt")
      call iout(ierror,"ierr")
      call vout(pertrb,"prtb")

      end if


c     if (nmax.lt.10) then
c     do kk=1,nt
c     call iout(kk,"**kk")
c     call aout(v(1,1,kk),"   v",nlat,nlon)
c     call aout(w(1,1,kk),"   w",nlat,nlon)
c     end do
c     end if

c
c     compare this v,w with ve,we
c
      err2v = 0.0
      err2w = 0.0
      do k=1,nt
	do j=1,nlon
	  do i=1,nlat
	    err2v = err2v + (v(i,j,k)-ve(i,j,k))**2
	    err2w = err2w + (w(i,j,k)-we(i,j,k))**2
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
