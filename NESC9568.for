
      program control !(input,xoutput,tty,tape59=tty,tape10=input,tape9
   !  1 =xoutput,tape7,tape8,tape2)
c
c
c **********************************************************************
c **********************************************************************
c
copyright,  1987,  the regents of the university of california.
c     this software was produced under a u. s. government contract
c     (w-7405-eng-36) by the los alamos national laboratory, which
c     is operated by the university of california for the u. s.
c     department of energy.  the u. s. government is licensed to use,
c     reproduce, and to distribute this software.  permission is
c     granted to the public to copy and use this software without
c     charge, provided that this notice and any statement of
c     authorship are reproduced on all copies.  neither the government
c     nor the university makes any warranty, express or implied,
c     or assumes any liability or responsibility for the use of this
c     software.
c
c **********************************************************************
c **********************************************************************
c
*ca slcom1
      include "param.fi"
c * * read problem input data
c
      call rinput
c
c * * skip all set-up routines for tape restart run
c
      if (td) 10,10,20
   10 continue
c
c * * read problem set-up data for the mesh, graphics, etc.
c
      call meshset
      call rgrafic
      call rcontur
      call setup
      call flmset
c
c * * the numerical solution alorithm and cyclic control
c     of the calculation is provided by subroutine sola
c
      call sola
!      call exita (405)
c
c * * tape dump restart
c
   20 call sola
!      call exita (100)
c
      end
*dk aset
      subroutine aset
*ca slcom1
      include "param.fi"
c ***
c *** conic fcn=oa2*x*x+oa1*x+ob2*y*y+ob1*y+oc2*x*y+oc1
c *** inside fcn=negative value
c *** ioh=1 add obs inside fcn, ioh=0 subtract obs inside fcn
c ***
      dimension iflg(5), dis(4), xm(5), ym(5)
c
      em6=1.0e-6
      do 20 k=1,kmax
      do 20 j=1,jmax
      do 20 i=1,imax
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      ar(ijk)=1.0
      at(ijk)=1.0
      abk(ijk)=1.0
      ac(ijk)=1.0
      beta(ijk)=0.0
      if (i.eq.im1.and.wr.le.2) ar(ijk)=0.0
      if (k.eq.km1.and.wt.le.2) at(ijk)=0.0
      if (j.eq.jm1.and.wbk.le.2) abk(ijk)=0.0
      if (i.eq.1.and.wl.le.2) go to 10
      if (i.eq.imax.and.wr.le.2) go to 10
      if (k.eq.1.and.wb.le.2) go to 10
      if (k.eq.kmax.and.wt.le.2) go to 10
      if (j.eq.1.and.wf.le.2) go to 10
      if (j.eq.jmax.and.wbk.le.2) go to 10
      go to 20
   10 ar(ijk)=0.0
      at(ijk)=0.0
      abk(ijk)=0.0
      ac(ijk)=0.0
      beta(ijk)=-1.0
   20 continue
c
      if (nobs.le.0) go to 300
c
      j=2
      do 250 n=1,nobs
      do 240 k=2,km1
      do 240 i=2,im1
      call ijkonly
      rdxdy=1.0/(delx(i)*delz(k))
      do 80 m=1,4
      go to (30,40,50,60), m
   30 x1=x(i)
      y1=z(k-1)
      dis(1)=delz(k)
      go to 70
   40 y1=z(k)
      x1=x(i)
      dis(2)=delx(i)
      go to 70
   50 x1=x(i-1)
      y1=z(k)
      dis(3)=delz(k)
      go to 70
   60 y1=z(k-1)
      x1=x(i-1)
      dis(4)=delx(i)
   70 iflg(m)=0
      fconic=oa2(n)*x1*x1+oa1(n)*x1+ob2(n)*y1*y1+ob1(n)*y1+oc2(n)*x1*y1
     1 +oc1(n)
      if (fconic.le.0.0) iflg(m)=1
      xm(m)=x1
      ym(m)=y1
   80 continue
      iflg(5)=iflg(1)
      xm(5)=xm(1)
      ym(5)=ym(1)
      iflgs=0
      do 90 m=1,4
   90 iflgs=iflgs+iflg(m)
      brij=0.0
      btij=0.0
      if (iflgs.eq.0) go to 240
      if (iflgs.lt.4) go to 100
      bij=1.0
      brij=1.0
      btij=1.0
      go to 220
  100 if (iflg(1).eq.1.and.iflg(2).eq.1) brij=1.0
      if (iflg(2).eq.1.and.iflg(3).eq.1) btij=1.0
      do 180 m=1,4
      if (iflg(m).eq.iflg(m+1)) go to 180
      x1=xm(m)
      y1=ym(m)
      x2=xm(m+1)
      y2=ym(m+1)
      if (iflg(m).eq.0) go to 110
      x2=xm(m)
      y2=ym(m)
      x1=xm(m+1)
      y1=ym(m+1)
  110 epsif=0.001*(abs(x2-x1)+abs(y2-y1))
      smn=0.0
      fmn=oa2(n)*x2*x2+oa1(n)*x2+ob2(n)*y2*y2+ob1(n)*y2+oc2(n)*x2*y2+oc1
     1 (n)
      smx=1.0
      fmx=oa2(n)*x1*x1+oa1(n)*x1+ob2(n)*y1*y1+ob1(n)*y1+oc2(n)*x1*y1+oc1
     1 (n)
      s=0.5
  120 xt=s*x1+(1.0-s)*x2
      yt=s*y1+(1.0-s)*y2
      fs=oa2(n)*xt*xt+oa1(n)*xt+ob2(n)*yt*yt+ob1(n)*yt+oc2(n)*xt*yt+oc1
     1 (n)
      if (abs(fs).lt.epsif) go to 150
      if (fs.ge.0.0) go to 130
      fden=abs(fs-fmn)+1.0e-10
      se=s-fs*(s-smn)/fden
      if (se.gt.smx) se=smx
      fmn=fs
      smn=s
      go to 140
  130 fden=abs(fmx-fs)+1.0e-10
      se=s-fs*(smx-s)/fden
      if (se.lt.smn) se=smn
      fmx=fs
      smx=s
  140 si=s-fs*(smx-smn)/(fmx-fmn)
      s=0.5*(se+si)
      go to 120
  150 dis(m)=sqrt((xt-x2)**2+(yt-y2)**2)
      go to (160,170,180,180), m
  160 brij=dis(1)/delz(k)
      go to 180
  170 btij=dis(2)/delx(i)
  180 continue
      m=0
      bij=0.0
  190 continue
      m=m+1
      if (m.eq.5) go to 210
      if (iflg(m).eq.0) go to 190
      mp1=m+1
      if (mp1.eq.5) mp1=1
      mm1=m-1
      if (mm1.eq.0) mm1=4
      bij=bij+dis(m)*dis(mm1)
      if (iflg(mp1).eq.1) go to 200
      dis2=dis(m)
  200 continue
      if (iflg(mm1).eq.1) go to 190
      dis1=dis(mm1)
      go to 190
  210 continue
      if (iflgs.eq.3) bij=bij-dis1*dis2
      bij=0.5*bij*rdxdy
      if (bij.gt.1.0) bij=1.0
  220 continue
      if (ioh(n).eq.0) go to 230
      bij=-bij
      brij=-brij
      btij=-btij
  230 ac(ijk)=ac(ijk)+bij
      if (ac(ijk).gt.0.99) ac(ijk)=1.0
      if (ac(ijk).lt.0.01) ac(ijk)=0.0
      abk(ijk)=ac(ijk)
      ar(ijk)=ar(ijk)+brij
      if (ar(ijk).gt.0.99) ar(ijk)=1.0
      if (ar(ijk).lt.0.01) ar(ijk)=0.0
      if (i.eq.im1.and.wr.le.2) ar(ijk)=0.0
      at(ijk)=at(ijk)+btij
      if (at(ijk).gt.0.99) at(ijk)=1.0
      if (at(ijk).lt.0.01) at(ijk)=0.0
      if (k.eq.km1.and.wt.le.2) at(ijk)=0.0
  240 continue
  250 continue
c
      do 270 k=1,kmax
      if (wl.le.2) go to 260
      ijkl=nq*(ii5*(k-1)+imax+1)+1
      ijklm=ijkl-nq
      ar(ijklm)=ar(ijkl)
      at(ijklm)=at(ijkl)
      abk(ijklm)=abk(ijkl)
      ac(ijklm)=ac(ijkl)
  260 if (wr.le.2) go to 270
      ijkr=nq*(ii5*(k-1)+imax+(im1-1))+1
      ijkrp=ijkr+nq
      ar(ijkrp)=ar(ijkr)
      at(ijkrp)=at(ijkr)
      abk(ijkrp)=abk(ijkr)
      ac(ijkrp)=ac(ijkr)
  270 continue
c
      do 290 i=1,imax
      if (wb.le.2) go to 280
      ijkb=nq*(ii5+imax+(i-1))+1
      ijkbm=ijkb-ii2
      at(ijkbm)=at(ijkb)
      ar(ijkbm)=ar(ijkb)
      abk(ijkbm)=abk(ijkb)
      ac(ijkbm)=ac(ijkb)
  280 if (wt.le.2) go to 290
      ijkt=nq*(ii5*(km1-1)+imax+(i-1))+1
      ijktp=ijkt+ii2
      at(ijktp)=at(ijkt)
      ar(ijktp)=ar(ijkt)
      abk(ijktp)=abk(ijkt)
      ac(ijktp)=ac(ijkt)
  290 continue
c
  300 continue
c
      j=2
      do 310 k=2,km1
      do 310 i=2,im1
      call ijkajct
      if (ac(ijk).gt.em6) go to 310
      ar(ijk)=0.0
      ar(imjk)=0.0
      at(ijk)=0.0
      at(ijkm)=0.0
      abk(ijk)=0.0
      abk(ijmk)=0.0
      beta(ijk)=-1.0
  310 continue
c
c     note: all jplanes have same values
c
      do 360 j=1,jmax
      if (j.eq.2) go to 360
      do 350 k=1,kmax
      do 340 i=1,imax
      i2k=nq*(ii5*(k-1)+imax+(i-1))+1
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      if (j.eq.1.and.wf.le.2) go to 330
      if (j.eq.jmax.and.wbk.le.2) go to 320
      ar(ijk)=ar(i2k)
      at(ijk)=at(i2k)
      abk(ijk)=abk(i2k)
      ac(ijk)=ac(i2k)
      beta(ijk)=beta(i2k)
      go to 340
  320 ijmk=ijk-ii1
      abk(ijmk)=0.0
  330 ar(ijk)=0.0
      at(ijk)=0.0
      abk(ijk)=0.0
      ac(ijk)=0.0
      beta(ijk)=-1.0
  340 continue
  350 continue
  360 continue
c ***
c *** set special values of ar,at,abk here
c ***
      return
      end
*dk bc
      subroutine bc
*ca slcom1
c     
c * * general-surface boundary conditions
c
c * * set left and right bdy conditions
c
      include "param.fi"
      do 150 k=1,kmax
      if (cyl.eq.0.0) go to 20
      sum1=0.0
      sum2=0.0
      sum3=0.0
      do 10 j=2,jm1
      ijk=nq*(ii5*(k-1)+imax*(j-1)+1)+1
      sum1=sum1+u(ijk)*cthj(j)-v(ijk)*sthjbk(j)
      sum2=sum2+u(ijk)*sthj(j)+v(ijk)*cthjbk(j)
      sum3=sum3+w(ijk)
   10 continue
      if (abs(sum1).lt.em6) sum1=0.0
      if (abs(sum2).lt.em6) sum2=0.0
      if (abs(sum3).lt.em6) sum3=0.0
      sum1=0.5*sum1/jbar
      sum2=0.5*sum2/jbar
      sum3=sum3/jbar
      phi0=1.570796327
      if (sum1.ne.0.0) phi0=atan(sum2/sum1)
      cphi0=cos(phi0)
      sphi0=sin(phi0)
      velxy=sqrt(sum1**2+sum2**2)
   20 continue
      do 140 j=1,jmax
      iljk=nq*(ii5*(k-1)+imax*(j-1)+1)+1
      ilmjk=iljk-nq
      ilpjk=iljk+nq
      irjk=iljk+ii0
      irmjk=irjk-nq
      f(ilmjk)=f(iljk)
      p(ilmjk)=p(iljk)
      f(irjk)=f(irmjk)
      p(irjk)=p(irmjk)
      if (cyl.eq.1.0) go to 70
      go to (30,40,50,60,50), wl
   30 u(ilmjk)=0.0
      v(ilmjk)=v(iljk)*(1.0-cyl+cyl*xi(2)*rxi(1))
      w(ilmjk)=w(iljk)
      go to 80
   40 u(ilmjk)=0.0
      v(ilmjk)=-v(iljk)*delxrl
      w(ilmjk)=-w(iljk)*delxrl
      go to 80
   50 if (iter.gt.0) go to 80
      u(ilmjk)=u(iljk)*rr(1)/rr(2)
      v(ilmjk)=v(iljk)
      w(ilmjk)=w(iljk)
      go to 80
   60 u(ilmjk)=u(irmjk)
      v(ilmjk)=v(irmjk)
      w(ilmjk)=w(irmjk)
      f(ilmjk)=f(irmjk)
      nf(ilmjk)=nf(irmjk)
      p(ilmjk)=p(irmjk)
      pn(ilmjk)=pn(irmjk)
      ps(ilmjk)=ps(irmjk)
      peta(ilmjk)=peta(irmjk)
      go to 80
   70 iopjk=nq*(ii5*(k-1)+imax*(jop(j)-1)+1)+1
      f(ilmjk)=f(iopjk)
      p(ilmjk)=p(iopjk)
      u(ilmjk)=velxy*(cthj(j)*cphi0+sthj(j)*sphi0)
      v(ilmjk)=-velxy*(sthjbk(j)*cphi0-cthjbk(j)*sphi0)
      w(ilmjk)=sum3
   80 go to (90,100,110,120,110), wr
   90 u(irmjk)=0.0
      v(irjk)=v(irmjk)*(1.0-cyl+cyl*xi(im1)*rxi(imax))
      w(irjk)=w(irmjk)
      go to 130
  100 u(irmjk)=0.0
      v(irjk)=-v(irmjk)*delxrr
      w(irjk)=-w(irmjk)*delxrr
      go to 130
  110 if (iter.gt.0) go to 130
      u(irjk)=u(irmjk)*rr(imax)/rr(im1)
      v(irjk)=v(irmjk)
      w(irjk)=w(irmjk)
      go to 130
  120 u(irjk)=u(ilmjk)
      v(irjk)=v(ilmjk)
      w(irjk)=w(ilmjk)
      f(irjk)=f(iljk)
      nf(irjk)=nf(iljk)
      p(irjk)=p(iljk)
      pn(irjk)=pn(iljk)
      ps(irjk)=ps(iljk)
      peta(irjk)=peta(iljk)
  130 continue
  140 continue
  150 continue
c
c * * set front and back bdy. conditions
c
      do 270 k=1,kmax
      do 260 i=1,imax
      ifjk=nq*(ii5*(k-1)+imax+i-1)+1
      ifjmk=ifjk-ii1
      ifjpk=ifjk+ii1
      ibkjk=ifjk+ii3
      ibkjmk=ibkjk-ii1
      f(ifjmk)=f(ifjk)
      p(ifjmk)=p(ifjk)
      f(ibkjk)=f(ibkjmk)
      p(ibkjk)=p(ibkjmk)
      go to (160,170,180,190,180), wf
  160 u(ifjmk)=u(ifjk)
      v(ifjmk)=0.0
      w(ifjmk)=w(ifjk)
      go to 200
  170 u(ifjmk)=-u(ifjk)*delyrf
      v(ifjmk)=0.0
      w(ifjmk)=-w(ifjk)*delyrf
      go to 200
  180 if (iter.gt.0) go to 200
      u(ifjmk)=u(ifjk)
      v(ifjmk)=v(ifjk)
      w(ifjmk)=w(ifjk)
      go to 200
  190 u(ifjmk)=u(ibkjmk)
      v(ifjmk)=0.5*(v(ifjmk)+v(ibkjmk))
      if (isor.eq.0) v(ifjmk)=v(ibkjmk)
      w(ifjmk)=w(ibkjmk)
      f(ifjmk)=f(ibkjmk)
      nf(ifjmk)=nf(ibkjmk)
      p(ifjmk)=p(ibkjmk)
      pn(ifjmk)=pn(ibkjmk)
      ps(ifjmk)=ps(ibkjmk)
      peta(ifjmk)=peta(ibkjmk)
  200 go to (210,220,230,240,230), wbk
  210 u(ibkjk)=u(ibkjmk)
      v(ibkjmk)=0.0
      w(ibkjk)=w(ibkjmk)
      go to 250
  220 u(ibkjk)=-u(ibkjmk)*delyrbk
      v(ibkjmk)=0.0
      w(ibkjk)=-w(ibkjmk)*delyrbk
      go to 250
  230 if (iter.gt.0) go to 250
      u(ibkjk)=u(ibkjmk)
      v(ibkjk)=v(ibkjmk)
      w(ibkjk)=w(ibkjmk)
      go to 250
  240 u(ibkjk)=u(ifjk)
      w(ibkjk)=w(ifjk)
      v(ibkjk)=v(ifjk)
      v(ibkjmk)=v(ifjmk)
      f(ibkjk)=f(ifjk)
      nf(ibkjk)=nf(ifjk)
      p(ibkjk)=p(ifjk)
      pn(ibkjk)=pn(ifjk)
      ps(ibkjk)=ps(ifjk)
      peta(ibkjk)=peta(ifjk)
  250 continue
  260 continue
  270 continue
c
c * * set bottom and top bdy. conditions
c
      do 390 j=2,jm1
      do 380 i=2,im1
      ibjk=nq*(ii5+imax*(j-1)+i-1)+1
      ibjkm=ibjk-ii2
      itjk=ibjk+ii4
      itjkm=itjk-ii2
      f(ibjkm)=f(ibjk)
      p(ibjkm)=p(ibjk)
      f(itjk)=f(itjkm)
      p(itjk)=p(itjkm)
      go to (280,290,300,310,300), wb
  280 u(ibjkm)=u(ibjk)
      v(ibjkm)=v(ibjk)
      w(ibjkm)=0.0
      go to 320
  290 u(ibjkm)=-u(ibjk)*delzrb
      v(ibjkm)=-v(ibjk)*delzrb
      w(ibjkm)=0.0
      go to 320
  300 if (iter.gt.0) go to 320
      u(ibjkm)=u(ibjk)
      v(ibjkm)=v(ibjk)
      w(ibjkm)=w(ibjk)
      go to 320
  310 u(ibjkm)=u(itjkm)
      v(ibjkm)=v(itjkm)
      w(ibjkm)=0.5*(w(ibjkm)+w(itjkm))
      f(ibjkm)=f(itjkm)
      nf(ibjkm)=nf(itjkm)
      p(ibjkm)=p(itjkm)
      pn(ibjkm)=pn(itjkm)
      ps(ibjkm)=ps(itjkm)
      peta(ibjkm)=peta(itjkm)
  320 go to (330,340,350,360,350), wt
  330 u(itjk)=u(itjkm)
      v(itjk)=v(itjkm)
      w(itjkm)=0.0
      go to 370
  340 u(itjk)=-u(itjkm)*delzrt
      v(itjk)=-v(itjkm)*delzrt
      w(itjkm)=0.0
      go to 370
  350 if (iter.gt.0) go to 370
      u(itjk)=u(itjkm)
      v(itjk)=v(itjkm)
      w(itjk)=w(itjkm)
      go to 370
  360 u(itjk)=u(ibjk)
      v(itjk)=v(ibjk)
      w(itjk)=w(ibjkm)
      f(itjk)=f(ibjk)
      nf(itjk)=nf(ibjk)
      p(itjk)=p(ibjk)
      pn(itjk)=pn(ibjk)
      ps(itjk)=ps(ibjk)
      peta(itjk)=ps(ibjk)
  370 continue
  380 continue
  390 continue
c
c * * set free surface boundary conditions
c
      call bcfs
c
c * * special boundary conditions
c * * (for specified in or out flow set beta=+1.0 in ficticious cells)
c
      return
      end
*dk bcfs
      subroutine bcfs
*ca slcom1
c
c * * set free surface boundary conditions
c
      include "param.fi"
      do 220 k=2,km1
      do 220 j=2,jm1
      do 220 i=2,im1
      knff=0
      call calcijk
c
c * * set f and p values in obstacle edge cells
c
      if (beta(ijk).ge.0.0) go to 10
      abr=1.0
      if (beta(ipjk).lt.0.0) abr=0.0
      abl=1.0
      if (beta(imjk).lt.0.0) abl=0.0
      abbk=1.0
      if (beta(ijpk).lt.0.0) abbk=0.0
      abf=1.0
      if (beta(ijmk).lt.0.0) abf=0.0
      abt=1.0
      if (beta(ijkp).lt.0.0) abt=0.0
      abb=1.0
      if (beta(ijkm).lt.0.0) abb=0.0
      absum=abr+abl+abbk+abf+abt+abb
      if (absum.lt.0.5) go to 10
      f(ijk)=(abr*f(ipjk)+abl*f(imjk)+abbk*f(ijpk)+abf*f(ijmk)+abt*f
     1 (ijkp)+abb*f(ijkm))/absum
      p(ijk)=(abr*p(ipjk)+abl*p(imjk)+abbk*p(ijpk)+abf*p(ijmk)+abt*p
     1 (ijkp)+abb*p(ijkm))/absum
      go to 220
   10 continue
      if (nmat.eq.2) go to 230
      nff=nf(ijk)
      if (nff.eq.0.or.nff.gt.7) go to 220
c
c * * set velocity on cell boundary between partially
      go to (20,30,40,50,60,70,80), nff
   20 if (ar(ijk).gt.em6) u(ijk)=u(imjk)*ar(imjk)*x(i-1)/(ar(ijk)*x(i))
      if (abk(ijk).gt.em6) v(ijk)=v(imjk)
      if (abk(ijmk).gt.em6) v(ijmk)=v(imjmk)
      if (at(ijk).gt.em6) w(ijk)=w(imjk)
      if (at(ijkm).gt.em6) w(ijkm)=w(imjkm)
      go to 80
   30 if (ar(imjk).gt.em6) u(imjk)=u(ijk)*ar(ijk)*x(i)/(ar(imjk)*x(i-1))
      if (abk(ijk).gt.em6) v(ijk)=v(ipjk)
      if (abk(ijmk).gt.em6) v(ijmk)=v(ipjmk)
      if (at(ijk).gt.em6) w(ijk)=w(ipjk)
      if (at(ijkm).gt.em6) w(ijkm)=w(ipjkm)
      go to 80
   40 if (abk(ijk).gt.em6) v(ijk)=v(ijmk)
      if (ar(ijk).gt.em6) u(ijk)=u(ijmk)
      if (ar(imjk).gt.em6) u(imjk)=u(imjmk)
      if (at(ijk).gt.em6) w(ijk)=w(ijmk)
      if (at(ijkm).gt.em6) w(ijkm)=w(ijmkm)
      go to 80
   50 if (abk(ijk).gt.em6) v(ijmk)=v(ijk)
      if (ar(ijk).gt.em6) u(ijk)=u(ijpk)
      if (ar(imjk).gt.em6) u(imjk)=u(imjpk)
      if (at(ijk).gt.em6) w(ijk)=w(ijpk)
      if (at(ijkm).gt.em6) w(ijkm)=w(ijpkm)
      go to 80
   60 if (at(ijk).gt.em6) w(ijk)=w(ijkm)
      if (ar(ijk).gt.em6) u(ijk)=u(ijkm)
      if (ar(imjk).gt.em6) u(imjk)=u(imjkm)
      if (abk(ijk).gt.em6) v(ijk)=v(ijkm)
      if (abk(ijkm).gt.em6) v(ijmk)=v(ijmkm)
      go to 80
   70 if (at(ijkm).gt.em6) w(ijkm)=w(ijk)
      if (ar(ijk).gt.em6) u(ijk)=u(ijkp)
      if (ar(imjk).gt.em6) u(imjk)=u(imjkp)
      if (abk(ijk).gt.em6) v(ijk)=v(ijkp)
      if (abk(ijmk).gt.em6) v(ijmk)=v(ijmkp)
      go to 80
c
c * * calculate velocity on nf selected cell boundary
c
   80 go to (90,100,110,120,130,140,140), nff
   90 if (nf(ipjk).le.7.or.ar(ijk).lt.em6) go to 150
      denom=-rri(i)*rdx(i)*ar(ijk)/rr(i)
      u(ijk)=(rri(i)*(rdx(i)*(-u(imjk)*ar(imjk)/rr(i-1))+rdy(j)*(v(ijk)
     1 *abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at(ijk)-w(ijkm)*at
     2 (ijkm)))/denom
      go to 160
  100 if (nf(imjk).le.7.or.ar(imjk).lt.em6) go to 150
      denom=rri(i)*rdx(i)*ar(imjk)/rr(i-1)
      u(imjk)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i))+rdy(j)*(v(ijk)*abk
     1 (ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at(ijk)-w(ijkm)*at(ijkm)
     2 ))/denom
      go to 160
  110 if (nf(ijpk).le.7.or.abk(ijk).lt.em6) go to 150
      denom=-rri(i)*rdy(j)*abk(ijk)
      v(ijk)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i
     1 -1))+rdy(j)*(-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at(ijk)-w(ijkm)
     2 *at(ijkm)))/denom
      go to 160
  120 if (nf(ijmk).le.7.or.abk(ijmk).lt.em6) go to 150
      denom=rri(i)*rdy(j)*abk(ijmk)
      v(ijmk)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr
     1 (i-1))+rdy(j)*(v(ijk)*abk(ijk)))+rdz(k)*(w(ijk)*at(ijk)-w(ijkm)
     2 *at(ijkm)))/denom
      go to 160
  130 if (nf(ijkp).le.7.or.at(ijk).lt.em6) go to 150
      denom=-rdz(k)*at(ijk)
      w(ijk)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i
     1 -1))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(-w(ijkm)
     2 *at(ijkm)))/denom
      go to 160
  140 if (nf(ijkm).le.7.or.at(ijkm).lt.em6) go to 150
      denom=rdz(k)*at(ijkm)
      w(ijkm)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr
     1 (i-1))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)
     2 *at(ijk)))/denom
      go to 160
  150 continue
c
c * * if nf designated cell is not empty try all neighboring
c     cells to find empty neighbor
c
      knff=knff+1
      nff=nff+1
      if (nff.gt.7) nff=1
      if (knff.gt.7) go to 160
      go to 80
  160 continue
c
c * * set velocities in empty cells adjacent to partial fluid cells
c
      if (flg.gt.0.5.and.iter.gt.0.and.isor.eq.1) go to 220
      if (nf(ipjk).le.7.or.beta(ipjk).lt.0.0) go to 170
      if (nf(ipjpk).gt.7.and.abk(ipjk).gt.em6) v(ipjk)=f(ijk)*v(ijk)*rri
     1 (i)/rri(i+1)
      if (nf(ipjmk).gt.7.and.abk(ipjmk).gt.em6) v(ipjmk)=f(ijk)*v(ijmk)
     1 *rri(i)/rri(i+1)
      if (nf(ipjkp).gt.7.and.at(ipjk).gt.em6) w(ipjk)=f(ijk)*w(ijk)
      if (nf(ipjkm).gt.7.and.at(ipjkm).gt.em6) w(ipjkm)=f(ijk)*w(ijkm)
  170 if (nf(imjk).le.7.or.beta(imjk).lt.0.0) go to 180
      if (nf(imjpk).gt.7.and.abk(imjk).gt.em6) v(imjk)=f(ijk)*v(ijk)*rri
     1 (i)/rri(i-1)
      if (nf(imjmk).gt.7.and.abk(imjmk).gt.em6) v(imjmk)=f(ijk)*v(ijmk)
     1 *rri(i)/rri(i-1)
      if (nf(imjkp).gt.7.and.at(imjk).gt.em6) w(imjk)=f(ijk)*w(ijk)
      if (nf(imjkm).gt.7.and.at(imjkm).gt.em6) w(imjkm)=f(ijk)*w(ijkm)
  180 if (nf(ijpk).le.7.or.beta(ijpk).lt.0.0) go to 190
      if (nf(imjpk).gt.7.and.ar(imjpk).gt.em6) u(imjpk)=f(ijk)*u(imjk)
      if (nf(ipjpk).gt.7.and.ar(ijpk).gt.em6) u(ijpk)=f(ijk)*u(ijk)
      if (nf(ijpkp).gt.7.and.at(ijpk).gt.em6) w(ijpk)=f(ijk)*w(ijk)
      if (nf(ijpkm).gt.7.and.at(ijpkm).gt.em6) w(ijpkm)=f(ijk)*w(ijkm)
  190 if (nf(ijmk).le.7.or.beta(ijmk).lt.0.0) go to 200
      if (nf(imjmk).gt.7.and.ar(imjmk).gt.em6) u(imjmk)=f(ijk)*u(imjk)
      if (nf(ipjmk).gt.7.and.ar(ijmk).gt.em6) u(ijmk)=f(ijk)*u(ijk)
      if (nf(ijmkp).gt.7.and.at(ijmk).gt.em6) w(ijmk)=f(ijk)*w(ijk)
      if (nf(ijmkm).gt.7.and.at(ijmkm).gt.em6) w(ijmkm)=f(ijk)*w(ijkm)
  200 if (nf(ijkp).le.7.or.beta(ijkp).lt.0.0) go to 210
      if (nf(ipjkp).gt.7.and.ar(ijkp).gt.em6) u(ijkp)=f(ijk)*u(ijk)
      if (nf(imjkp).gt.7.and.ar(imjkp).gt.em6) u(imjkp)=f(ijk)*u(imjk)
      if (nf(ijpkp).gt.7.and.abk(ijkp).gt.em6) v(ijkp)=f(ijk)*v(ijk)
      if (nf(ijmkp).gt.7.and.abk(ijmkp).gt.em6) v(ijmkp)=f(ijk)*v(ijmk)
  210 if (nf(ijkm).le.7.or.beta(ijkm).lt.0.0) go to 220
      if (nf(ipjkm).gt.7.and.ar(ijkm).gt.em6) u(ijkm)=f(ijk)*u(ijk)
      if (nf(imjkm).gt.7.and.ar(imjkm).gt.em6) u(imjkm)=f(ijk)*u(imjk)
      if (nf(ijpkm).gt.7.and.abk(ijkm).gt.em6) v(ijkm)=f(ijk)*v(ijk)
      if (nf(ijmkm).gt.7.and.abk(ijmkm).gt.em6) v(ijmkm)=f(ijk)*v(ijmk)
  220 continue
  230 continue
      return
      end
*dk betacal
      subroutine betacal
*ca slcom1
c
c * * calculate beta(ijk) for mesh interior
c * * note that all ficticious cells have beta set to zero in aset
c * * note that all obstacle cells have beta set to -1.0 in aset
c
c * * calculate beta(ijk) for non-obstacle cells
c
      include "param.fi"
      rijk=0.0
      do 10 k=2,km1
      do 10 j=2,jm1
      do 10 i=2,im1
      call ijkajct
      if (beta(ijk).lt.0.0) go to 10
      rijk=rijk+1.
      abr=ar(ijk)
      if (beta(ipjk).lt.0.0.or.beta(ipjk).eq.1.0) abr=0.0
      abl=ar(imjk)
      if (beta(imjk).lt.0.0.or.beta(imjk).eq.1.0) abl=0.0
      abbk=abk(ijk)
      if (beta(ijpk).lt.0.0.or.beta(ijpk).eq.1.0) abbk=0.0
      abf=abk(ijmk)
      if (beta(ijmk).lt.0.0.or.beta(ijmk).eq.1.0) abf=0.0
      abt=at(ijk)
      if (beta(ijkp).lt.0.0.or.beta(ijkp).eq.1.0) abt=0.0
      abb=at(ijkm)
      if (beta(ijkm).lt.0.0.or.beta(ijkm).eq.1.0) abb=0.0
      xx=2.0*delt*(rdx(i)*(abr/(delx(i)+delx(i+1))+abl/(delx(i)+delx(i-1
     1 )))+rdy(j)*rri(i)*(abbk*rri(i)/(dely(j)+dely(j+1))+abf*rri(i)/
     2 (dely(j)+dely(j-1)))+rdz(k)*(abt/(delz(k)+delz(k+1))+abb/(delz(k)
     3 +delz(k-1)))+cyl*0.5*(abr/(delx(i)+delx(i+1))-abl/(delx(i-1)+delx
     4 (i)))*rri(i)/x(im1))
      beta(ijk)=omg/xx*ac(ijk)
   10 continue
c
      rijk=1.0/rijk
      return
      end
*dk calcijk
      subroutine calcijk
*ca slcom1
c
c * * calculate "ijk" and other indices for cell (i,j,k)
c
      include "param.fi"
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      imjk=ijk-nq
      ipjk=ijk+nq
      ijmk=ijk-ii1
      ijpk=ijk+ii1
      ijkm=ijk-ii2
      ijkp=ijk+ii2
      ipjmk=ijk+nq-ii1
      ipjkm=ijk+nq-ii2
      imjpk=ijk-nq+ii1
      imjmk=ijk-nq-ii1
      ijpkm=ijk+ii1-ii2
      imjkp=ijk-nq+ii2
      ijmkp=ijk-ii1+ii2
      ijmkm=ijk-ii1-ii2
      ipjkp=ijk+nq+ii2
      imjkm=ijk-nq-ii2
      ipjpk=ijk+nq+ii1
      ijpkp=ijk+ii1+ii2
      imjpkp=ijk-nq+ii1+ii2
      ipjmkp=ijk+nq-ii1+ii2
      ipjpkm=ijk+nq+ii1-ii2
      ipjpkp=ijk+nq+ii1+ii2
c
      return
c
      entry ijkonly
c
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      return
c
      entry ijkajct
c
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      imjk=ijk-nq
      ipjk=ijk+nq
      ijmk=ijk-ii1
      ijpk=ijk+ii1
      ijkm=ijk-ii2
      ijkp=ijk+ii2
      return
c
      end
*dk cntr
      subroutine cntr (i1,i2,j1,j2,k1,k2,na,iper)
*ca slcom1
c
c * * contour plotting routine (set-up in subroutine rcontur)
c
      common /sscm10/ mplnc, ifmc, jfmc, kfmc, jnmc, datc, clkc, tc,
     1 iterc, cyclec, namec(10), i1c, i2c, j1c, j2c, k1c, k2c, iperk
     
      include "param.fi"
c
      integer cyclec
c     do 10 i=1,50
c      do 10 j=1,50
c   10 q(i,j)=0.0
c
c * * determine constant plane, mpln
c
c      ifm=0
c      jfm=0
c      kfm=0
c      if (i1.ne.i2) go to 20
c      ifm=i1
c      mpln=1
c   20 if (j1.ne.j2) go to 30
c      jfm=j1
c      mpln=2
c   30 if (k1.ne.k2) go to 40
c      kfm=k1
c      mpln=3
c   40 continue
c
c * * set cell(2,2,2) reference index
c
c      i=2
c      j=2
c      k=2
c      call calcijk
c      ijkc=ijk
c
c * * calculate pressure values, q-array, to be plotted
c
c      do 80 k=2,km1
c      do 80 j=2,jm1
c      do 80 i=2,im1
c      call calcijk
c      iq=i-1
c      jq=j-1
c      kq=k-1
c      go to (50,60,70), mpln
c   50 if (i.ne.i1) go to 80
c      if (j.lt.j1) go to 80
c      if (j.gt.j2) go to 80
c      if (k.lt.k1) go to 80
c      if (k.gt.k2) go to 80
c      q(jq,kq)=p(ijk)-p(ijkc)
c      go to 80
c   60 if (j.ne.j1) go to 80
c      if (k.lt.k1) go to 80
c      if (k.gt.k2) go to 80
c      if (i.lt.i1) go to 80
c      if (i.gt.i2) go to 80
c      q(iq,kq)=p(ijk)-p(ijkc)
c      go to 80
c   70 if (k.ne.k1) go to 80
c      if (j.lt.j1) go to 80
c      if (j.gt.j2) go to 80
c      if (i.lt.i1) go to 80
c      if (i.gt.i2) go to 80
c      q(iq,jq)=p(ijk)-p(ijkc)
c   80 continue
c      mplnc=mpln
c      ifmc=ifm
c      jfmc=jfm
c      kfmc=kfm
c      jnmc=jnm
c      datc=dat
c      clkc=clk
c      tc=t
c      iterc=iter
c      cyclec=cycle
c      do 90 n=1,10
c   90 namec(n)=name(n)
c      i1c=i1
c      i2c=i2
c      j1c=j1
c      j2c=j2
c      k1c=k1
c      k2c=k2
c
c * * call contour plotting routines
c
c      nym=nny+1
c      go to (100,110,120), mpln
c  100 call contrjb (yj(2),nym,zk(2),nnz,q,nzx,nzy,nc,zmn,zmx,dlz,zcq
c     1 ,dmpy,dmpz,igrd,ititle1,ntitle1,ylable,nylbl,zlable,nzlbl)
c      go to 130
c  110 call contrjb (xi(2),nnx,zk(2),nnz,q,nzx,nzy,nc,zmn,zmx,dlz,zcq
c    1 ,dmpx,dmpz,igrd,ititle1,ntitle1,xlable,nxlbl,zlable,nzlbl)
c      go to 130
c  120 call contrjb (xi(2),nnx,yj(2),nym,q,nzx,nzy,nc,zmn,zmx,dlz,zcq
c     1 ,dmpx,dmpy,igrd,ititle1,ntitle1,xlable,nxlbl,ylable,nylbl)
c  130 continue
c
      return
      end
*dk contrjb
      subroutine contrjb (x,nnx,y,nny,z,nzx,nzy,nc,zmn,zmx,dlz,zplan
     1 ,dmpx,dmpy,igrd,ititle,ntitle,xlabel,nxlbl,ylabel,nylbl)
c     include "param.fi"
      common /cntrcom/ isym(50)
      common /sscm10/ mplnc, ifmc, jfmc, kfmc, jnmc, datc, clkc, tc,
     1 iterc, cyclec, namec(10), i1c, i2c, j1c, j2c, k1c, k2c, iperk
      dimension xscale(2), yscale(2)
      equivalence (xmin,xscale(1)), (xmax,xscale(2))
      equivalence (ymin,yscale(1)), (ymax,yscale(2))
      dimension x(1), y(1), z(nzx,1), zplan(1)
      dimension zt(4)
      integer cycle, cyclec
c      noc=min0(iabs(nc),50)
c      zmin=zmn
c      zmax=zmx
c      dlzz=dlz
c      dmapx=dmpx
c      dmapy=dmpy
c      nox=iabs(nnx)
c      noy=iabs(nny)
c      do 10 i=1,50
c   10 isym(i)=0
c
c * * establish scales
c
c      call minv (x,1,nox,i,xmin)
c      call maxv (y,1,noy,i,ymax)
c      call maxv (x,1,nox,i,xmax)
c      call minv (y,1,noy,i,ymin)
c      fgrd=0.
c      if (igrd.gt.0) fgrd=-igrd
c      call plojb (xscale,yscale,2,1,1,1,fgrd,dmapx,dmapy,ititle,-ntitle
c     1 ,xlabel,nxlbl,ylabel,nylbl)
c      if (nc.lt.0) go to 50
c      if (nnx.le.0) call minm (z,nzx,nox,noy,i,j,zmin)
c      if (nny.le.0) call maxm (z,nzx,nox,noy,i,j,zmax)
c      if (dlzz.gt.0) go to 20
c      dlzz=(zmax-zmin)/(noc-1.)
c   20 if (nzy.gt.0) go to 30
c      zmax=zmax-amod(zmax,dlzz)
c      zmin=zmin-amod(zmin,dlzz)
c      noc=min0(noc,ifix((zmax-zmin)/dlzz+1.01))
c   30 zplan(1)=zmin
c      do 40 i=2,noc
c   40 zplan(i)=zplan(i-1)+dlzz
c   50 do 90 ny=2,noy
c      ix=mod(ny,2)
c      dy=y(ny)-y(ny-1)
c      do 80 inx=2,nox
c      nx=inx
c      if (ix.ne.0) nx=nox-inx+2
c      zt(1)=z(nx-1,ny-1)
c      zt(2)=z(nx,ny-1)
c      zt(3)=z(nx,ny)
c      zt(4)=z(nx-1,ny)
c      dx=x(nx)-x(nx-1)
c      if (abs(zt(3)-zt(1))-abs(zt(4)-zt(2))) 70,60,60
c   60 call tricjb (x(nx),y(ny),-dx,-dy,noc,zplan,zt(4),zt(3),zt(2))
c      call tricjb (x(nx-1),y(ny-1),dx,dy,noc,zplan,zt(2),zt(1),zt(4))
c      go to 80
c   70 call tricjb (x(nx-1),y(ny),dx,-dy,noc,zplan,zt(3),zt(4),zt(1))
c      call tricjb (x(nx),y(ny-1),-dx,dy,noc,zplan,zt(1),zt(2),zt(3))
c   80 continue
c   90 continue
c      call lincnt (1)
c      write (12,230) (namec(i),i=1,5),jnmc,datc,clkc
c      call lincnt (59)
c      write (12,240) iterc,tc,cyclec
c      go to (100,110,120), mplnc
c  100 write (12,250) ifmc,j1c,j2c,k1c,k2c
c      go to 130
c  110 write (12,260) jfmc,i1c,i2c,k1c,k2c
c      go to 130
c  120 write (12,270) kfmc,i1c,i2c,j1c,j2c
c  130 continue
c      n=1
c      do 140 i=1,noc
c      ztt=abs(zplan(i))
c      if (ztt.eq.0.or.(ztt.gt..01.and.(ztt.lt.10000))) go to 140
c      n=2
c      go to 150
c  140 continue
c  150 i1=1
c      i2=min0(noc,30)
c  160 call fadv (1)
c      call lincnt (56)
c      write (12,230) (namec(i),i=1,5),jnmc,datc,clkc
c      write (12,240) iterc,tc,cyclec
c     go to (170,180,190), mplnc
c  170 write (12,250) ifmc,j1c,j2c,k1c,k2c
c      go to 200
c  180 write (12,260) jfmc,i1c,i2c,k1c,k2c
c      go to 200
c  190 write (12,270) kfmc,i1c,i2c,j1c,j2c
c  200 continue
c      kys=20
c      call exl
c      call dlch (max0(0,460-ntitle*6),kys,ntitle,ititle,1)
c      kys=kys+60
c      call dlch (250,kys,35,35hidentification        contour value,1)
c      do 220 i=i1,i2
c      kys=kys+30
c      if (n.eq.1) encode (8,280,bcda) zplan(i)
c      if (n.eq.2) encode (8,290,bcda) zplan(i)
c      call exh
c      call drvec (250,kys,400,kys)
c      call exl
c      do 210 k=1,4
c  210 call dlch (k*50+200,kys,0,i,1)
c  220 call dlch (550,kys-6,10,bcda,1)
c      call exh
c      if (i2.eq.noc) return
c      i1=i2+1
c      i2=noc
c      go to 160
c
c  230 format (1h ,2x,5a8,1x,a8,2(1x,a8))
c  240 format (1x,5hiter=,i4,18x,6htime= ,1pe12.5,6x,7hcycle= ,i4)
c  250 format (13h  constant i=,i4,17h   surface     j=,i3,3h to,i3,7h
c     1  k=,i3,4h  to,i3)
c  260 format (13h  constant j=,i4,17h   surface     i=,i3,3h to,i3,7h
c     1  k=,i3,4h  to,i3)
c  270 format (13h  constant k=,i4,17h   surface     i=,i3,3h to,i3,7h
c     1  j=,i3,4h  to,i3)
c  280 format (f8.3)
c  290 format (1pe8.1)
      end
*dk deltadj
      subroutine deltadj
*ca slcom1
c
      include "param.fi"
      data con, itmin, itmost, liter, itcrmx /0.40,15,30,50,200/
c * * adjust time step (delt)
c
      deltn=delt
      if (flgc.lt.0.5) go to 20
c
c * * f convection limit exceeded
c
      t=t-delt
      cycle=cycle-1
      delt=0.5*delt
      do 10 k=1,kmax
      do 10 j=1,jmax
      do 10 i=1,imax
      call calcijk
      f(ijk)=fn(ijk)
      u(ijk)=un(ijk)
      v(ijk)=vn(ijk)
      w(ijk)=wn(ijk)
   10 continue
      nflgc=nflgc+1
   20 continue
      if (autot.lt.0.5.and.fnoc.lt.0.5) go to 90
c
c * * determine maximum velocity term
c
      idumx=1
      jdumx=1
      kdumx=1
      idvmx=1
      jdvmx=1
      kdvmx=1
      idwmx=1
      jdwmx=1
      kdwmx=1
      imx=1
      jmx=1
      kmx=1
      ijkmx=1
      ijkumx=1
      ijkvmx=1
      ijkwmx=1
      dumx=1.0e-10*delx(2)/delt
      dvmx=1.0e-10*dely(2)/rri(2)/delt
      dwmx=1.0e-10*delz(2)/delt
      duvw=amax1(dumx,dvmx,dwmx)
      if (fnoc.gt.0.5) delt=0.5*delt
      do 60 k=2,km1
      do 60 j=2,jm1
      do 60 i=2,im1
      call ijkonly
      if (nf(ijk).ne.0) go to 60
      if (beta(ijk).lt.0.0) go to 60
      udm=abs(un(ijk))/(xi(i+1)-xi(i))
      vdm=abs(vn(ijk))/(yj(j+1)-yj(j))*rri(i)
      wdm=abs(wn(ijk))/(zk(k+1)-zk(k))
      if (udm.lt.dumx) go to 30
      dumx=udm
      idumx=i
      jdumx=j
      kdumx=k
      ijkumx=ijk
   30 if (vdm.lt.dvmx) go to 40
      dvmx=vdm
      jdvmx=j
      idvmx=i
      kdvmx=k
      ijkvmx=ijk
   40 if (wdm.lt.dwmx) go to 50
      dwmx=wdm
      kdwmx=k
      idwmx=i
      jdwmx=j
      ijkwmx=ijk
   50 uvws=amax1(udm,vdm,wdm)
      if (duvw.gt.uvws) go to 60
      duvw=uvws
      imx=i
      jmx=j
      kmx=k
      ijkmx=ijk
   60 continue
c
c * * adjust delt to reflect number of iterations
c
      dtmp=1.02
      if (iter.lt.itmin.and.liter.lt.itmin) dtmp=1.04
      if (isor.ne.1) go to 70
      if (iter.gt.itmost.and.liter.gt.itmost) dtmp=0.99
      go to 80
   70 if (iter.gt.itcrmx.and.liter.gt.itcrmx) dtmp=0.98
   80 delto=delt*dtmp
c
c * * select smallest delt
c
      delt=amin1(delto,con/dumx,con/dvmx,con/dwmx,dtvis)
   90 if (delt.eq.deltn) go to 110
c
c * * adjust beta array to reflect delt change
c
      if (isor.ne.1) go to 110
      do 100 k=1,kmax
      do 100 j=1,jmax
      do 100 i=1,imax
      call ijkonly
      if (beta(ijk).eq.1.0) go to 100
      if (beta(ijk).lt.0.0) go to 100
      beta(ijk)=beta(ijk)*deltn/delt
      if (beta(ijk).eq.1.0) beta(ijk)=1.10
  100 continue
  110 continue
c
      liter=iter
      if (mod(cycle,20).ne.0) return
      write (9,120) cycle,imx,jmx,kmx,ijkmx,u(ijkmx),v(ijkmx),w(ijkmx)
      write (12,120) cycle,imx,jmx,kmx,ijkmx,u(ijkmx),v(ijkmx),w(ijkmx)
      return
c
  120 format (2x,"deltadj",i5,3i4,i6,3(3x,1pe12.5)) !!!!!
      end
*dk dgap
      subroutine dgap (xe,ye,ze,xc1,yc1,zc1)
*ca slcom2
*ca slcom1
c
c * * define graph area for perspective view plots (no film produced)
c
      include "param.fi"
c      xcc=xc1
c      ycc=yc1
c      zcc=zc1
c      yemyc=ye-yc1+1.0e-20
c      xxsq=(xe-xc1)**2
c      yysq=yemyc**2
c      zzsq=(ze-zc1)**2
c      rtxy=sqrt(xxsq+yysq)
c      rhp=1.0/rtxy
c      rhp2=1.0/sqrt(xxsq+yysq+zzsq)
c      csth=-yemyc*rhp
c      snth=(xe-xc1)*rhp
c      trcr=yemyc*csth-(xe-xc1)*snth
c      csphi=rtxy*rhp2
c      snphi=(ze-zc1)*rhp2
c      if (abs(snphi).gt.0.866) go to 10
c      csphi=1.0
c      snphi=0.0
c   10 continue
c      xeb=(xe-xc1)*csth+(ye-yc1)*snth
c      yeb=(zc1-ze)*snphi+trcr*csphi
c      zeb=(ze-zc1)*csphi+trcr*snphi
c      if (izoom.eq.1) go to 30
c      ximx=-1.0e+15
c      etamx=-1.0e+15
c      ximn=+1.0e+15
c      etamn=+1.0e+15
c      do 20 ixr=1,23,3
c      xbnd=grdbn(ixr)
c      ybnd=grdbn(ixr+1)
c      zbnd=grdbn(ixr+2)
c      trcr=(ybnd-yc1)*csth-(xbnd-xc1)*snth
c      xb=(xbnd-xc1)*csth+(ybnd-yc1)*snth
c      yb=(zc1-zbnd)*snphi+trcr*csphi
c      zzb=(zbnd-zc1)*csphi+trcr*snphi
c      yden=amax1((yb-yeb),1.0e-6)
c      yrat=yeb/yden
c      xibnd=xeb-(xb-xeb)*yrat
c      etabnd=zeb-(zzb-zeb)*yrat
c      ximx=amax1(ximx,xibnd)
c      ximn=amin1(ximn,xibnd)
c      etamx=amax1(etamx,etabnd)
c      etamn=amin1(etamn,etabnd)
c   20 continue
c      dxi=ximx-ximn
c      deta=etamx-etamn
c      xietamx=amax1(dxi,deta)
c      gsz=900.0
c      if (lpr.eq.0) gsz=1020.0
c      gdrat=gsz/xietamx
c   30 continue
c
      return
      end
*dk draw
      subroutine draw
*ca slcom1
c
c * * plot velocity vectors, particles, etc. on film
c
c * * grid plot
c
c * * velocity vector and constant plane free surface plot
c
c******* plot a two dimensional vector and surface j=2 and j=jm1
c
      include "param.fi"
c      if (cyl.gt.0.5.and.(abs(y(jm1)/x(im1)-6.2831853).lt.0.01.or.abs(y
c     1 (jm1)/x(im1)-3.14159265).lt.0.01)) call drawq
c
c      ivewo=1
c      call dgap (xea(ivewo),yea(ivewo),zea(ivewo),xca(ivewo),yca(ivewo)
c     1 ,zca(ivewo))
c      if (nvplts.eq.0) go to 30
c      do 20 ni=1,nvplts
c      if (ivvew(ni).eq.ivewo) go to 10
c      ivewo=ivvew(ni)
c      call dgap (xea(ivewo),yea(ivewo),zea(ivewo),xca(ivewo),yca(ivewo)
c     1 ,zca(ivewo))
c   10 call velv (iv1(ni),iv2(ni),jv1(ni),jv2(ni),kv1(ni),kv2(ni),nav(ni)
c     1 ,iperv(ni))
c   20 continue
c   30 continue
c
c * * contour plots
c
c      ivewo=1
c      if (ncplts.eq.0) go to 60
c      do 50 ni=1,ncplts
c      if (icvew(ni).eq.ivewo) go to 40
c      ivewo=icvew(ni)
c      call dgap (xea(ivewo),yea(ivewo),zea(ivewo),xca(ivewo),yca(ivewo)
c     1 ,zca(ivewo))
c   40 call cntr (ic1(ni),ic2(ni),jc1(ni),jc2(ni),kc1(ni),kc2(ni),nac(ni)
c     1 ,iperc(ni))
c   50 continue
c   60 continue
c
c * * perspective surface plots (no hidden lines in this version)
c
c      ivewo=1
c      if (nsplts.eq.0) go to 90
c      do 80 ni=1,nsplts
c      if (isvew(ni).eq.ivewo) go to 70
c      ivewo=isvew(ni)
c      call dgap (xea(ivewo),yea(ivewo),zea(ivewo),xca(ivewo),yca(ivewo)
c     1 ,zca(ivewo))
c   70 call surfplt (is1(ni),is2(ni),js1(ni),js2(ni),ks1(ni),ks2(ni),nas
c     1 (ni),ipers(ni))
c   80 continue
c   90 continue
      return
      end
*dk drf
      subroutine drf (ifm,jfm,kfm)
*ca slcom1
c
c * * drf - draw frame (plane view)
c     drfp -draw frame in perspective (entry point)
c
c * * select constant plane, mpln
c
      include "param.fi"
c      if (ifm.ne.0) mpln=1
c      if (jfm.ne.0) mpln=2
c      if (kfm.ne.0) mpln=3
c      go to (10,20,30), mpln
c
c * * constant i plane
c
c   10 yleft=0.0
c      yright=ybbk-ybf
c      zbot=0.0
c      ztop=zbt-zbb
c      ixl=fixl(mpln)+yleft*xconv(mpln)
c      ixr=fixl(mpln)+yright*xconv(mpln)
c      ixb=fiyb+zbot*yconv(mpln)
c      ixt=fiyb+ztop*yconv(mpln)
c      call frame (ixl,ixr,ixt,ixb)
c      call frame (ixl,ixr,ixt,ixb)
c      go to 50
c
c * * constant j plane
c
c   20 xleft=0.0
c      xright=xbr-xbl
c     zbot=0.0
c      ztop=zbt-zbb
c      ixl=fixl(mpln)+xleft*xconv(mpln)
c      ixr=fixl(mpln)+xright*xconv(mpln)
c      ixb=fiyb+zbot*yconv(mpln)
c      ixt=fiyb+ztop*yconv(mpln)
c      call frame (ixl,ixr,ixt,ixb)
c      call frame (ixl,ixr,ixt,ixb)
c      go to 50
c
c * * constant k plane
c
c   30 if (cyl.eq.1.0) go to 40
c      xleft=0.0
c      xright=xbr-xbl
c      ybot=0.0
c      ytop=ybbk-ybf
c      ixl=fixl(mpln)+xleft*xconv(mpln)
c      ixr=fixl(mpln)+xright*xconv(mpln)
c      ixb=fiyb+ybot*yconv(mpln)
c      ixt=fiyb+ytop*yconv(mpln)
c      call frame (ixl,ixr,ixt,ixb)
c      call frame (ixl,ixr,ixt,ixb)
c      go to 50
c
c * * draw frame for constant k-plane in cylindrical coordinates
c
c   40 continue
c      call drfcyl (mpln)
c
c * * label coordinate axes
c
c   50 call lincnt (59)
c      ix1=amax1(0.0,(fixl(mpln)-8))
c      iy1=fiyb-32
c      ix2=fixl(mpln)+16
c
c      draw obstacle boundaries for jplane plot
c
c      call drwobs (mpln)
c
c     go to (60,70,80), mpln
c
c * * constant i plane
c
c   60 call plt (ix1,iy1,57)
c      call plt (ix2,iy2,56)
c      return
c
c * * constant j plane
c
c   70 call plt (ix1,iy1,57)
c      call plt (ix2,iy2,55)
      return
c
c * * constant k plane
c
c   80 continue
      return
c
      end
*dk drfcyl
      subroutine drfcyl (mpln)
*ca slcom1
c
c * * draw frame for constant k-plane in cylindrical coordinates
c
      include "param.fi"
c      x1=x(1)
c      x2=x(im1)
c      crx=cyl/x2
c      xleft=xblc
c      yleft=ybfc
c      jpl=jm1
c      do 20 j=2,jpl
c      th1=crx*y(j-1)
c      th2=crx*y(j)
c      sth1=sin(th1)
c      cth1=cos(th1)
c      sth2=sin(th2)
c      cth2=cos(th2)
c      xx1=x1*cth1
c      yy1=x1*sth1+y(j-1)*(1.0-cyl)
c      xx2=x1*cth2
c      yy2=x1*sth2+y(j)*(1.0-cyl)
c      xright1=xx1-xleft
c      yright1=yy1-yleft
c      ix1=fixl(mpln)+xright1*xconv(mpln)
c      iy1=fiyb+yright1*yconv(mpln)
c      xright2=xx2-xleft
c      yright2=yy2-yleft
c      ix2=fixl(mpln)+xright2*xconv(mpln)
c      iy2=fiyb+yright2*yconv(mpln)
c      call drv (ix1,iy1,ix2,iy2)
c      dth21=(th2-th1)/10.0
c      th0=th1
c      do 10 l=1,10
c      th01=th0+(l-1)*dth21
c      th02=th0+l*dth21
c      cth1=cos(th01)
c      cth2=cos(th02)
c      sth1=sin(th01)
c      sth2=sin(th02)
c      xx1=x2*cth1
c      xx2=x2*cth2
c      yy1=x2*sth1+y(j-1)*(1.0-cyl)
c      yy2=x2*sth2+y(j)*(1.0-cyl)
c      xright1=xx1-xleft
c      yright1=yy1-yleft
c      ix1=fixl(mpln)+xright1*xconv(mpln)
c     iy1=fiyb+yright1*yconv(mpln)
c      xright2=xx2-xleft
c      yright2=yy2-yleft
c      ix2=fixl(mpln)+xright2*xconv(mpln)
c      iy2=fiyb+yright2*yconv(mpln)
c   10 call drv (ix1,iy1,ix2,iy2)
c   20 continue
c      if (abs(y(jpl)*crx-6.28327853).lt.0.0001) return
c      y1=y(1)
c      y2=y(jm1)
c      th1=crx*y1
c      th2=crx*y2
c      sth1=sin(th1)
c      cth1=cos(th1)
c      sth2=sin(th2)
c      cth2=cos(th2)
c      xx1=x1*cth1
c      yy1=x1*sth1+y1*(1.0-cyl)
c      xx2=x2*cth1
c      yy2=x2*sth1+y1*(1.0-cyl)
c      xright1=xx1-xleft
c      yright1=yy1-yleft
c      ix1=fixl(mpln)+xright1*xconv(mpln)
c      iy1=fiyb+yright1*yconv(mpln)
c      xright2=xx2-xleft
c      yright2=yy2-yleft
c      ix2=fixl(mpln)+xright2*xconv(mpln)
c      iy2=fiyb+yright2*yconv(mpln)
c      call drv (ix1,iy1,ix2,iy2)
c      xx1=x1*cth2
c      yy1=x1*sth2+y2*(1.0-cyl)
c      xx2=x2*cth2
c      yy2=x2*sth2+y2*(1.0-cyl)
c      xright1=xx1-xleft
c      yright1=yy1-yleft
c      ix1=fixl(mpln)+xright1*xconv(mpln)
c      iy1=fiyb+yright1*yconv(mpln)
c      xright2=xx2-xleft
c      yright2=yy2-yleft
c      ix2=fixl(mpln)+xright2*xconv(mpln)
c      iy2=fiyb+yright2*yconv(mpln)
c     call drv (ix1,iy1,ix2,iy2)
c
      return
      end
*dk drfp
      subroutine drfp
*ca slcom1
c
c * * draw frame in perspective
c
      include "param.fi"
c      x1=x(1)
c      x2=x(im1)
c      y1=y(1)
c      y2=y(jm1)
c      z1=z(1)
c      z2=z(km1)
c      crx=cyl/x2
c      th1=crx*y1
c      th2=crx*y2
c      cth1=cos(th1)
c      cth2=cos(th2)
c      sth1=sin(th1)
c      sth2=sin(th2)
c      xx1=x1*cth1
c      yy1=x1*sth1+y1*(1.0-cyl)
c      xx2=x2*cth1
c      yy2=x2*sth1+y1*(1.0-cyl)
c      call pcnv (ixi1,ieta1,xx1,yy1,z1)
c      call pcnv (ixi2,ieta2,xx2,yy2,z1)
c      call drvec (ixi1,ieta1,ixi2,ieta2)
c      call pcnv (ixi3,ieta3,xx1,yy1,z2)
c      call pcnv (ixi4,ieta4,xx2,yy2,z2)
c      call drvec (ixi3,ieta3,ixi4,ieta4)
c      call pcnv (ixi5,ieta5,xx1,yy1,z1)
c      call pcnv (ixi6,ieta6,xx1,yy1,z2)
c      call drvec (ixi5,ieta5,ixi6,ieta6)
c      call pcnv (ixi7,ieta7,xx2,yy2,z1)
c      call pcnv (ixi8,ieta8,xx2,yy2,z2)
c      call drvec (ixi7,ieta7,ixi8,ieta8)
c      xx1=x1*cth2
c      yy1=x1*sth2+y2*(1.0-cyl)
c      xx2=x2*cth2
c      yy2=x2*sth2+y2*(1.0-cyl)
c      call pcnv (ixi1,ieta1,xx1,yy1,z1)
c      call pcnv (ixi2,ieta2,xx2,yy2,z1)
c      call drvec (ixi1,ieta1,ixi2,ieta2)
c      call pcnv (ixi3,ieta3,xx1,yy1,z2)
c      call pcnv (ixi4,ieta4,xx2,yy2,z2)
c      call drvec (ixi3,ieta3,ixi4,ieta4)
c      call pcnv (ixi5,ieta5,xx1,yy1,z1)
c      call pcnv (ixi6,ieta6,xx1,yy1,z2)
c      call drvec (ixi5,ieta5,ixi6,ieta6)
c      call pcnv (ixi7,ieta7,xx2,yy2,z1)
c      call pcnv (ixi8,ieta8,xx2,yy2,z2)
c      call drvec (ixi7,ieta7,ixi8,ieta8)
c      jpl=jm1
c      do 10 j=2,jpl
c      th1=crx*y(j-1)
c      th2=crx*y(j)
c      sth1=sin(th1)
c      cth1=cos(th1)
c      sth2=sin(th2)
c      cth2=cos(th2)
c      xx1=x1*cth1
c      yy1=x1*sth1+y(j-1)*(1.0-cyl)
c      xx2=x1*cth2
c      yy2=x1*sth2+y(j)*(1.0-cyl)
c      call pcnv (ixi1,ieta1,xx1,yy1,z1)
c      call pcnv (ixi2,ieta2,xx2,yy2,z1)
c      call drvec (ixi1,ieta1,ixi2,ieta2)
c      call pcnv (ixi3,ieta3,xx1,yy1,z2)
c      call pcnv (ixi4,ieta4,xx2,yy2,z2)
c      call drvec (ixi3,ieta3,ixi4,ieta4)
c      xx1=x2*cth1
c      xx2=x2*cth2
c      yy1=x2*sth1+y(j-1)*(1.0-cyl)
c      yy2=x2*sth2+y(j)*(1.0-cyl)
c      call pcnv (ixi1,ieta1,xx1,yy1,z1)
c      call pcnv (ixi2,ieta2,xx2,yy2,z1)
c      call drvec (ixi1,ieta1,ixi2,ieta2)
c      call pcnv (ixi3,ieta3,xx1,yy1,z2)
c      call pcnv (ixi4,ieta4,xx2,yy2,z2)
c      call drvec (ixi3,ieta3,ixi4,ieta4)
c   10 continue
c
c * * draw obstacle in perspective
c
c     call drobsp
c
      return
      end
*dk drobsp
      subroutine drobsp
*ca slcom1
c
c * * draw obstacles in perspective
c
      include "param.fi"
c      ipl=im1
c      jpl=jm1
c      kpl=km1
c
c      do 40 i=1,ipl
c      do 40 j=1,jpl
c      do 40 k=1,kpl
c      acs=1.0
c      ars=1.0
c      abks=1.0
c      ats=1.0
c      if (beta(ijk).lt.0.0) acs=0.0
c      if (beta(ipjk).lt.0.0) ars=0.0
c      if (beta(ijpk).lt.0.0) abks=0.0
c      if (beta(ijkp).lt.0.0) ats=0.0
c
c      iplx=0
c      if (acs.eq.0.0) iplx=1
c      if (beta(ijpkp).lt.0.0) iplx=iplx+1
c      ipld=0
c      if (ats.eq.0.0) ipld=1
c      if (abks.eq.0.0) ipld=ipld+1
c      iplx=iplx-ipld
c
c      iply=0
c      if (acs.eq.0.0) iply=1
c      if (beta(ipjkp).lt.0.0) iply=iply+1
c      ipld=0
c      if (ars.eq.0.0) ipld=1
c      if (ats.eq.0.0) ipld=ipld+1
c      iply=iply-ipld
c
c      iplz=0
c      if (acs.eq.0.0) iplz=1
c      if (beta(ipjpk).lt.0.0) iplz=iplz+1
c      ipld=0
c      if (ars.eq.0.0) ipld=1
c      if (abks.eq.0.0) ipld=ipld+1
c      iplz=iplz-ipld
c
c      if (iplx.eq.0.and.iply.eq.0.and.iplz.eq.0) go to 40
c      xb1=x(i)
c      yb1=y(j)
c      zz1=z(k)
c      crx=cyl/x(im1)
c      th1=crx*yb1
c      cth1=cos(th1)
c      sth1=sin(th1)
c      xx1=xb1*cth1
c      yy1=xb1*sth1+yb1*(1.0-cyl)
c      call pcnv (ixi1,ieta1,xx1,yy1,zz1)
c      iret=1
c      if (iplx.eq.0.or.i.eq.1) go to 10
c      xx2=x(i-1)*cth1
c      yy2=x(i-1)*sth1+yb1*(1.0-cyl)
c      zz2=zz1
c      go to 30
c   10 iret=2
c      if (iplz.eq.0.or.k.eq.1.) go to 20
c      xx2=xx1
c      yy2=yy1
c      zz2=z(k-1)
c      go to 30
c   20 iret=3
c      if (iply.eq.0.or.j.eq.1) go to 40
c      yb2=y(j-1)
c      th2=crx*yb2
c      cth2=cos(th2)
c      sth2=sin(th2)
c      xx2=xb1*cth2
c      yy2=xb1*sth2+yb2*(1.0-cyl)
c      zz2=zz1
c   30 call pcnv (ixi2,ieta2,xx2,yy2,zz2)
c      call drvec (ixi1,ieta1,ixi2,ieta2)
c      go to (10,20,40), iret
c   40 continue
c
      return
      end
*dk drwobs
      subroutine drwobs (mpln)
*ca slcom1
c
c     *** draw around all obstacles
c
      include "param.fi"
c      tpi=2.0*pi
c
c      go to (200,10,200), mpln
c
c   10 continue
c      j=2
c      do 190 i=2,im1
c      atr=1.0-cyl+cyl*(1.0-em6)
c      atl=1.0-cyl+cyl*(1.0-em6)
c     atc=1.0-cyl+cyl*(1.0-em6)
c      do 190 k=2,km1
c      call ijkajct
c      if (ac(ijk).lt.em6) go to 190
c      afr=1.0
c      aft=1.0
c      afl=1.0
c      afb=1.0
c      if (ar(ijk).lt.atr) afr=ar(ijk)/atr
c      if (at(ijk).lt.atc) aft=at(ijk)/atc
c      if (ar(imjk).lt.atl) afl=ar(imjk)/atl
c      if (at(ijkm).lt.atc) afb=at(ijkm)/atc
c      if (ac(ijk).ge.atc) go to 130
c      if (i.eq.2) afl=afr-em6
c      if (i.eq.im1) afr=afl-em6
c      if (k.eq.2) afb=aft-em6
c      if (k.eq.km1) aft=afb-em6
c      if ((aft+afb).lt.em6.or.(afl+afr).lt.em6) go to 190
c      m=1
c      amn=afb+afr
c      if ((afr+aft).gt.amn) go to 20
c      m=2
c      amn=afr+aft
c   20 if ((aft+afl).gt.amn) go to 30
c      m=3
c      amn=aft+afl
c   30 if ((afl+afb).gt.amn) go to 40
c      m=4
c   40 go to (50,70,90,110), m
c   50 x1=x(i-1)+aft*delx(i)
c      y1=z(k)
c      if (aft.lt.1.0) go to 60
c      y1=y1-afr*delz(k)
c   60 x2=x(i-1)
c      y2=z(k)-afl*delz(k)
c      if (afl.lt.1.0) go to 170
c      x2=x2+afb*delx(i)
c      go to 170
c   70 x1=x(i-1)
c      y1=z(k-1)+afl*delz(k)
c      if (afl.lt.1.0) go to 80
c      x1=x1+aft*delx(i)
c   80 x2=x(i-1)+afb*delx(i)
c      y2=z(k-1)
c      if (afb.lt.1.0) go to 170
c      y2=y2+afr*delz(k)
c      go to 170
c   90 x1=x(i)-afb*delx(i)
c      y1=z(k-1)
c      if (afb.lt.1.0) go to 100
c      y1=y1+afl*delz(k)
c  100 x2=x(i)
c      y2=z(k-1)+afr*delz(k)
c      if (afr.lt.1.0) go to 170
c      x2=x2-aft*delx(i)
c      go to 170
c  110 x1=x(i)
c      y1=z(k)-afr*delz(k)
c      if (afr.lt.1.0) go to 120
c      x1=x1-afb*delx(i)
c  120 x2=x(i)-aft*delx(i)
c      y2=z(k)
c      if (aft.lt.1.0) go to 170
c      y2=y2-afl*delz(k)
c      go to 170
c  130 if (afr.gt.em6) go to 140
c      x1=x(i)
c      y1=z(k-1)
c      x2=x1
c      y2=z(k)
c      assign 140 to kr1
c      go to 180
c  140 if (aft.gt.em6) go to 150
c      x1=x(i-1)
c      y1=z(k)
c      x2=x(i)
c      y2=y1
c      assign 150 to kr1
c      go to 180
c  150 if (afl.gt.em6) go to 160
c      x1=x(i-1)
c      y1=z(k)
c      x2=x1
c      y2=z(k-1)
c      assign 160 to kr1
c      go to 180
c  160 if (afb.gt.em6) go to 190
c      x1=x(i-1)
c      y1=z(k-1)
c      x2=x(i)
c      y2=y1
c      assign 170 to kr1
c      go to 180
c  170 assign 190 to kr1
c  180 ix1=fixl(mpln)+x1*xconv(mpln)
c      ix2=fixl(mpln)+x2*xconv(mpln)
c      iy1=fiyb+y1*yconv(mpln)
c      iy2=fiyb+y2*yconv(mpln)
c      call drv (ix1,iy1,ix2,iy2)
c      go to kr1, (140,150,160,170,190)
c  190 continue
c  200 continue
      return
      end
*dk drvec
      subroutine drvec (ixi1,ieta1,ixi2,ieta2)
*ca slcom1
c
c * * after clipping,draw vector (a part of zoom feature)
c
c * * working in floating point numbers
c
      include "param.fi"
c      if (iclip.eq.0) go to 90
c      if ((ixi1+ieta1).eq.0.or.(ixi2+ieta2).eq.0) go to 100
c      xil=123.0
c      xir=1023.0
c      if (lpr.eq.0) xil=2.0
c      if (lpr.eq.0) xir=1022.0
c      etat=1.0
c      etab=901.0
c      if (lpr.eq.0) etat=2.0
c      if (lpr.eq.0) etab=1022.0
c      cx1=ixi1
c      x11=cx1
c      cx2=ixi2
c      x22=cx2
c      cy1=ieta1
c      y11=cy1
c      cy2=ieta2
c      y22=cy2
c      dx=x22-x11
c      dy=y22-y11
c      if (abs(dx).lt.1.0e-3) go to 40
c      sm=dy/dx
c      if (cx1.ge.xil) go to 10
c      cy1=y11+sm*(xil-x11)
c      cx1=xil
c   10 if (cx1.le.xir) go to 20
c      cy1=y11+sm*(xir-x11)
c      cx1=xir
c   20 if (cx2.ge.xil) go to 30
c      cy2=y11+sm*(xil-x11)
c      cx2=xil
c   30 if (cx2.le.xir) go to 40
c      cy2=y11+sm*(xir-x11)
c      cx2=xir
c   40 if (cx1.lt.xil.or.cx1.gt.xir) go to 100
c      if (abs(dy).lt.1.0e-3) go to 80
c      sm=dx/dy
c      if (cy1.le.etab) go to 50
c      cx1=x11+sm*(etab-y11)
c      cy1=etab
c   50 if (cy1.ge.etat) go to 60
c      cx1=x11+sm*(etat-y11)
c      cy1=etat
c   60 if (cy2.le.etab) go to 70
c      cx2=x11+sm*(etab-y11)
c      cy2=etab
c   70 if (cy2.ge.etat) go to 80
c      cx2=x11+sm*(etat-y11)
c      cy2=etat
c   80 continue
c      if (cy1.lt.etat.or.cy1.gt.etab) go to 100
c      ds=abs(cx2-cx1)+abs(cy2-cy1)
c      if (ds.lt.1.0e-3) go to 100
c      ixi1=cx1
c      ixi2=cx2
c      ieta1=cy1
c      ieta2=cy2
c   90 call drv (ixi1,ieta1,ixi2,ieta2)
  100 return
      end
*dk equib
      subroutine equib (y,z,nx1,bond,cangle,cyl)
      dimension y(1), z(1)
      data pi /3.141592654/, eps /1.e-04/
c +++
c +++ useful intermediate quantities
c +++
      rangle=cangle
      costst=cos(rangle)
      cosca=(1.+cyl)*cos(rangle)
      omcyl=1.-cyl
      nx=nx1-1
      dr=1./float(nx)
      iter=0
c +++
c +++ set initial guess for y(1)
c +++
      y(1)=0.
      if (abs(bond).gt.0.) y(1)=-cosca/(2.*abs(bond))
c +++
c +++ numerical integration
c +++
   10 z(1)=0.
      iter=iter+1
      z(nx1+1)=0.
      write (2,100) y(1)
      do 20 j=2,nx1
      rjm=dr*float(j-2)*cyl+omcyl
      rj=dr*float(j-1)*cyl+omcyl
      rjh=0.5*(rjm+rj)
      z(j)=z(j-1)*rjm/rj+dr*(rjh/rj)*(cosca-bond*(y(j-1)+0.5*dr*z(j-1)
     1 /sqrt(1.-z(j-1)**2)))
      if (z(j).ge.1.) go to 30
      slope=z(j)/sqrt(1.-z(j)**2)
      slopf=z(j-1)/sqrt(1.-z(j-1)**2)
      z(nx1+j)=(slope-slopf)/(dr*(1.+slope**2)**1.5)
      y(j)=y(j-1)+0.5*dr*(z(j)/sqrt(1.-z(j)**2)+z(j-1)/sqrt(1.-z(j-1)**2
     1 ))
   20 continue
      go to 40
   30 continue
      write (2,110) j,z(j)
      y(1)=y(1)*1.05
      write (2,120) z(nx1)
      if (iter.gt.400) go to 60
      go to 10
c +++
c +++ check constraint and convergence
c +++
   40 continue
      write (2,130)
      ysum=0.5*dr*y(nx1)
      do 50 j=2,nx
      ysum=ysum+(float(j-1)*dr*cyl+omcyl)*dr*y(j)
   50 continue
      write (2,140) ysum,y(1),z(nx1),costst
      if (abs(z(nx1)-cos(rangle)).le.eps.and.bond.ne.0.) go to 70
      if (abs(ysum).lt.eps.and.bond.eq.0.) go to 70
      y(1)=y(1)-ysum/(omcyl+cyl*2.)
      if (iter.gt.400) go to 60
      go to 10
c +++
c +++ exit if convergence fails
c +++
   60 continue
      write (59,150)
      write (12,150)
      write (9,150)
      call exit
   70 continue
      write (2,80) (i,y(i),z(i),z(nx1+i),i=1,nx1)
      write (59,90)
      write (2,90)
      return
c
   80 format (///3x,1hi,8x,1hy,11x,1hz,9x,3hkxy/(1x,i3,1p3e12.4))
   90 format (15h finished equib)
  100 format (1x,17hstep b with y(1)=,1pe14.6)
  110 format (1x,9hz too big,i5,1pe12.4)
  120 format (1x,17hgo to step b   z=,1pe12.4)
  130 format (1x,23hparts b and c converged)
  140 format (1x,5hysum=,1pe13.5,4h y1=,e13.5,4h zn=e13.5,5h cos=,e13.5)
  150 format (1x,28h*** equib failed to converge)
      end
*dk flmset
      subroutine flmset
*ca slcom1
c
c * * calculate mesh boundaries for plotting
c
c * * set boundaries for rectangular mesh
c
      include "param.fi"
c      xbr=x(im1)
c      y2=y(jm1)
c      ybbk=y2
c      ybf=y(1)
c      xbl=x(1)
c      zbt=z(km1)
c      zbb=z(1)
c      if (cyl.lt.0.5) go to 10
c
c * * reset boundaries for cyclindrical coordinates
c
c      xbl=xbr*cos(y2/xbr)
c      if (y2/xbr.gt.0.5*pi) xbl=xbr*cos(y2/xbr)
c      if ((y2/xbr).gt.pi) xbl=-xbr
c      ybbk=xbr*sin(y2/xbr)
c      if ((y2/xbr).gt.(0.5*pi)) ybbk=xbr*sin(y(1)/xbr)
c      if (y2/xbr.ge.pi) ybbk=xbr
c      ybf=x(1)*sin(y(1)/xbr)
c      if ((y2/xbr).lt.0.5*pi) ybf=x(1)*sin(y2/xbr)
c      if (y2/xbr.gt.pi) ybf=-ybbk
c   10 continue
c
c * * set boundary values into grdbn-array for use in dgap subroutine
c
c      grdbn(1)=xbl
c      grdbn(2)=ybf
c      grdbn(3)=zbb
c      grdbn(4)=xbr
c      grdbn(5)=ybf
c      grdbn(6)=zbb
c      grdbn(7)=xbl
c      grdbn(8)=ybbk
c      grdbn(9)=zbb
c      grdbn(10)=xbr
c      grdbn(11)=ybbk
c      grdbn(12)=zbb
c      grdbn(13)=xbl
c      grdbn(14)=ybf
c      grdbn(15)=zbt
c      grdbn(16)=xbr
c      grdbn(17)=ybf
c      grdbn(18)=zbt
c      grdbn(19)=xbl
c      grdbn(20)=ybbk
c      grdbn(21)=zbt
c      grdbn(22)=xbr
c      grdbn(23)=ybbk
c      grdbn(24)=zbt
c
c * * scale mesh boundary values for crt plotting device
c
c * * for k plane
c
c      fiyb=916.0
c      xd=(xbr-xbl)/(ybbk-ybf)
c      yy=0.0
c      if (xd.lt.1.13556) yy=1.0
c      fixl(3)=amax1(0.0,(511.0-450.0*xd)*yy)
c      fixr=(511.0+450.0*xd)*yy+1022.0*(1.0-yy)
c      fiyt=16.0*yy+(916.0-1022.0/xd)*(1.0-yy)
c      xconv(3)=(fixr-fixl(3))/(xbr-xbl)
c      yconv(3)=(fiyt-fiyb)/(ybbk-ybf)
c
c * * reset limits for rectangular mesh for i and j plane
c
c      xblc=xbl
c      ybfc=ybf
c      xbr=x(im1)
c      ybbk=y(jm1)
c      xbl=x(1)
c      ybf=y(1)
c
c * * for j plane
c
c      xd=(xbr-xbl)/(zbt-zbb)
c      yy=0.0
c      if (xd.lt.1.13556) yy=1.0
c      fixl(2)=amax1(0.0,(511.0-450.0*xd)*yy)
c      fixr=(511.0+450.0*xd)*yy+1022.0*(1.0-yy)
c      fiyt=16.0*yy+(916.0-1022.0/xd)*(1.0-yy)
c      xconv(2)=(fixr-fixl(2))/(xbr-xbl)
c      yconv(2)=(fiyt-fiyb)/(zbt-zbb)
c
c * * for i plane
c
c      xd=(ybbk-ybf)/(zbt-zbb)
c      yy=0.0
c      if (xd.lt.1.13556) yy=1.0
c      fixl(1)=amax1(0.0,(511.0-450.0*xd)*yy)
c      fixr=(511.0+450.0*xd)*yy+1022.0*(1.0-yy)
c      fiyt=16.0*yy+(916.0-1022.0/xd)*(1.0-yy)
c      xconv(1)=(fixr-fixl(1))/(ybbk-ybf)
c      yconv(1)=(fiyt-fiyb)/(zbt-zbb)
c
      return
      end
*dk lprt
      subroutine lprt
*ca slcom1
      include "param.fi"
c      if (lpr.eq.0) return
c
c * * print data on film
c
c      call lincnt (64)
c      call fadv (1)
c      assign 30 to kr1
c      go to 20
c   10 lines=lines+1
c      if (lines.lt.55) go to 40 !!!!
c   20 lines=0
c      write (12,70)
c      write (12,100) name,jnm,dat,clk
c      write (12,110) iter,t,delt,cycle
c      write (12,90)
c      go to kr1, (30,40)
c   30 assign 40 to kr1
c      do 50 j=1,jmax
c      do 50 i=1,imax
c      do 50 k=1,kmax
c      call ijkonly
c      write (12,80) i,j,k,u(ijk),v(ijk),w(ijk),p(ijk),ps(ijk),f(ijk),nf
c     1 (ijk),nfp(ijk),nfs(ijk),nfo(ijk),beta(ijk),peta(ijk),d(ijk)
c      go to 10
c   40 continue
c   50 continue
c      if (lpr.le.2) return
c
c * * print data on paper
c
c      entry lprt2
c      write (9,70)
c      write (9,100) name,jnm,dat,clk
c      write (9,110) iter,t,delt,cycle
c      write (9,90)
c      do 60 j=1,jmax
c      do 60 i=1,imax
c      do 60 k=1,kmax
c      call ijkonly
c      write (9,80) i,j,k,u(ijk),v(ijk),w(ijk),p(ijk),ps(ijk),f(ijk),nf
c     1 (ijk),nfp(ijk),nfs(ijk),nfo(ijk),beta(ijk),peta(ijk),d(ijk)
c   60 continue
      return
c
c   70 format (1h1)
c   80 format (3i3,6(1x,1pe11.4),2i3,2i4,3(1x,1pe10.3))
c   90 format (10h  i  j  k,6x,1hu,11x,1hv,11x,1hw,11x,1hp,10x,2hps,10x,1
c     1 hf,6x,2hnf,1x,3hnfp,1x,3hnfs,1x,3hnfo,4x,4hbeta,7x,4hpeta,8x,1hd)
c  100 format (1h ,18x,8a8,1x,a8,2(1x,a8))
c  110 format (5x,6hiter= ,i5,5x,6htime= ,1pe12.5,5x,6hdelt= ,1pe12.5,5x,
c     1 7hcycle= ,i4)
      end
*dk meshset
      subroutine meshset
*ca slcom1
c
c * * mesh generator for x,y,z-coordinates
c
      include "param.fi"
      namelist /meshgn/ nkx,xl,xc,nxl,nxr,dxmn,nky,yl,yc,nyl,nyr,dymn
     1 ,nkz,zl,zc,nzl,nzr,dzmn,nobs,oa2,oa1,ob2,ob1,oc2,oc1,ioh
      data nobs /0/
c
c * * read in data from namelist / meshgn /
c
      read (10,meshgn)
c
c * * write input data for mesh generator onto tapes 9 and 12.
c
      write (9,220) nkx
      write (12,220) nkx
      do 10 l=1,nkx
      write (9,230) l,xl(l),xc(l),xl(l+1),nxl(l),nxr(l),dxmn(l)
      write (12,230) l,xl(l),xc(l),xl(l+1),nxl(l),nxr(l),dxmn(l)
   10 continue
c
      write (9,240) nky
      write (12,240) nky
      do 20 l=1,nky
      write (9,250) l,yl(l),yc(l),yl(l+1),nyl(l),nyr(l),dymn(l)
      write (12,250) l,yl(l),yc(l),yl(l+1),nyl(l),nyr(l),dymn(l)
   20 continue
c
      write (9,260) nkz
      write (12,260) nkz
      do 30 l=1,nkz
      write (9,270) l,zl(l),zc(l),zl(l+1),nzl(l),nzr(l),dzmn(l)
      write (12,270) l,zl(l),zc(l),zl(l+1),nzl(l),nzr(l),dzmn(l)
   30 continue
c
c * * compute x-coordinate values
c
      call meshx (numxp1,nkx)
c
c * * compute y-coordinate values
c
      call meshy (numyp1,nky)
c
c * * compute z-coordinate values
c
      call meshz (numzp1,nkz)
c
      write (9,280) ibar,jbar,kbar
      write (12,280) ibar,jbar,kbar
      write (59,280) ibar,jbar,kbar
c
c * * compute constant terms
c
      imax=ibar+2
      jmax=jbar+2
      kmax=kbar+2
      im1=imax-1
      jm1=jmax-1
      km1=kmax-1
      im2=imax-2
      jm2=jmax-2
      km2=kmax-2
      ii5=imax*jmax
      ii0=nq*ibar
      ii1=nq*imax
      ii2=nq*ii5
      ii3=nq*imax*jbar
      ii4=nq*ii5*kbar
      ii6=imax*jbar
      ii7=imax*jm1
c
c * * set arrays to indefinite
c
      call xindf
c
!      call dateh (dat)
!      call clock1 (clk)
c
c * * initialize radii for cyl. coordinates
c     note, y=x(im1)*theta
c
      do 40 l=2,imax
      rr(l)=1.0-cyl+cyl*x(im1)/x(l)
   40 rri(l)=1.0-cyl+cyl*x(im1)/xi(l)
      rr(1)=1.0
      rri(1)=1.0
      if (x(1).eq.0.0) go to 50
      rr(1)=x(im1)/x(1)
      rri(1)=x(im1)/xi(1)
   50 continue
c
      do 60 j=1,jmax
      sthj(j)=sin(yj(j)/x(im1))
      sthjbk(j)=sin(y(j)/x(im1))
      cthj(j)=cos(yj(j)/x(im1))
      cthjbk(j)=cos(y(j)/x(im1))
   60 continue
      jc2pi=0
      if (cyl.lt.0.5) go to 110
      if (abs(y(jm1)/x(im1)-6.2831853).gt.001) go to 90
      jc2pi=1
      inc=jbar/2
      jlast=inc+1
      jstart=jlast+1
      do 70 j=1,jlast
      jop(j)=j+inc
   70 continue
      do 80 j=jstart,jmax
      jop(j)=j-inc
   80 continue
      go to 110
   90 do 100 j=1,jmax
      jop(j)=jmax-j+1
  100 continue
  110 continue
c * * write mesh generator output
c
      write (9,290)
      write (12,290)
      do 120 n=1,numxp1
      write (9,300) n,x(n),n,rx(n),n,delx(n),n,rdx(n),n,xi(n),n,rxi(n)
      write (12,300) n,x(n),n,rx(n),n,delx(n),n,rdx(n),n,xi(n),n,rxi(n)
  120 continue
      write (9,290)
      write (12,290)
      do 130 n=1,numyp1
      write (9,310) n,y(n),n,dely(n),n,rdy(n),n,yj(n),n,ryj(n)
      write (12,310) n,y(n),n,dely(n),n,rdy(n),n,yj(n),n,ryj(n)
  130 continue
      write (9,290)
      write (12,290)
      do 140 n=1,numzp1
      write (9,320) n,z(n),n,delz(n),n,rdz(n),n,zk(n),n,rzk(n)
      write (12,320) n,z(n),n,delz(n),n,rdz(n),n,zk(n),n,rzk(n)
  140 continue
c
c      calculate partial areas
c
      call aset
c
c     print obstacle data
c
      if (lpr.eq.0) go to 160
      if (nobs.le.0) go to 160
      write (9,210)
      write (12,210)
      do 150 i=1,nobs
      write (9,200) i,oa2(i),oa1(i),ob2(i),ob1(i),oc2(i),oc1(i),ioh(i)
      write (12,200) i,oa2(i),oa1(i),ob2(i),ob1(i),oc2(i),oc1(i),ioh(i)
  150 continue
  160 continue
c * * calc ratios for use in no slip bdy. conditions in subroutine bc
c
      delxrl=delx(1)/delx(2)
      delxrr=delx(imax)/delx(im1)
      delyrf=dely(1)/dely(2)
      delyrbk=dely(jmax)/dely(jm1)
      delzrb=delz(1)/delz(2)
      delzrt=delz(kmax)/delz(km1)
c
c * * determine minimum cell size (length)
c
      delxmn=1.0e+20
      delymn=1.0e+20
      delzmn=1.0e+20
      do 170 n=1,imax
  170 delxmn=amin1(delx(n),delxmn)
      do 180 n=1,jmax
  180 delymn=amin1(dely(n),delymn)
      do 190 n=1,kmax
  190 delzmn=amin1(delz(n),delzmn)
      delmn=amin1(delxmn,delymn,delzmn)
      return
c
  200 format (2x,2hi=,i1,2x,5hoa2= ,1pe12.5,2x,5hoa1= ,e12.5,2x,5hob2=
     1 ,e12.5,2x,5hob1= ,e12.5,2x,5hoc2= ,e12.5,2x,5hoc1= ,e12.5,2x,5hio
     2h= ,i2)
  210 format (///20x,25hconic obstacle parameters)
  220 format (2x,5hnkx= ,i4)
  230 format (2x,8hmesh-x= ,i4,3x,4hxl= ,1pe12.5,3x,4hxc= ,e12.5,3x,4hxr
     1= ,e12.5,3x,5hnxl= ,i4,3x,5hnxr= ,i4,3x,6hdxmn= ,e12.5)
  240 format (2x,5hnky= ,i4)
  250 format (2x,8hmesh-y= ,i4,3x,4hyl= ,1pe12.5,3x,4hyc= ,e12.5,3x,4hyr
     1= ,e12.5,3x,5hnyl= ,i4,3x,5hnyr= ,i4,3x,6hdymn= ,e12.5)
  260 format (2x,5hnkz= ,i4)
  270 format (2x,8hmesh-z= ,i4,3x,4hzl= ,1pe12.5,3x,4hzc= ,e12.5,3x,4hzr
     1= ,e12.5,3x,5hnzl= ,i4,3x,5hnzr= ,i4,3x,6hdzmn= ,e12.5)
  280 format (/,7h  ibar=,i3,8h   jbar=,i3,8h   kbar=,i3)
  290 format (1h1)
  300 format (1x,2hx(,i2,2h)=,1pe12.5,2x,3hrx(,i2,2h)=,1pe12.5,2x,5hdelx
     1(,i2,2h)=,1pe12.5,1x,4hrdx(,i2,2h)=,1pe12.5,2x,3hxi(,i2,2h)=,1pe12
     2 .5,2x,4hrxi(,i2,2h)=,1pe12.5)
  310 format (1x,2hy(,i2,2h)=,1pe12.5,3x,5hdely(,i2,2h)=,1pe12.5,3x,4hrd
     1y(,i2,2h)=,1pe12.5,3x,3hyj(,i2,2h)=,1pe12.5,3x,4hryj(,i2,2h)=,1pe1
     2 2.5)
  320 format (1x,2hz(,i2,2h)=,1pe12.5,3x,5hdelz(,i2,2h)=,1pe12.5,3x,4hrd
     1z(,i2,2h)=,1pe12.5,3x,3hzk(,i2,2h)=,1pe12.5,3x,4hrzk(,i2,2h)=,1pe1
     2 2.5)
      end
*dk meshx
      subroutine meshx (numxp1,nkx)
*ca slcom1
c
c * * compute x-coordinate values and reciprocals
c
      include "param.fi"
      x(1)=xl(1)
      i11=1
      do 30 l=1,nkx
      dxml=(xc(l)-xl(l))/nxl(l)
      dxmn1=dxmn(l)
      nt=nxl(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      dxmn(l)=amin1(dxmn1,dxml)
      cmc=(xc(l)-xl(l)-tn*dxmn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=xc(l)-xl(l)-cmc
      do 10 l1=1,nt
      i11=i11+1
      rln=(float(l1)-tn)/tn
   10 x(i11)=xc(l)+bmc*rln-cmc*rln*rln
      dxmr=(xl(l+1)-xc(l))/nxr(l)
      dxmn(l)=amin1(dxmn1,dxmr)
      nt=nxr(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      cmc=(xl(l+1)-xc(l)-tn*dxmn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=xl(l+1)-xc(l)-cmc
      do 20 l1=1,nt
      i11=i11+1
      rln=float(l1)/tn
   20 x(i11)=xc(l)+bmc*rln+cmc*rln*rln
   30 continue
      numx=i11
      numxm1=numx-1
      numxp1=numx+1
      ibar=numx-1
      do 40 i11=2,numx
   40 delx(i11)=x(i11)-x(i11-1)
      delx(1)=delx(2)
      delx(numxp1)=delx(numx)
      x(numxp1)=x(numx)+delx(numxp1)
c
      do 60 i11=1,numxp1
      if (x(i11).eq.0.0) go to 50
      rx(i11)=1.0/x(i11)
      go to 60
   50 rx(i11)=0.0
   60 continue
      do 70 i11=2,numxp1
      xi(i11)=0.5*(x(i11-1)+x(i11))
      rxi(i11)=1.0/xi(i11)
      if (i11.lt.numxp1) rdxp(i11)=2.0/(delx(i11)+delx(i11+1))
   70 rdx(i11)=1.0/delx(i11)
      xi(1)=xi(2)-delx(2)
      rxi(1)=1.0/xi(1)
      rdx(1)=1.0/delx(1)
      rdxp(1)=2.0/(delx(1)+delx(2))
c
      return
      end
*dk meshy
      subroutine meshy (numyp1,nky)
*ca slcom1
c
c * * compute y-coordinate values and reciprocals
c
      include "param.fi"
      if (nky.eq.1.and.(nyl(1).eq.0.or.nyr(1).eq.0)) go to 40
      y(1)=yl(1)
      j11=1
      do 30 l=1,nky
      dyml=(yc(l)-yl(l))/nyl(l)
      dymn1=dymn(l)
      nt=nyl(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      dymn(l)=amin1(dymn1,dyml)
      cmc=(yc(l)-yl(l)-tn*dymn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=yc(l)-yl(l)-cmc
      do 10 l1=1,nt
      j11=j11+1
      rln=(float(l1)-tn)/tn
   10 y(j11)=yc(l)+bmc*rln-cmc*rln*rln
      dymr=(yl(l+1)-yc(l))/nyr(l)
      dymn(l)=amin1(dymn1,dymr)
      nt=nyr(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      cmc=(yl(l+1)-yc(l)-tn*dymn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=yl(l+1)-yc(l)-cmc
      do 20 l1=1,nt
      j11=j11+1
      rln=float(l1)/tn
   20 y(j11)=yc(l)+bmc*rln+cmc*rln*rln
   30 continue
      go to 50
   40 j11=2
      y(j11)=1.0
      y(j11-1)=0.0
   50 continue
      numy=j11
      numym1=numy-1
      numyp1=numy+1
      jbar=numy-1
      do 60 j11=2,numy
   60 dely(j11)=y(j11)-y(j11-1)
      dely(1)=dely(2)
      dely(numyp1)=dely(numy)
      y(numyp1)=y(numy)+dely(numyp1)
c
      do 80 j11=1,numyp1
      if (y(j11).eq.0.0) go to 70
      ry(j11)=1.0/y(j11)
      go to 80
   70 ry(j11)=0.0
   80 continue
      do 90 j11=2,numyp1
      yj(j11)=0.5*(y(j11-1)+y(j11))
      ryj(j11)=1.0/yj(j11)
      if (j11.lt.numyp1) rdyp(j11)=2.0/(dely(j11)+dely(j11+1))
   90 rdy(j11)=1.0/dely(j11)
      yj(1)=yj(2)-dely(2)
      ryj(1)=1.0/yj(1)
      rdy(1)=1.0/dely(1)
      rdyp(1)=2.0/(dely(1)+dely(2))
c
      return
      end
*dk meshz
      subroutine meshz (numzp1,nkz)
*ca slcom1
c
c * * compute z-coordinate values and reciprocals
c
      include "param.fi"
      z(1)=zl(1)
      k11=1
      do 30 l=1,nkz
      dzml=(zc(l)-zl(l))/nzl(l)
      dzmn1=dzmn(l)
      nt=nzl(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      dzmn(l)=amin1(dzmn1,dzml)
      cmc=(zc(l)-zl(l)-tn*dzmn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=zc(l)-zl(l)-cmc
      do 10 l1=1,nt
      k11=k11+1
      rln=(float(l1)-tn)/tn
   10 z(k11)=zc(l)+bmc*rln-cmc*rln*rln
      dzmr=(zl(l+1)-zc(l))/nzr(l)
      dzmn(l)=amin1(dzmn1,dzmr)
      nt=nzr(l)
      tn=nt
      tn=amax1(tn,1.0+1.0e-14)
      cmc=(zl(l+1)-zc(l)-tn*dzmn(l))*tn/(tn-1.0)
      if (nt.eq.1) cmc=0.0
      bmc=zl(l+1)-zc(l)-cmc
      do 20 l1=1,nt
      k11=k11+1
      rln=float(l1)/tn
   20 z(k11)=zc(l)+bmc*rln+cmc*rln*rln
   30 continue
      numz=k11
      numzm1=numz-1
      numzp1=numz+1
      kbar=numz-1
      do 40 k11=2,numz
   40 delz(k11)=z(k11)-z(k11-1)
      delz(1)=delz(2)
      delz(numzp1)=delz(numz)
      z(numzp1)=z(numz)+delz(numzp1)
c
      do 60 k11=1,numzp1
      if (z(k11).eq.0.0) go to 50
      rz(k11)=1.0/z(k11)
      go to 60
   50 rz(k11)=0.0
   60 continue
      do 70 k11=2,numzp1
      zk(k11)=0.5*(z(k11-1)+z(k11))
      rzk(k11)=1.0/zk(k11)
      if (k11.lt.numzp1) rdzp(k11)=2.0/(delz(k11)+delz(k11+1))
   70 rdz(k11)=1.0/delz(k11)
      zk(1)=zk(2)-delz(2)
      rzk(1)=1.0/zk(1)
      rdz(1)=1.0/delz(1)
      rdzp(1)=2.0/(delz(1)+delz(2))
c
      return
      end
*dk pcnv
      subroutine pcnv (ixi,ieta,xx1,yy1,zz1)
*ca slcom1
*ca slcom2
c
c * * perspective convert subr - produces 4020 coordinates for
c     perspective views
c
      include "param.fi"
      ixi=0
      ieta=0
      trcr=(yy1-ycc)*csth-(xx1-xcc)*snth
      xb=(xx1-xcc)*csth+(yy1-ycc)*snth
      yb=(zcc-zz1)*snphi+trcr*csphi
      zzb=(zz1-zcc)*csphi+trcr*snphi
      yden=yb-yeb
      if (yden.lt.1.0e-6) go to 10
      yrat=yeb/yden
      xid=xeb-(xb-xeb)*yrat
      eta=zeb-(zzb-zeb)*yrat
      ish=61
      if (lpr.eq.0) ish=0
      xia=0.5*(ximx+ximn)
      etaa=0.5*(etamx+etamn)
      ieta=512-(eta-etaa)*gdrat-ish
      ixi=512+(xid-xia)*gdrat+ish
   10 continue
c
      return
      end
*dk pdfcalc
*dk petacal
      subroutine petacal
*ca slcom1
c
c * * determine provisional nf values, surface tension pressures, and
c * * pressure interpolation factor peta for surface and neighbor cells
c
      include "param.fi"
      do 10 k=1,kmax
      do 10 j=1,jmax
      do 10 i=1,imax
      call calcijk
      nfo(ijk)=nf(ijk)
      nf(ijk)=0
      peta(ijk)=1.0
      ps(ijk)=0.0
      nfp(ijk)=0
      nfs(ijk)=0
   10 continue
c
      call pcal
c     calc surface tension
      if (isurft.eq.1) call surf10n
c
      do 100 k=2,km1
      do 100 j=2,jm1
      do 100 i=2,im1
      call ijkajct
      if (beta(ijk).lt.0.0) go to 100
      l=i
      m=j
      n=k
      nfc=nf(ijk)+1
      go to (100,20,30,40,50,60,70,80,100), nfc
   20 dsur=delx(i)
      dnbr=delx(i-1)
      l=i-1
      go to 90
   30 dsur=delx(i)
      dnbr=delx(i+1)
      l=i+1
      go to 90
   40 dsur=dely(j)
      dnbr=dely(j-1)
      m=j-1
      go to 90
c      radial distance factor cancels
   50 dsur=dely(j)
      dnbr=dely(j+1)
      m=j+1
      go to 90
   60 dsur=delz(k)
      dnbr=delz(k-1)
      n=k-1
      go to 90
   70 dsur=delz(k)
      dnbr=delz(k+1)
      n=k+1
      go to 90
   80 p(ijk)=0.1666667*(p(ipjk)+p(ijpk)+p(ijkp)+p(imjk)+p(ijmk)+p(ijkm))
      go to 100
   90 continue
c
c * * calculate peta
c
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      sdis=f(ijk)*dsur+f(lmn)*dnbr*0.5
      dcnt=0.5*(dsur+dnbr)
      sdis=amax1(sdis,0.5*dcnt)
      peta(ijk)=dcnt/sdis
      if (beta(lmn).lt.0.0.or.nf(lmn).ne.0) peta(ijk)=1.0
      if (peta(ijk).gt.2.0) peta(ijk)=2.0
      if (peta(ijk).lt.0.0) peta(ijk)=0.0
  100 continue
c
c * * calculate peta in neighboring interpolation cells
c
      call petaset
c
      if (cycle.gt.0) call presck
c
      return
      end
*dk pcal
      subroutine pcal
*ca slcom1
c
c * * determine provisional nf values
c
      include "param.fi"
      do 10 k=1,kmax
      do 10 j=1,jmax
      do 10 i=1,imax
      call ijkonly
      nf(ijk)=7
      if (i.eq.1.or.i.eq.imax.or.j.eq.1.or.j.eq.jmax.or.k.eq.1.or.k.eq
     1 .kmax.or.beta(ijk).lt.0.0) nf(ijk)=0
   10 continue
      do 160 k=2,km1
      do 160 j=2,jm1
      do 160 i=2,im1
      call calcijk
      if (beta(ijk).lt.0.0) go to 20
      if (f(ijk).gt.emf1.or.ac(ijk).lt.emf) go to 20
      if (f(ijk).gt.emf) go to 30
      nf(ijk)=8
      nfp(ijk)=nf(ijk)
      go to 160
   20 nf(ijk)=0
      nfp(ijk)=nf(ijk)
      go to 160
c
c     looking only at partially filled cells to find surface cells
c      determine nf value for provisional surface orientation
c
   30 fxm=f(ijk)
      fxp=f(ijk)
      fzm=f(ijk)
      fzp=f(ijk)
      fyp=f(ijk)
      fym=f(ijk)
      if (cyl.eq.1.0.and.i.eq.2) imjk=nq*(ii5*(k-1)+imax*(jop(j)-1)+(i-1
     1 ))+1
      if (ar(imjk).gt.emf) fxm=f(imjk)
      if (ar(ijk).gt.emf) fxp=f(ipjk)
      if (at(ijkm).gt.emf) fzm=f(ijkm)
      if (at(ijk).gt.emf) fzp=f(ijkp)
      if (abk(ijk).gt.emf) fyp=f(ijpk)
      if (abk(ijmk).gt.emf) fym=f(ijmk)
c
c      initialize 3 indices for nf decisions and 3 variables for
c      volume of fluid calculations
c
      mobs=1
      inf=1
      iobs=1
      vf=0.0
      vfxm=0.0
      vfxp=0.0
c
c      set up do loops for volume of fluid calculations at xm and xp
c
      do 50 kk=1,3
      n=k-2+kk
      do 40 jj=1,3
      m=j-2+jj
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(i-1))+1
      lmnm=lmn-nq
      lmnp=lmn+nq
      vfxm=vfxm+f(lmnm)
      vfxp=vfxp+f(lmnp)
   40 continue
   50 continue
c
c     calculate nf=1 and reset indices if appropiate
c
      if (beta(imjk).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 60
      if (fxm.lt.emf) go to 60
      inf=inf+1
      if (vfxm.gt.vf) nf(ijk)=1
      if (nf(ijk).eq.1) vf=vfxm
   60 continue
      iobs=1
c
c      calculate nf=2 and reset indices if appropiate
c
      if (beta(ipjk).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 70
      if (fxp.lt.emf) go to 70
      inf=inf+1
      if (vfxp.gt.vf) nf(ijk)=2
      if (nf(ijk).eq.2) vf=vfxp
   70 continue
      iobs=1
c
c      initialize 2 variables for volume of fluid calculations
c
      vfzm=0.0
      vfzp=0.0
c
c      set up do loops for volume of fluid calculations at zm and zp
c
      do 90 ii=1,3
      l=i-2+ii
      do 80 jj=1,3
      m=j-2+jj
      lmn=nq*(imax*jmax*(k-1)+imax*(m-1)+(l-1))+1
      lmnm=lmn-ii2
      lmnp=lmn+ii2
      vfzm=vfzm+f(lmnm)
      vfzp=vfzp+f(lmnp)
   80 continue
   90 continue
c
c     calculate nf=5 and reset indices if appropiate
c
      if (beta(ijkm).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 100
      if (fzm.lt.emf) go to 100
      inf=inf+1
      if (vfzm.gt.vf) nf(ijk)=5
      if (nf(ijk).eq.5) vf=vfzm
  100 continue
      iobs=1
c
c      calculate nf=6 and reset indices if appropiate
c
      if (beta(ijkp).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 110
      if (fzp.lt.emf) go to 110
      inf=inf+1
      if (vfzp.gt.vf) nf(ijk)=6
      if (nf(ijk).eq.6) vf=vfzp
  110 continue
      iobs=1
c
c      initialize 2 variables for volume of fluid calculations
c
      vfym=0.0
      vfyp=0.0
c
c      set up do loops for volume of fluid calculations at ym and yp
c
      do 130 kk=1,3
      n=k-2+kk
      do 120 ii=1,3
      l=i-2+ii
      lmn=nq*(imax*jmax*(n-1)+imax*(j-1)+(l-1))+1
      lmnm=lmn-ii1
      lmnp=lmn+ii1
      vfym=vfym+f(lmnm)
      vfyp=vfyp+f(lmnp)
  120 continue
  130 continue
c
c     calculate nf=3 and reset indices if appropiate
c
      if (beta(ijmk).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 140
      if (fym.lt.emf) go to 140
      inf=inf+1
      if (vfym.gt.vf) nf(ijk)=3
      if (nf(ijk).eq.3) vf=vfym
  140 continue
      iobs=1
c
c      calculate nf=4 and reset indices if appropiate
c
      if (beta(ijpk).lt.0.0) iobs=2
      mobs=mobs+(iobs-1)
      if (iobs.eq.2) go to 150
      if (fyp.lt.emf) go to 150
      inf=inf+1
      if (vfyp.gt.vf) nf(ijk)=4
      if (nf(ijk).eq.4) vf=vfyp
  150 continue
c
c      if we do not reset nf value default correctly gives isolated
c      cell : nf=7
c
c      if we do reset nf value we will have correct nf value - except
c      when there is fluid in each nhbr cell that is not an obstacle
c      cell
c
c      in this exceptional case we should have nf=0
c
c      this is done by defining a critical value for inf and
c      testing on that value
c
      infcr=8-mobs
      if (inf.eq.infcr.and.infcr.gt.1) nf(ijk)=0
c
      nfp(ijk)=nf(ijk)
c
  160 continue
      return
      end
*dk petaset
      subroutine petaset
*ca slcom1
c
c * * calculate peta in neighboring interpolation cells
c
      include "param.fi"
      do 90 k=1,kmax
      do 90 j=1,jmax
      do 90 i=1,imax
      call calcijk
      nff=nf(ijk)
      if (nff.eq.0.or.beta(ijk).lt.0.0) go to 90
      if (nff.gt.7) go to 80
      l=i
      m=j
      n=k
      go to (10,20,30,40,50,60,90), nff
   10 l=i-1
      amn=ar(imjk)
      dnbr=delx(l)
      dcnt=0.5*(dnbr+delx(i))
      go to 70
   20 l=i+1
      amn=ar(ijk)
      dnbr=delx(l)
      dcnt=0.5*(dnbr+delx(i))
      go to 70
   30 m=j-1
      amn=abk(ijmk)
      dnbr=dely(m)/rri(i)
      dcnt=0.5*(dnbr+dely(j)/rri(i))
      go to 70
   40 m=j+1
      amn=abk(ijk)
      dnbr=dely(m)/rri(i)
      dcnt=0.5*(dnbr+dely(j)/rri(i))
      go to 70
   50 n=k-1
      amn=at(ijkm)
      dnbr=delz(n)
      dcnt=0.5*(dnbr+delz(k))
      go to 70
   60 n=k+1
      amn=at(ijk)
      dnbr=delz(n)
      dcnt=0.5*(dnbr+delz(k))
   70 continue
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).gt.0.or.amn.lt.em6) go to 90
      bpd=1.0/peta(lmn)-beta(lmn)*(1.0-peta(ijk))*amn/ac(lmn)*delt/(dnbr
     1 *dcnt)
      peta(lmn)=amin1(1.0/bpd,1.98/omg)
      go to 90
   80 continue
c
c * * set void region nff pressure into cell ijk
c
      p(ijk)=pr(nff)
   90 continue
c
      return
      end
*dk presck
      subroutine presck
c
*ca slcom1
      include "param.fi"
      data psat /0.0/
c
      do 150 k=2,km1
      do 150 j=2,jm1
      do 150 i=2,im1
      call ijkajct
      if (ac(ijk).lt.em6.or.beta(ijk).le.0.0) go to 130
      nff=nf(ijk)
      nffo=nfo(ijk)
      if (nff.eq.nffo) go to 150
      if (nff.eq.0) go to 140
      if (nff.gt.7) go to 100
      l=i
      m=j
      n=k
      go to (10,20,30,40,50,60,70), nff
   10 l=i-1
      go to 80
   20 l=i+1
      go to 80
   30 m=j-1
      go to 80
   40 m=j+1
      go to 80
   50 n=k-1
      go to 80
   60 n=k+1
      go to 80
   70 go to 90
   80 lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      nfel=nf(imjk)
      nfer=nf(ipjk)
      nfef=nf(ijmk)
      nfebk=nf(ijpk)
      nfeb=nf(ijkm)
      nfet=nf(ijkp)
      nfe=max0(nfel,nfer,nfebk,nfef,nfeb,nfet,1)
      psurf=ps(ijk)+pr(nfe)
      plmn=p(lmn)
      if (nf(lmn).ne.0.and.beta(ijk).gt.0.0) plmn=psurf
      p(ijk)=(1.0-peta(ijk))*plmn+peta(ijk)*psurf
      go to 150
   90 if (psat.le.0.0) go to 120
      if (f(ijk).lt.emf1) go to 110
      pmps=0.0
      if (pn(ijk).lt.psat) pmps=p(ijk)-psat
      dijk=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i-1
     1 ))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at
     2 (ijk)-w(ijkm)*at(ijkm)))/ac(ijk)
      dijk=dijk-pmps**2/delt
      p(ijk)=pn(ijk)-dijk*beta(ijk)/(1.0-2.0*pmps*beta(ijk)/delt)
      go to 150
  100 p(ijk)=pr(nff)
      go to 150
  110 p(ijk)=psat
      go to 150
  120 p(ijk)=pn(ijk)
      go to 150
  130 p(ijk)=0.0
  140 continue
  150 continue
c
      if (cycle.lt.9999) return
c
      do 160 k=2,km1
      do 160 j=2,jm1
      do 160 i=2,im1
      call ijkajct
      if (ac(ijk).lt.em6.or.beta(ijk).le.0.0) go to 160
      nff=nf(ijk)
      nffo=nfo(ijk)
      if (nff.eq.nffo) go to 160
      write (9,170) cycle,i,j,k,nf(ijk),nfo(ijk),d(ijk),p(ijk),pn(ijk),u
     1 (ijk),u(imjk),v(ijk),v(ijmk),w(ijk),w(ijkm)
  160 continue
      return
c
c  170 format (* nf switch  after petacal  *,6i4/9(2x,1pe12.5))
  170 format (" nf switch  after petacal  ",6i4/9(2x,1pe12.5))
      end
*dk pltpt
      subroutine pltpt (xone,yone,ichar,isym)
*ca slcom1
      include "param.fi"
      parameter (ibar2q=2*ibar2-2)
      common betaq(ibar2q,kbar2), fq(ibar2q,kbar2), acq(ibar2q,kbar2),
     1 arq(ibar2q,kbar2), atq(ibar2q,kbar2), uq(ibar2q,kbar2), vq(ibar2q
     2 ,kbar2), xq(ibar2q), xiq(ibar2q), yq(kbar2), yjq(kbar2), im1q,
     3 jm1q, jbar2q, jmaxq, imaxq, sf, xshft, yshft, delxq(ibar2q),
     4 delyq(kbar2), xminq, xmaxq, yminq, ymaxq
c
c +++ plot (draw) a point
c +++ provides a system dependant call
c
      ic=0
      x1=xone
      y1=yone
      x01=(x1-xminq)*sf+xshft
      y01=(y1-yminq)*sf+yshft
      ix1=16.+900.0*x01
      iy1=16.+900.0*(1.0-y01)
!      call plt (ix1,iy1,42)
      if (abs(x1).le.em6) go to 10
   10 return
      end
*dk prescr
      subroutine prescr
*ca slcom1
c
      include "param.fi"
      dimension ff(ibar2*jbar2*kbar2), af(ibar2*jbar2*kbar2), ad(ibar2
     1 *jbar2*kbar2), qq(ibar2*jbar2*kbar2), vlm(ibar2*jbar2*kbar2), ssm
     2 (ibar2*jbar2*kbar2), rxr(ibar2*jbar2*kbar2), rxl(ibar2*jbar2
     3 *kbar2), rzt(ibar2*jbar2*kbar2), rzb(ibar2*jbar2*kbar2), rybk
     4 (ibar2*jbar2*kbar2), ryf(ibar2*jbar2*kbar2), div(ibar2*jbar2
     5 *kbar2), dp(ibar2*jbar2*kbar2), cq(ibar2*jbar2*kbar2)
c
      data adefm, bdefm /100.0,0.10/
      data itmax, itmin, itmost /500,5,90/
      data ktran /1/
      data psat /0.0/
c
c     this subroutine must be modified for boundary conditions
c     other than normal velocity equal to zero
c
      go to (10,200), ktran
   10 ktran=2
c
      lvec=ii5*kmax
      if (lvec.ne.ibar2*jbar2*kbar2) go to 760
      epsip=2.0*epsi
      epsim=0.5*epsi
      do 20 k=1,kmax
      do 20 j=1,jmax
      do 20 i=1,imax
      call ijkonly
      lijk=ii5*(k-1)+imax*(j-1)+i
      vlm(lijk)=ac(ijk)*delx(i)*delz(k)*dely(j)/rri(i)
   20 continue
c
      if (cycle.gt.1) go to 200
      if (gz.eq.0.0) go to 40
      do 30 j=2,jm1
      do 30 i=2,im1
      psadd=0.0
      do 30 kk=2,km1
      k=km1-kk+2
      call ijkonly
      if (ps(ijk).ne.0.0) psadd=ps(ijk)-p(ijk)
      p(ijk)=p(ijk)+psadd
   30 continue
c
   40 do 190 k=1,kmax
      do 190 j=1,jmax
      do 190 i=1,imax
      call ijkajct
      lijk=ii5*(k-1)+imax*(j-1)+i
      vlm(lijk)=ac(ijk)*delx(i)*delz(k)*dely(j)/rri(i)
      if (ac(ijk).lt.em6.or.beta(ijk).le.0.0) go to 170
      nff=nf(ijk)
      if (nff.eq.0) go to 180
      if (nff.gt.7) go to 140
      l=i
      m=j
      n=k
      go to (50,60,70,80,90,100,110), nff
   50 l=i-1
      go to 120
   60 l=i+1
      go to 120
   70 m=j-1
      go to 120
   80 m=j+1
      go to 120
   90 n=k-1
      go to 120
  100 n=k+1
      go to 120
  110 go to 130
  120 lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      nfel=nf(imjk)
      nfer=nf(ipjk)
      nfef=nf(ijmk)
      nfebk=nf(ijpk)
      nfeb=nf(ijkm)
      nfet=nf(ijkp)
      nfe=max0(nfel,nfer,nfebk,nfef,nfeb,nfet,1)
      psurf=ps(ijk)+pr(nfe)
      plmn=p(lmn)
      if (nf(lmn).ne.0.and.beta(ijk).gt.0.0) plmn=psurf
      p(ijk)=(1.0-peta(ijk))*plmn+peta(ijk)*psurf
      go to 180
  130 if (psat.le.0.0) go to 160
      if (f(ijk).lt.emf1) go to 150
      pmps=0.0
      if (pn(ijk).lt.psat) pmps=p(ijk)-psat
      dijk=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i-1
     1 ))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at
     2 (ijk)-w(ijkm)*at(ijkm)))/ac(ijk)
      dijk=dijk-pmps**2/delt
      p(ijk)=pn(ijk)-dijk*beta(ijk)/(1.0-2.0*pmps*beta(ijk)/delt)
      go to 180
  140 p(ijk)=pr(nff)
      go to 180
  150 p(ijk)=psat
      go to 180
  160 p(ijk)=pn(ijk)
      go to 180
  170 p(ijk)=0.0
  180 pn(ijk)=p(ijk)
  190 continue
c
      call bc
c
  200 if (cycle.lt.1) return
c
      c1=1.0
      dd=1.0
      do 360 k=1,kmax
      do 360 j=1,jmax
      do 360 i=1,imax
      call ijkajct
      lijk=ii5*(k-1)+imax*(j-1)+i
      if (ac(ijk).lt.em6.or.beta(ijk).le.0.0) go to 340
      if (nf(ijk).gt.0) go to 350
      dijk=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i-1
     1 ))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at
     2 (ijk)-w(ijkm)*at(ijkm)))/ac(ijk)
c
      if (idefm.eq.0.or.f(ijk).ge.emf1) go to 210
c
      dijk=dijk+amin1(adefm*epsi,bdefm*(1.0-f(ijk))/delt)
c
  210 div(lijk)=dijk*vlm(lijk)
c
      ssm(lijk)=0.0
      rxr(lijk)=0.0
      rxra=0.0
      rxrb=0.0
      if (ar(ijk).lt.em6) go to 230
      rxr(lijk)=delz(k)*dely(j)*rdxp(i)*ar(ijk)/rr(i)
      ssm(lijk)=ssm(lijk)+delt*rxr(lijk)
      if (nf(ipjk).eq.0) go to 230
      if (nf(ipjk).gt.7) go to 220
      if (nf(ipjk).ne.1) go to 220
      rxra=(1.0-peta(ipjk))*rxr(lijk)
      rxrb=((1.0-peta(ipjk))*pn(ijk)+peta(ipjk)*(ps(ipjk)+psat)-pn(ipjk)
     1 )*rxr(lijk)
      go to 220
  220 rxr(lijk)=0.0
  230 rxl(lijk)=0.0
      rxla=0.0
      rxlb=0.0
      if (ar(imjk).lt.em6) go to 250
      rxl(lijk)=delz(k)*dely(j)*rdxp(i-1)*ar(imjk)/rr(i-1)
      ssm(lijk)=ssm(lijk)+delt*rxl(lijk)
      if (nf(imjk).eq.0) go to 250
      if (nf(imjk).gt.7) go to 240
      if (nf(imjk).ne.2) go to 240
      rxla=(1.0-peta(imjk))*rxl(lijk)
      rxlb=((1.0-peta(imjk))*pn(ijk)+peta(imjk)*(ps(imjk)+psat)-pn(imjk)
     1 )*rxl(lijk)
      go to 240
  240 rxl(lijk)=0.0
  250 rzt(lijk)=0.0
      rzta=0.0
      rztb=0.0
      if (at(ijk).lt.em6) go to 270
      rzt(lijk)=delx(i)*dely(j)/rri(i)*rdzp(k)*at(ijk)
      ssm(lijk)=ssm(lijk)+delt*rzt(lijk)
      if (nf(ijkp).eq.0) go to 270
      if (nf(ijkp).gt.7) go to 260
      if (nf(ijkp).ne.5) go to 260
      rzta=(1.0-peta(ijkp))*rzt(lijk)
      rztb=((1.0-peta(ijkp))*pn(ijk)+peta(ijkp)*(ps(ijkp)+psat)-pn(ijkp)
     1 )*rzt(lijk)
      go to 260
  260 rzt(lijk)=0.0
  270 rzb(lijk)=0.0
      rzba=0.0
      rzbb=0.0
      if (at(ijkm).lt.em6) go to 290
      rzb(lijk)=delx(i)*dely(j)/rri(i)*rdzp(k-1)*at(ijkm)
      ssm(lijk)=ssm(lijk)+delt*rzb(lijk)
      if (nf(ijkm).eq.0) go to 290
      if (nf(ijkm).gt.7) go to 280
      if (nf(ijkm).ne.6) go to 280
      rzba=(1.0-peta(ijkm))*rzb(lijk)
      rzbb=((1.0-peta(ijkm))*pn(ijk)+peta(ijkm)*(ps(ijkm)+psat)-pn(ijkm)
     1 )*rzb(lijk)
      go to 280
  280 rzb(lijk)=0.0
  290 rybk(lijk)=0.0
      rybka=0.0
      rybkb=0.0
      if (abk(ijk).lt.em6) go to 310
      rybk(lijk)=delx(i)*delz(k)*rdyp(j)*rri(i)*abk(ijk)
      ssm(lijk)=ssm(lijk)+delt*rybk(lijk)
      if (nf(ijpk).eq.0) go to 310
      if (nf(ijpk).gt.7) go to 300
      if (nf(ijpk).ne.3) go to 300
      rybka=(1.0-peta(ijpk))*rybk(lijk)
      rybkb=((1.0-peta(ijpk))*pn(ijk)+peta(ijpk)*(ps(ijpk)+psat)-pn(ijpk
     1 ))*rybk(lijk)
      go to 300
  300 rybk(lijk)=0.0
  310 ryf(lijk)=0.0
      ryfa=0.0
      ryfb=0.0
      if (abk(ijmk).lt.em6) go to 330
      ryf(lijk)=delx(i)*delz(k)*rdyp(j-1)*rri(i)*abk(ijmk)
      ssm(lijk)=ssm(lijk)+delt*ryf(lijk)
      if (nf(ijmk).eq.0) go to 330
      if (nf(ijmk).gt.7) go to 320
      if (nf(ijmk).ne.4) go to 320
      ryfa=(1.0-peta(ijmk))*ryf(lijk)
      ryfb=((1.0-peta(ijmk))*pn(ijk)+peta(ijmk)*(ps(ijmk)+psat)-pn(ijmk)
     1 )*ryf(lijk)
      go to 320
  320 ryf(lijk)=0.0
  330 ssm(lijk)=ssm(lijk)-delt*(rxra+rxla+rzta+rzba+rybka+ryfa)
      div(lijk)=(div(lijk)-delt*(rxrb+rxlb+rztb+rzbb+rybkb+ryfb))/ssm
     1 (lijk)
      ff(lijk)=0.0
      af(lijk)=0.0
      qq(lijk)=0.0
      ad(lijk)=0.0
      dp(lijk)=0.0
      go to 360
  340 vlm(lijk)=1.0
  350 div(lijk)=0.0
      ff(lijk)=0.0
      af(lijk)=0.0
      ad(lijk)=0.0
      qq(lijk)=0.0
      ssm(lijk)=1.0
      rxr(lijk)=0.0
      rxl(lijk)=0.0
      rzt(lijk)=0.0
      rzb(lijk)=0.0
      rybk(lijk)=0.0
      ryf(lijk)=0.0
      dp(lijk)=0.0
  360 continue
c
      if (wbk.ne.4) go to 380
      do 370 k=1,kmax
      do 370 i=1,imax
      lijkf=ii5*(k-1)+i
      lijkfp=lijkf+imax
      lijkbk=lijkfp+imax*jbar
      lijkbkm=lijkbk-imax
      ssm(lijkf)=ssm(lijkbkm)
      vlm(lijkf)=vlm(lijkbkm)
      div(lijkf)=div(lijkbkm)
      rxr(lijkf)=rxr(lijkbkm)
      rxl(lijkf)=rxl(lijkbkm)
      rzt(lijkf)=rzt(lijkbkm)
      rzb(lijkf)=rzb(lijkbkm)
      rybk(lijkf)=rybk(lijkbkm)
      ryf(lijkf)=ryf(lijkbkm)
      ssm(lijkbk)=ssm(lijkfp)
      vlm(lijkbk)=vlm(lijkfp)
      div(lijkbk)=div(lijkfp)
      rxr(lijkbk)=rxr(lijkfp)
      rxl(lijkbk)=rxl(lijkfp)
      rzt(lijkbk)=rzt(lijkfp)
      rzb(lijkbk)=rzb(lijkfp)
      rybk(lijkbk)=rybk(lijkfp)
      ryf(lijkbk)=ryf(lijkfp)
  370 continue
c
c  **********
  380 continue
c  **********
  390 a=c1/dd
!      call saxpy (lvec,a,ff,1,dp,1)
!      call saxpy (lvec,-a,qq,1,div,1)
c
      do 400 k=1,kmax
      kkk=ii5*(k-1)
      do 400 j=1,jmax
      jjj=kkk+imax*(j-1)
      do 400 i=1,imax
c
c     the following commented statements replaced in order to vectorize
c     ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
c     lijk=ii5*(k-1)+imax*(j-1)+i
c     cq(lijk)=div(lijk)*ssm(lijk)/vlm(lijk)
c     if((nf(ijk).ne.0).or.(beta(ijk).lt.0.0))cq(lijk)=0.0
c
      cq(jjj+i)=div(jjj+i)*ssm(jjj+i)/vlm(jjj+i)
      if ((nf(nq*(jjj+i-1)+1).ne.0).or.(beta(nq*(jjj+i-1)+1).lt.0.0)) cq
     1 (jjj+i)=0.0
  400 continue
c
!      ijmax=isamax(lvec,cq,1)
      dmax=abs(cq(ijmax))
      if (dmax.le.epsi) go to 450
c
      do 410 k=1,kmax
      kkk=ii5*(k-1)
      do 410 j=1,jmax
      jjj=kkk+imax*(j-1)
      do 410 i=1,imax
c
c     the following commented statements replaced in order to vectorize
c     ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
c     lijk=ii5*(k-1)+imax*(j-1)+i
c     lipjk=lijk+1
c     limjk=lijk-1
c     lijpk=lijk+imax
c     lijmk=lijk-imax
c     lijkp=lijk+ii5
c     lijkm=lijk-ii5
c     ad(lijk)=delt*(rxr(lijk)*div(lipjk)+rxl(lijk)*div(limjk)+rzt(lijk)
c    1 *div(lijkp)+rzb(lijk)*div(lijkm)+rybk(lijk)*div(lijpk)+ryf(lijk)*
c    2 div(lijmk))-ssm(lijk)*div(lijk)
c     if((nf(ijk).ne.0).or.(beta(ijk).lt.0.0))ad(lijk)=0.0
c
      ad(jjj+i)=delt*(rxr(jjj+i)*div(jjj+i+1)+rxl(jjj+i)*div(jjj+i-1)
     1 +rzt(jjj+i)*div(jjj+i+ii5)+rzb(jjj+i)*div(jjj+i-ii5)+rybk(jjj+i)
     2 *div(jjj+i+imax)+ryf(jjj+i)*div(jjj+i-imax))-ssm(jjj+i)*div(jjj+i
     3 )
      if ((nf(nq*(jjj+i-1)+1).ne.0).or.(beta(nq*(jjj+i-1)+1).lt.0.0)) ad
     1 (jjj+i)=0.0
  410 continue
c
      if (wbk.ne.4) go to 430
      do 420 k=1,kmax
      kkk=ii5*(k-1)
cdir$ ivdep
      do 420 i=1,imax
c
c     the following commented statements replaced in order to vectorize
c     lijkf=ii5*(k-1)+i
c     lijkfp=lijkf+imax
c     lijkbk=lijkfp+imax*jbar
c     lijkbkm=lijkbk-imax
c     ad(lijkf)=ad(lijkbkm)
c     ad(lijkbk)=ad(lijkfp)
c
      ad(kkk+i)=ad(kkk+ii6+i)
      ad(kkk+ii7+i)=ad(kkk+imax+i)
  420 continue
c
  430 c2=1 !sdot(lvec,div,1,ad,1)
      b=c2/c1
!      call sscal (lvec,b,ff,1)
!      call saxpy (lvec,1.0,div,1,ff,1)
c
!      call sscal (lvec,b,af,1)
!      call saxpy (lvec,1.0,ad,1,af,1)
c
      do 440 k=1,kmax
      do 440 j=1,jmax
      do 440 i=1,imax
      lijk=ii5*(k-1)+imax*(j-1)+i
      qq(lijk)=af(lijk)/ssm(lijk)
  440 continue
c
      dd=1 !sdot(lvec,qq,1,af,1)
      c1=c2
      iter=iter+1
      if (iter.gt.itmax) go to 750
      go to 390
c
  450 if (iter.gt.itmost) epsi=1.05*epsi
      if (iter.lt.itmin) epsi=0.95*epsi
      epsi=amin1(epsi,epsip)
      epsi=amax1(epsi,epsim)
c
c  **********
c
c +++ unfold pressure for interior fluid cells
c
      do 460 k=1,kmax
      kkk=ii5*(k-1)
      do 460 j=1,jmax
      jjj=kkk+imax*(j-1)
      do 460 i=1,imax
      d(nq*(jjj+i-1)+1)=cq(jjj+i)
      if ((nf(nq*(jjj+i-1)+1).le.0).and.(ac(nq*(jjj+i-1)+1).ge.em6).and.
     1 (beta(nq*(jjj+i-1)+1).ge.0.0)) p(nq*(jjj+i-1)+1)=pn(nq*(jjj+i-1)+
     2 1)+dp(jjj+i)
  460 continue
c
c +++ calculate pressure changes for surface cells
c
      do 610 k=1,kmax
      do 610 j=1,jmax
      do 610 i=1,imax
      call ijkajct
      lijk=ii5*(k-1)+imax*(j-1)+i
      if (ac(ijk).lt.em6.or.beta(ijk).lt.0.0) go to 590
      nff=nf(ijk)
      if (nff.eq.0) go to 600
      if (nff.gt.7) go to 560
      l=i
      m=j
      n=k
      go to (470,480,490,500,510,520,530), nff
  470 l=i-1
      go to 540
  480 l=i+1
      go to 540
  490 m=j-1
      go to 540
  500 m=j+1
      go to 540
  510 n=k-1
      go to 540
  520 n=k+1
      go to 540
  530 go to 550
  540 lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      nfel=nf(imjk)
      nfer=nf(ipjk)
      nfef=nf(ijmk)
      nfebk=nf(ijpk)
      nfeb=nf(ijkm)
      nfet=nf(ijkp)
      nfe=max0(nfel,nfer,nfebk,nfef,nfeb,nfet,1)
      psurf=ps(ijk)+pr(nfe)
      plmn=p(lmn)
      if (nf(lmn).ne.0.and.beta(ijk).gt.0.0) plmn=psurf
      dp(lijk)=(1.0-peta(ijk))*plmn+peta(ijk)*psurf-pn(ijk)
      go to 600
  550 if (psat.le.0.0) go to 580
      if (f(ijk).lt.emf1) go to 570
      pmps=0.0
      if (pn(ijk).lt.psat) pmps=p(ijk)-psat
      dijk=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i-1
     1 ))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)*at
     2 (ijk)-w(ijkm)*at(ijkm)))/ac(ijk)
      dijk=dijk-pmps**2/delt
      dp(lijk)=-dijk*beta(ijk)/(1.0-2.0*pmps*beta(ijk)/delt)
      go to 600
  560 dp(lijk)=pr(nff)-pn(ijk)
      go to 600
  570 dp(lijk)=psat-pn(ijk)
      go to 600
  580 dp(lijk)=0.0
      go to 600
  590 dp(lijk)=0.0
  600 continue
  610 continue
c
c      calc new velocities consistent with surface pressure assumptions
c     note: yields div=0.0 for interior cells adjacent to surface cells
c
      do 730 k=1,km1
      do 730 j=1,jm1
      do 730 i=1,im1
      call ijkajct
      lijk=ii5*(k-1)+imax*(j-1)+i
      lipjk=lijk+1
      lijpk=lijk+imax
      lijkp=lijk+ii5
      if (ac(ijk).lt.em6.or.beta(ijk).lt.0.0) go to 730
      if (ar(ijk).lt.em6) go to 650
      if (nf(ijk).eq.0.and.nf(ipjk).eq.0) go to 620
      if (nf(ijk).eq.0.and.nf(ipjk).ne.1) go to 630
      if (nf(ijk).ne.2.and.nf(ipjk).eq.0) go to 640
  620 u(ijk)=u(ijk)+delt*(1.0*(dp(lijk)-dp(lipjk)))*rdxp(i)
      go to 650
  630 u(ijk)=u(ijk)+delt*dp(lijk)*rdxp(i)
      go to 650
  640 u(ijk)=u(ijk)+delt*(-dp(lipjk))*rdxp(i)
      go to 650
  650 if (abk(ijk).lt.em6) go to 690
      if (nf(ijk).eq.0.and.nf(ijpk).eq.0) go to 660
      if (nf(ijk).eq.0.and.nf(ijpk).ne.3) go to 670
      if (nf(ijk).ne.4.and.nf(ijpk).eq.0) go to 680
  660 v(ijk)=v(ijk)+delt*(1.0*(dp(lijk)-dp(lijpk)))*rdyp(j)*rri(i)
      go to 690
  670 v(ijk)=v(ijk)+delt*dp(lijk)*rdyp(j)*rri(i)
      go to 690
  680 v(ijk)=v(ijk)+delt*(-dp(lijpk))*rdyp(j)*rri(i)
  690 if (at(ijk).lt.em6) go to 730
      if (nf(ijk).eq.0.and.nf(ijkp).eq.0) go to 700
      if (nf(ijk).eq.0.and.nf(ijkp).ne.5) go to 710
      if (nf(ijk).ne.6.and.nf(ijkp).eq.0) go to 720
  700 w(ijk)=w(ijk)+delt*(1.0*(dp(lijk)-dp(lijkp)))*rdzp(k)
      go to 730
  710 w(ijk)=w(ijk)+delt*dp(lijk)*rdzp(k)
      go to 730
  720 w(ijk)=w(ijk)+delt*(-dp(lijkp))*rdzp(k)
  730 continue
c
c +++ put in pressures for surface cells, isolated cells, voids
c
      do 740 k=1,kmax
      kkk=ii5*(k-1)
      do 740 j=1,jmax
      jjj=kkk+imax*(j-1)
      do 740 i=1,imax
      if ((nf(nq*(jjj+i-1)+1).eq.0).or.(ac(nq*(jjj+i-1)+1).lt.em6).or.
     1 (beta(nq*(jjj+i-1)+1).lt.0.0)) go to 740
      p(nq*(jjj+i-1)+1)=pn(nq*(jjj+i-1)+1)+dp(jjj+i)
  740 continue
c
      return
c
  750 write (59,770) cycle
      write (12,770) cycle
      write (9,770) cycle
      call exit
c
  760 write (59,780) cycle
      write (12,780) cycle
      write (9,780) cycle
      call exit
c
  770 format (1x,29htoo many cr iterations, cycle,i7)
  780 format (1x,44hmesh dimensions and storage not equal, cycle,i7)
      end
*dk pressit
      subroutine pressit
*ca slcom1
      include "param.fi"
      data itmax /1000/
c
c * * pressure iteration
c
   10 if (iorder.eq.2) alpha=1.0
      if (flg.eq.0.0) go to 180
      iter=iter+1
      if (iter.lt.itmax) go to 20
      fnoc=1.0
      nocon=nocon+1
!      call lprt2
      go to 180
   20 flg=0.0
c
c * * compute updated cell pressure and velocities
c
      do 170 k=2,km1
      do 170 j=2,jm1
      do 170 i=2,im1
      call calcijk
      if (beta(ijk).lt.0.0) go to 170
      if (f(ijk).lt.emf) go to 170
      if (nf(ijk).eq.0) go to 100
c
c * * calculate pressure change for surface cells
c
      nff=nf(ijk)
      l=i
      m=j
      n=k
      go to (30,40,50,60,70,80,170), nff
   30 l=i-1
      go to 90
   40 l=i+1
      go to 90
   50 m=j-1
      go to 90
   60 m=j+1
      go to 90
   70 n=k-1
      go to 90
   80 n=k+1
   90 continue
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      nfer=nf(ipjk)
      nfel=nf(imjk)
      nfebk=nf(ijpk)
      nfef=nf(ijmk)
      nfeb=nf(ijkm)
      nfet=nf(ijkp)
      nfe=max0(nfel,nfer,nfebk,nfef,nfeb,nfet,1)
      psurf=pr(nfe)+ps(ijk)
      plmn=p(lmn)
      if (nf(lmn).ne.0.and.beta(ijk).gt.0.0) plmn=psurf
      delp=(1.0-peta(ijk))*plmn+peta(ijk)*psurf-p(ijk)
      go to 110
  100 continue
c
c * * calculate pressure change for full cells
c     and test for convergence
c
      d(ijk)=(rri(i)*(rdx(i)*(u(ijk)*ar(ijk)/rr(i)-u(imjk)*ar(imjk)/rr(i
     1 -1))+rdy(j)*(v(ijk)*abk(ijk)-v(ijmk)*abk(ijmk)))+rdz(k)*(w(ijk)
     2 *at(ijk)-w(ijkm)*at(ijkm)))/ac(ijk)
      dijk=d(ijk)
      if (abs(dijk).ge.epsi) flg=1.0
      delp=-beta(ijk)*dijk*peta(ijk)
c
c * * update cell pressure and velocity
c
  110 dtdelp=delt*delp
      p(ijk)=p(ijk)+delp
      if (beta(ipjk).lt.0.0.or.beta(ipjk).eq.1.0) go to 120
      u(ijk)=u(ijk)+dtdelp*rdxp(i)
  120 if (beta(imjk).lt.0.0.or.beta(imjk).eq.1.0) go to 130
      u(imjk)=u(imjk)-dtdelp*rdxp(i-1)
  130 if (beta(ijpk).lt.0.0.or.beta(ijpk).eq.1.0) go to 140
      v(ijk)=v(ijk)+dtdelp*rdyp(j)*rri(i)
  140 if (beta(ijmk).lt.0.0.or.beta(ijmk).eq.1.0) go to 150
      v(ijmk)=v(ijmk)-dtdelp*rdyp(j-1)*rri(i)
  150 if (beta(ijkp).lt.0.0.or.beta(ijkp).eq.1.0) go to 160
      w(ijk)=w(ijk)+dtdelp*rdzp(k)
  160 if (beta(ijkm).lt.0.0.or.beta(ijkm).eq.1.0) go to 170
      w(ijkm)=w(ijkm)-dtdelp*rdzp(k-1)
  170 continue
c
c * * impose boundary conditions
c
      call bc
c
      go to 10
  180 continue
c
      return
      end
*dk rcontur
      subroutine rcontur
*ca slcom1
c
c * * contour set-up
c
c * * this set up is for cell centered quantities
c
      include "param.fi"
      namelist /contur/ nnx,nny,nnz,nzx,nzy,nzz,nc,zmn,zmx,dlz,zcq,dmpx
     1 ,dmpy,dmpz,igrd,ititle1,ntitle1,ititle2,ntitle2,xlable,ylable
     2 ,zlable,nxlbl,nylbl,nzlbl
c
c * * default values for contur
c
      nnx=-ibar
      nny=-jbar
      nnz=-kbar
      nzx=50
      nzy=50
      nzz=50
      nc=10
      zmn=-1.0
      zmx=-1.0
      dlz=0.0
      zcq(1)=1.0
      dmpx=x(im1)
      dmpy=y(jm1)
      dmpz=z(km1)
      igrd=0
      ititle1(1)=0   !8hpressure
      ititle1(2)=0   !8hcontours
      ntitle1=16
      ititle2(1)=0   !8hvorticit
      ititle2(2)=0   !8hy contou
      ititle2(3)=0   !8hrs
      ntitle2=18
      xlable=0 !6hx-axis
      ylable=0   !6hy-axis
      zlable=0   !6hz-axis
      nxlbl=6
      nylbl=6
      nzlbl=6
c
c * * read values for contour from namelist input / contur /
c
      read (10,contur)
c
c * * write contour input values onto tape 9 and tape 12.
c
c * * write effective data from "contur" to tapes 9 and 12.
c
      if (lpr.eq.0) go to 30
c
      write (9,40)
      write (12,40)
      write (9,300)
      write (12,300)
c
      kpr=9
      assign 20 to kret
   10 write (kpr,50) nnx
      write (kpr,60) nny
      write (kpr,70) nnz
      write (kpr,80) nzx
      write (kpr,90) nzy
      write (kpr,100) nzz
      write (kpr,110) nc
      write (kpr,120) zmn
      write (kpr,130) zmx
      write (kpr,140) dlz
      write (kpr,150) zcq
      write (kpr,160) dmpx
      write (kpr,170) dmpy
      write (kpr,180) dmpz
      write (kpr,190) igrd
      write (kpr,200) (ititle1(l),l=1,3)
      write (kpr,210) ntitle1
      write (kpr,220) (ititle2(l),l=1,3)
      write (kpr,230) ntitle2
      write (kpr,240) xlable
      write (kpr,250) ylable
      write (kpr,260) zlable
      write (kpr,270) nxlbl
      write (kpr,280) nylbl
      write (kpr,290) nzlbl
      go to kret, (20,30)
   20 assign 30 to kret
      kpr=12
      go to 10
   30 continue
      return
c
   40 format (1h1)
   50 format (10x,10h    nnx=  ,i3)
   60 format (10x,10h    nny=  ,i3)
   70 format (10x,10h    nnz=  ,i3)
   80 format (10x,10h    nzx=  ,i3)
   90 format (10x,10h    nzy=  ,i3)
  100 format (10x,10h    nzz=  ,i3)
  110 format (10x,10h     nc=  ,i3)
  120 format (10x,10h    zmn=  ,1pe12.5)
  130 format (10x,10h    zmx=  ,1pe12.5)
  140 format (10x,10h    dlz=  ,1pe12.5)
  150 format (10x,10h    zcq=  ,1pe12.5)
  160 format (10x,10h   dmpx=  ,1pe12.5)
  170 format (10x,10h   dmpy=  ,1pe12.5)
  180 format (10x,10h   dmpz=  ,1pe12.5)
  190 format (10x,10h   igro=  ,i3)
  200 format (10x,10hititle1=  ,3a8)
  210 format (10x,10hntitle1=  ,i3)
  220 format (10x,10hititle2=  ,3a8)
  230 format (10x,10hntitle2=  ,i3)
  240 format (10x,10h xlable=  ,a8)
  250 format (10x,10h ylable=  ,a8)
  260 format (10x,10h zlable=  ,a8)
  270 format (10x,10h  nxlbl=  ,i3)
  280 format (10x,10h  nylbl=  ,i3)
  290 format (10x,10h  nxlbl=  ,i3)
  300 format (1h0,/,25x,26hcontour plotting specifics,/)
      end
*dk rdtape
      subroutine rdtape
*ca slcom1
c
c * * read input tape
c     
      include "param.fi"
      rewind 7
      jtd=td
      jnsc=1 !locf(zlast)-locf(aa)+1
   10 read (7) (aa(n),n=1,jnsc)
      read (7) (basc(n),n=1,ncr2)
      read (7) (basc1(n),n=1,ncr2)
      read (7) (basc2(n),n=1,ncr2)
      read (7) (basc3(n),n=1,ncr2)
      if (jtd.ne.numtd) go to 10
      close(7)       !call close (7)
!      call second (stim)
      trl=stim
      t1=trl
      t2=t1
!      call dateh (dat)
!      call clock1 (clk)
!      call getjtl (tl)
      tlm=tl-150.0+(1.0-tlimd)*1.0e+10
      write (9,20) jtd
      write (9,30) name,jnm,dat,clk
      write (9,40) iter,t,cycle
      write (12,20) jtd
      write (12,30) name,jnm,dat,clk
      write (12,40) iter,t,cycle
      write (59,20) jtd
      write (59,30) name,jnm,dat,clk
      write (59,40) iter,t,cycle
c
c * * changes to tape  dump inserted here
c
c      changes to allow input parameter modification on restart
c
      call rinput
      call rgrafic
      return
c
   20 format (22h0 restarting from td= ,i3)
   30 format (1h ,18x,10a8,1x,a8,2(1x,a8))
   40 format (6x,6hiter= ,i5,18x,6htime= ,1pe12.5,12x,7hcycle= ,i4)
      end
*dk rgrafic
      subroutine rgrafic
*ca slcom1
c
      include "param.fi"
      namelist /grafic/ xea,yea,zea,xca,yca,zca,nvews,nvplts,ivvew,iv1
     1 ,iv2,jv1,jv2,kv1,kv2,nav,iperv,ncplts,icvew,ic1,ic2,jc1,jc2,kc1
     2 ,kc2,nac,iperc,nsplts,isvew,is1,is2,js1,js2,ks1,ks2,nas,ipers
c*******
c******* default values for namelist /grafic/
c*******
c******* Default values are for one velocity  vector view.  With  full
c******* range on  I and full range on J with k=10 constant.
c*******
      if (td.ge.1) go to 10
      xea(1)=50.0
      yea(1)=50.0
      zea(1)=50.0
      xca(1)=0.05
      yca(1)=0.05
      zca(1)=0.05
      nvews=1
      nvplts=1
      iv1(1)=2
      iv2(1)=99
      jv1(1)=2
      jv2(1)=99
      kv1(1)=10
      kv2(1)=10
      ncplts=0
      nsplts=0
c
c * * read in plot-output specifics
c
   10 read (10,grafic)
c
      ipmx=im1
      jpmx=jm1
      kpmx=km1
      if (nvplts.eq.0) go to 30
      do 20 n=1,nvplts
      if (iv1(n).gt.imax) iv1(n)=ipmx
      if (iv2(n).gt.imax) iv2(n)=ipmx
      if (jv1(n).gt.jmax) jv1(n)=jpmx
      if (jv2(n).gt.jmax) jv2(n)=jpmx
      if (kv1(n).gt.kmax) kv1(n)=kpmx
   20 if (kv2(n).gt.kmax) kv2(n)=kpmx
   30 if (ncplts.eq.0) go to 50
      do 40 n=1,ncplts
      if (ic1(n).gt.imax) ic1(n)=ipmx
      if (ic2(n).gt.imax) ic2(n)=ipmx
      if (jc1(n).gt.jmax) jc1(n)=jpmx
      if (jc2(n).gt.jmax) jc2(n)=jpmx
      if (kc1(n).gt.kmax) kc1(n)=kpmx
   40 if (kc2(n).gt.kmax) kc2(n)=kpmx
   50 if (nsplts.eq.0) go to 70
      do 60 n=1,nsplts
      if (is1(n).gt.imax) is1(n)=ipmx
      if (is2(n).gt.imax) is2(n)=ipmx
      if (js1(n).gt.jmax) js1(n)=jpmx
      if (js2(n).gt.jmax) js2(n)=jpmx
      if (ks1(n).gt.kmax) ks1(n)=kpmx
   60 if (ks2(n).gt.kmax) ks2(n)=kpmx
   70 continue
c
c * * print/write(4020) plotting input data
c
      kpr=9
      assign 90 to kret
   80 write (kpr,110)
      write (kpr,120) (xea(n),yea(n),zea(n),xca(n),yca(n),zca(n),n=1
     1 ,nvews)
      write (kpr,130) (ivvew(n),iv1(n),iv2(n),jv1(n),jv2(n),kv1(n),kv2(n
     1 ),nav(n),iperv(n),n=1,nvplts)
      write (kpr,140) (icvew(n),ic1(n),ic2(n),jc1(n),jc2(n),kc1(n),kc2(n
     1 ),nac(n),iperc(n),n=1,ncplts)
      write (kpr,150) (isvew(n),is1(n),is2(n),js1(n),js2(n),ks1(n),ks2(n
     1 ),nas(n),ipers(n),n=1,nsplts)
      go to kret, (90,100)
   90 assign 100 to kret
      kpr=12
      go to 80
  100 continue
      return
c
c  110 format (1h0,30x,32h  plotting output specifications//)
c  120 format (38h0           eye point-view plane table/(5h  xe=,1pe12.5
c     1 ,5h  ye=,e12.5,5h  ze=,e12.5,5h  xc=,e12.5,5h  yc=,e12.5,5h  zc=
c     2 ,e12.5))
c  130 format (32h0          velocity vector plots/(8h  ivvew=,i3,6h  iv1
c     1=,i3,6h  iv2=,i3,6h  jv1=,i3,6h  jv2=,i3,6h  kv1=,i3,6h  kv2=,i3,6
c     2 h  nav=,i3,8h  iperv=,i3))
c  140 format (24h0         contour plots /(8h  icvew=,i3,6h  ic1=,i3,6h
c     1 ic2=,i3,6h  jc1=,i3,6h  jc2=,i3,6h  kc1=,i3,6h  kc2=,i3,6h  nac=
c     2 ,i3,8h  iperc=,i3))
c  150 format (35h0         perspective surface plots/(8h  isvew=,i3,6h
c     1is1=,i3,6h  is2=,i3,6h  js1=,i3,6h  js2=,i3,6h  ks1=,i3,6h  ks2=
c     2 ,i3,6h  nas=,i3,8h  ipers=,i3))


  110 format (1h0,30x,'plotting output specifications'//)
  120 format ('eye point-view plane table'/(' xe=',1pe12.5
     1 ,'ye=',e12.5,'ze=',e12.5,'xc=',e12.5,'yc=',e12.5,'zc='
     2 ,e12.5))
  130 format ('velocity vector plots'/('ivvew=',i3,'iv1
     1=',i3,'iv2=',i3,'jv1=',i3,'jv2=',i3,'kv1=',i3,'kv2=',i3,
     2 'nav=',i3,'iperv=',i3))
  140 format ('contour plots' /('icvew=',i3,'ic1=',i3,
     1 'ic2=',i3,'jc1=',i3,'jc2=',i3,'kc1=',i3,'kc2=',i3,'nac='
     2 ,i3,'iperc=',i3))
  150 format ('perspective surface plots'/('isvew=',i3,
     1'is1=',i3,'is2=',i3,'js1=',i3,'js2=',i3,'ks1=',i3,'ks2='
     2 ,i3,'nas=',i3,'ipers=',i3))

      end
*dk rinput
      subroutine rinput
*ca slcom1
c
c * * read problem input data
c
      include "param.fi"
      namelist /xput/ alpha,pltdt,prtdt,cyl,delt,epsi,gx,gy,gz,icsurf
     1 ,iequib,jnm,lpr,name,nfcal,nu,omg,tddt,td,tlimd,autot,sigma
     2 ,cangle,idefm,isor,isurft,t,radps,nowall,twfin,ui,velmx,vi,wb,wbk
     3 ,wf,wi,wl,wr,wt,flht,iorder,iclip,izoom,rhof
c
      data ktran /1/
c * * default values for namelist / xput / .
c
      go to (10,20), ktran
   10 continue
      alpha=1.0
      autot=1.0
      pltdt=0.1
      prtdt=0.1
      cyl=0.0
      delt=0.01
      epsi=0.001
      gx=0.0
      gy=-9.81
      gz=0.0
      iclip=0
      icsurf=0
      iequib=0
      iorder=1
      idefm=0
      isor=0
      izoom=0
      jnm=0      !8h run 1
      lpr=2
      name(1)=0  !8hprob. no
      name(2)=0      !8hname
      nfcal=3
      nowall=1
      nu=0.0
      omg=1.0
      flht=1.0
      radps=0.0
      rhof=1.0
      sigma=0.0
      isurft=0
      cangle=0.0
      t=0.0
      tddt=+10
      td=-1.0
      tlimd=1.0
      twfin=10.0
      ui=0.0
      velmx=2.0
      vi=0.0
      wb=2
      wbk=1 !symmetry
      wf=1 !symmetry
      wi=0.0
      wl=2
      wr=2
      wt=2
   20 continue
c
c * * read input data from namelist /xput /.
c
!      read (10,xput)
c
      if (cangle.eq.90.0) cangle=cangle-em6
      cangle=cangle*0.0174532925
      tanca=tan(cangle)
      sang=sin(cangle)
      csang=cos(cangle)
c * * set up film (cgscft) package
c
      go to (30,40), ktran
   30 ktran=2
!      call grphcft
!      call gplot (1hu,14ht3mdt sola-3dv,14)
!      call lib4020
!      call grphlun (12,4hplot)
!      call setflsh
   40 continue
c
c * * write effective data from xput to tapes 9 and 12.
c
      if (lpr.eq.0) go to 70
c
      kpr=9
      assign 60 to kret
   50 write (kpr,80) (name(i),i=1,8)
      write (kpr,90) alpha
      write (kpr,430) autot
      write (kpr,480) cangle
      write (kpr,120) cyl
      write (kpr,130) delt
      write (kpr,140) epsi
      write (kpr,360) flht
      write (kpr,150) gx
      write (kpr,160) gy
      write (kpr,170) gz
      write (kpr,380) iclip
      write (kpr,390) icsurf
      write (kpr,490) idefm
      write (kpr,400) iequib
      write (kpr,500) isor
      write (kpr,370) iorder
      write (kpr,510) isurft
      write (kpr,420) izoom
      write (kpr,180) jnm
      write (kpr,190) lpr
      write (kpr,410) nfcal
      write (kpr,450) nowall
      write (kpr,200) nu
      write (kpr,210) omg
      write (kpr,100) pltdt
      write (kpr,110) prtdt
      write (kpr,440) radps
      write (kpr,460) rhof
      write (kpr,470) sigma
      write (kpr,520) t
      write (kpr,220) tddt
      write (kpr,230) td
      write (kpr,240) tlimd
      write (kpr,250) twfin
      write (kpr,260) ui
      write (kpr,270) velmx
      write (kpr,280) vi
      write (kpr,290) wb
      write (kpr,300) wbk
      write (kpr,310) wf
      write (kpr,320) wi
      write (kpr,330) wl
      write (kpr,340) wr
      write (kpr,350) wt
      go to kret, (60,70)
   60 assign 70 to kret
      kpr=12
      go to 50
   70 continue
      return
c
   80 format (1h1,10a8,//)
   90 format (10x,10h  alpha=  ,1pe12.5)
  100 format (10x,10h  pltdt=  ,1pe12.5)
  110 format (10x,10h  prtdt=  ,1pe12.5)
  120 format (10x,10h    cyl=  ,1pe12.5)
  130 format (10x,10h   delt=  ,1pe12.5)
  140 format (10x,10h   epsi=  ,1pe12.5)
  150 format (10x,10h     gx=  ,1pe12.5)
  160 format (10x,10h     gy=  ,1pe12.5)
  170 format (10x,10h     gz=  ,1pe12.5)
  180 format (10x,10h    jnm=  ,a8)
  190 format (10x,10h    lpr=  ,i2)
  200 format (10x,10h     nu=  ,1pe12.5)
  210 format (10x,10h    omg=  ,1pe12.5)
  220 format (10x,10h   tddt=  ,1pe12.5)
  230 format (10x,10h     td=  ,i4)
  240 format (10x,10h  tlimd=  ,1pe12.5)
  250 format (10x,10h  twfin=  ,1pe12.5)
  260 format (10x,10h     ui=  ,1pe12.5)
  270 format (10x,10h  velmx=  ,1pe12.5)
  280 format (10x,10h     vi=  ,1pe12.5)
  290 format (10x,10h     wb=  ,i2)
  300 format (10x,10h    wbk=  ,i2)
  310 format (10x,10h     wf=  ,i2)
  320 format (10x,10h     wi=  ,1pe12.5)
  330 format (10x,10h     wl=  ,i2)
  340 format (10x,10h     wr=  ,i2)
  350 format (10x,10h     wt=  ,i2)
  360 format (10x,10h   flht=  ,1pe12.5)
  370 format (10x,10h iorder=  ,i2)
  380 format (10x,10h  iclip=  ,i2)
  390 format (10x,10h  icsurf= ,i2)
  400 format (10x,10h  iequib= ,i2)
  410 format (10x,10h  nfcal=  ,i2)
  420 format (10x,10h  izoom=  ,i2)
  430 format (10x,10h  autot=  ,1pe12.5)
  440 format (10x,10h  radps=  ,1pe12.5)
  450 format (10x,10h nowall=  ,i2)
  460 format (10x,10h   rhof=  ,1pe12.5)
  470 format (10x,10h  sigma=  ,1pe12.5)
  480 format (10x,10h cangle=  ,1pe12.5)
  490 format (10x,10h  idefm=  ,i2)
  500 format (10x,10h   isor=  ,i2)
  510 format (10x,10h isurft=  ,i2)
  520 format (10x,10h      t=  ,1pe12.5)
      end
*dk secord
      subroutine secord
*ca slcom1
c
c * * second order option
c
      include "param.fi"
      do 20 j=1,jmax
      do 20 i=1,imax
      kt1=kt1+1
      do 10 k=1,kmax
      call calcijk
      udum=un(ijk)
      un(ijk)=u(ijk)
      u(ijk)=udum
      vdum=vn(ijk)
      vn(ijk)=v(ijk)
      v(ijk)=vdum
      wdum=wn(ijk)
      wn(ijk)=w(ijk)
      w(ijk)=wdum
   10 continue
   20 continue
      alpha=-1.0
      ave=0.5
c
      return
      end
*dk setacc
      subroutine setacc
c
*ca slcom1
      include "param.fi"
      ij=i+(j-1)*imax
      gxa(ij)=gx*cthj(j)+gy*sthj(j)
      gya(ij)=gx*sthjbk(j)+gy*cthjbk(j)
      return
      end
*dk setfs
      subroutine setfs
c
*ca slcom1
      include "param.fi"
      dimension iflg(5), dis(4), xm(5), zm(5)
c
      namelist /fluidgn/ nqbs,qa3,qa2,qa1,qb3,qb2,qb1,qc3,qc2,qc1,qd3
     1 ,qd2,qd1,iqh
      data nqbs /0/
      data qa2, qa1, qb2, qb1, qc2, qc1, iqh /10*0.0,10*0.0,10*0.0,10*0.
     1 0,10*0.0,10*0.0,10*0/
      data qa3, qb3, qc3, qd3, qd2, qd1 /10*0.0,10*0.0,10*0.0,10*0.0,10*
     1 0.0,10*0.0/
c
c * * read in data from namelist / fluidgn /
c
      read (10,fluidgn)
c
c * * write input data for fluid generator onto tapes 9 and 12.
c
      if (nqbs.le.0) go to 410
      if (lpr.eq.0) go to 20
      write (9,450)
      write (12,450)
      write (9,420) nqbs
      write (12,420) nqbs
      do 10 i=1,nqbs
      write (9,430) i,qa2(i),qa1(i),qb2(i),qb1(i),qc2(i),qc1(i),iqh(i)
      write (9,440) i,qa3(i),qb3(i),qc3(i),qd3(i),qd2(i),qd1(i)
      write (12,430) i,qa2(i),qa1(i),qb2(i),qb1(i),qc2(i),qc1(i),iqh(i)
      write (12,440) i,qa3(i),qb3(i),qc3(i),qd3(i),qd2(i),qd1(i)
   10 continue
c
   20 if (icsurf.gt.1) go to 300
c
c +++  this portion for axisymmetric surfaces
c +++ conic fcn=qa2*x*x+qa1*x+qb2*z*z+qb1*z+qc2*x*z+qc1
c +++ inside fcn=negative value
c +++ iqh=1 add fluid inside fcn, iqh=0 subtract fluid inside fcn
c +++
      j=2
      do 250 l=1,nqbs
      do 240 k=2,km1
      do 240 i=2,im1
      call ijkonly
      rdxdz=1.0/(delx(i)*delz(k))
      do 80 m=1,4
      go to (30,40,50,60), m
   30 x1=x(i)
      z1=z(k-1)
      dis(1)=delz(k)
      go to 70
   40 z1=z(k)
      x1=x(i)
      dis(2)=delx(i)
      go to 70
   50 x1=x(i-1)
      z1=z(k)
      dis(3)=delz(k)
      go to 70
   60 z1=z(k-1)
      x1=x(i-1)
      dis(4)=delx(i)
   70 iflg(m)=0
      fconic=qa2(l)*x1*x1+qa1(l)*x1+qb2(l)*z1*z1+qb1(l)*z1+qc2(l)*x1*z1
     1 +qc1(l)
      if (fconic.le.0.0) iflg(m)=1
      xm(m)=x1
      zm(m)=z1
   80 continue
      iflg(5)=iflg(1)
      xm(5)=xm(1)
      zm(5)=zm(1)
      iflgs=0
      do 90 m=1,4
   90 iflgs=iflgs+iflg(m)
      brik=0.0
      btik=0.0
      if (iflgs.eq.0) go to 240
      if (iflgs.lt.4) go to 100
      bik=1.0
      brik=1.0
      btik=1.0
      go to 220
  100 if (iflg(1).eq.1.and.iflg(2).eq.1) brik=1.0
      if (iflg(2).eq.1.and.iflg(3).eq.1) btik=1.0
      do 180 m=1,4
      if (iflg(m).eq.iflg(m+1)) go to 180
      x1=xm(m)
      z1=zm(m)
      x2=xm(m+1)
      z2=zm(m+1)
      if (iflg(m).eq.0) go to 110
      x2=xm(m)
      z2=zm(m)
      x1=xm(m+1)
      z1=zm(m+1)
  110 epsif=0.001*(abs(x2-x1)+abs(z2-z1))
      smn=0.0
      fmn=qa2(l)*x2*x2+qa1(l)*x2+qb2(l)*z2*z2+qb1(l)*z2+qc2(l)*x2*z2+qc1
     1 (l)
      smx=1.0
      fmx=qa2(l)*x1*x1+qa1(l)*x1+qb2(l)*z1*z1+qb1(l)*z1+qc2(l)*x1*z1+qc1
     1 (l)
      s=0.5
  120 xt=s*x1+(1.0-s)*x2
      zt=s*z1+(1.0-s)*z2
      fs=qa2(l)*xt*xt+qa1(l)*xt+qb2(l)*zt*zt+qb1(l)*zt+qc2(l)*xt*zt+qc1
     1 (l)
      if (abs(fs).lt.epsif) go to 150
      if (fs.ge.0.0) go to 130
      fden=abs(fs-fmn)+1.0e-10
      se=s-fs*(s-smn)/fden
      if (se.gt.smx) se=smx
      fmn=fs
      smn=s
      go to 140
  130 fden=abs(fmx-fs)+1.0e-10
      se=s-fs*(smx-s)/fden
      if (se.lt.smn) se=smn
      fmx=fs
      smx=s
  140 si=s-fs*(smx-smn)/(fmx-fmn)
      s=0.5*(se+si)
      go to 120
  150 dis(m)=sqrt((xt-x2)**2+(zt-z2)**2)
      go to (160,170,180,180), m
  160 brik=dis(1)/delz(k)
      go to 180
  170 btik=dis(2)/delx(i)
  180 continue
      m=0
      bik=0.0
  190 continue
      m=m+1
      if (m.eq.5) go to 210
      if (iflg(m).eq.0) go to 190
      mp1=m+1
      if (mp1.eq.5) mp1=1
      mm1=m-1
      if (mm1.eq.0) mm1=4
      bik=bik+dis(m)*dis(mm1)
      if (iflg(mp1).eq.1) go to 200
      dis2=dis(m)
  200 continue
      if (iflg(mm1).eq.1) go to 190
      dis1=dis(mm1)
      go to 190
  210 continue
      if (iflgs.eq.3) bik=bik-dis1*dis2
      bik=0.5*bik*rdxdz
      if (bik.gt.1.0) bik=1.0
  220 continue
      if (iqh(l).eq.1) go to 230
      bik=-bik
  230 f(ijk)=f(ijk)+bik
      if (f(ijk).gt.0.99) f(ijk)=1.0
      if (f(ijk).lt.0.01) f(ijk)=0.0
  240 continue
  250 continue
c
c     note: all jplanes have same values
c
      do 290 j=1,jmax
      if (j.eq.2) go to 290
      do 280 k=1,kmax
      do 270 i=1,imax
      i2k=nq*(ii5*(k-1)+imax+(i-1))+1
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      if (j.eq.1.and.wf.le.2) go to 260
      if (j.eq.jmax.and.wbk.le.2) go to 260
      f(ijk)=f(i2k)
      go to 270
  260 f(ijk)=0.0
  270 continue
  280 continue
  290 continue
c
      go to 410
c
c +++  this portion for non-axisymmetric surfaces
c
c +++  fcn=(x-qa1)**2+(y-qb1)**2+(z-qc1)**2+qd1
c
c      for general surfaces change definition of fcn
c
c       qa2,qa3 are used in this formulation to denote approximate
c       upper and lower limits for i index of mesh containg surface for
c       efficiency purposes.  qb3,qb2 for j index, qc3,qc2 for k index
c
c       qd2 is used to indicate if mesh is initially fluid or void
c       =1.0 fluid, =0.0 void
c
c +++ inside fcn=negative value
c +++ iqh=1 add fluid inside fcn, iqh=0 subtract fluid inside fcn
c +++
c
      data lmax, mmax, nmax /20,40,20/
      data sumvoid, sumliqd, sumvol /0.0,0.0,0.0/
c
  300 rxim1=1.0/x(im1)
      flmax=1.0/float(lmax)
      fmmax=1.0/float(mmax)
      fnmax=1.0/float(nmax)
      do 400 ns=1,nqbs
      ibot=qa2(ns)
      itop=qa3(ns)
      jbot=qb2(ns)
      jtop=qb3(ns)
      kbot=qc2(ns)
      ktop=qc3(ns)
      fiqh=iqh(ns)
c
      do 390 k=1,kmax
      do 390 j=1,jmax
      do 390 i=1,imax
      call ijkonly
      if (ns.gt.1) go to 310
      f(ijk)=1.0
      if (qd2(ns).eq.0.0) f(ijk)=0.0
  310 if (i.eq.1.or.i.eq.imax) go to 390
      if (j.eq.1.or.j.eq.jmax) go to 390
      if (k.eq.1.or.k.eq.kmax) go to 390
      vollc=delx(i)*delz(k)*dely(j)/rri(i)
      if (i.lt.ibot.or.i.gt.itop) go to 380
      if (j.lt.jbot.or.j.gt.jtop) go to 380
      if (k.lt.kbot.or.k.gt.ktop) go to 380
c
      ctest=0.0
      do 320 n=k-1,k
      do 320 m=j-1,j
      do 320 l=i-1,i
      x1=x(l)*cthjbk(m)
      y1=x(l)*sthjbk(m)
      z1=z(n)
      fnc=(x1-qa1(ns))**2+(y1-qb1(ns))**2+(z1-qc1(ns))**2+qd1(ns)
!      ctest=ctest+cvmgp(1.0,0.0,fnc)
  320 continue
      if (ctest.gt.7.5) go to 380
      if (ctest.gt.0.5) go to 330
      f(ijk)=fiqh
      go to 380
c
  330 svolsc=0.0
      deltay=dely(j)*fmmax*rxim1
      deltax=delx(i)*flmax
      deltaz=delz(k)*fnmax
      rsc=x(i-1)-0.5*deltax
      rsd2=x(i-1)
      do 360 l=1,lmax
      rsc=rsc+deltax
      rsd1=rsd2
      rsd2=rsd1+deltax
      ra=sqrt(0.5*(rsd2**2+rsd1**2))
      deltas=deltay*rsc
      volsc=deltax*deltaz*deltas
      ysc=y(j-1)*rsc*rxim1-0.5*deltas
      do 350 m=1,mmax
      ysc=ysc+deltas
      thsc=ysc/rsc
      x1=ra*cos(thsc)
      y1=ra*sin(thsc)
      zsc=z(k-1)-0.5*deltaz
      do 340 n=1,nmax
      zsc=zsc+deltaz
      z1=zsc
c
      fnc=(x1-qa1(ns))**2+(y1-qb1(ns))**2+(z1-qc1(ns))**2+qd1(ns)
c
!      svolsc=svolsc+cvmgp(0.0,volsc,fnc)
c
  340 continue
  350 continue
  360 continue
c
      fvol=svolsc/vollc
      if (iqh(ns).eq.1) f(ijk)=f(ijk)+fvol
      if (iqh(ns).eq.0) f(ijk)=f(ijk)-fvol
      if (f(ijk).le.1.0) go to 370
      f(ijk)=1.0
  370 if (f(ijk).ge.0.0) go to 380
      f(ijk)=0.0
  380 continue
      sumvoid=sumvoid+(1.0-f(ijk))*vollc
      sumliqd=sumliqd+f(ijk)*vollc
      sumvol=sumvol+vollc
c
  390 continue
c
  400 continue
c
  410 return
c
  420 format (2x,6hnqbs= ,i2)
  430 format (2x,2hi=,i1,2x,5hqa2= ,1pe12.5,2x,5hqa1= ,e12.5,2x,5hqb2=
     1 ,e12.5,2x,5hqb1= ,e12.5,2x,5hqc2= ,e12.5,2x,5hqc1= ,e12.5,2x,5hiq
     2h= ,i2)
  440 format (2x,2hi=,i1,2x,5hqa3= ,1pe12.5,2x,5hqb3= ,e12.5,2x,5hqc3=
     1 ,e12.5,2x,5hqd3= ,e12.5,2x,5hqd2= ,e12.5,2x,5hqd1= ,e12.5)
  450 format (///20x,25hconic fluid parameters   )
      end
*dk setup
      subroutine setup
*ca slcom1
c
c * * set constant parameters
c
      include "param.fi"
      pi=3.141592654
      emf=1.0e-6
      emf1=1.0-emf
      fnoc=0.0
      flgc=0.0
      nflgc=0.0
      nocon=0
      flg=0.0
      vchgt=0.0
      numtd=1
      t=0.0
      iter=0
      cycle=0
      if (iorder.eq.2) alpha=1.0
      twtd=tddt
      twplt=pltdt
      twprt=prtdt
      sigma=sigma/rhof
c
c * * special input data
c
c * * set pressure field and other arrays to zero
c
      do 10 k=1,kmax
      do 10 j=1,jmax
      do 10 i=1,imax
      call ijkonly
      p(ijk)=0.0
      pn(ijk)=0.0
      f(ijk)=0.0
      nf(ijk)=0
      peta(ijk)=0.0
      fn(ijk)=0.0
      d(ijk)=0.0
   10 continue
      do 20 l=1,100
   20 pr(l)=0.0
c
c * * compute initial volume fraction function f in cells
c
      if (icsurf.gt.0) go to 70
      if (iequib.eq.0) go to 30
      bond=gz*x(im1)**2/sigma
      call equib (basc(1),basc(1001),1000,bond,cangle,cyl)
      sflht=flht
   30 do 60 i=1,imax
      do 60 j=1,jmax
      ibjk=nq*(ii5+imax*(j-1)+i-1)+1
      ibjkm=ibjk-ii2
      do 50 k=2,kmax
      call ijkonly
      f(ijk)=1.0
      if (iequib.eq.0) go to 40
      ldck=1000.*xi(i)/x(im1)+1.000001
      ldck=min0(1000,ldck)
      ldck=max0(1,ldck)
      flht=sflht+basc(ldck)*x(im1)
   40 if (flht.gt.z(k-1).and.flht.lt.z(k)) f(ijk)=rdz(k)*(flht-z(k-1))
      if (z(k-1).ge.flht) f(ijk)=0.0
   50 continue
      f(ibjkm)=f(ibjk)
   60 continue
      go to 80
c * * generate special f-function (fluid) configuration
   70 call setfs
c
c * * calculate hydrostatic pressure
c
   80 do 90 i=2,im1
      do 90 j=2,jm1
      do 90 kk=2,km1
      k=km1-kk+2
      call ijkajct
      p(ijk)=p(ijkp)-gz*(amin1(f(ijkp),0.5)*delz(k+1)+amax1(f(ijk)-0.5,0
     1 .0)*delz(k))
   90 continue
c
c * * calculate delt to ensure viscous diffusion stability
c
      dtvis=delt*1.0e+10
      do 100 k=2,km1
      do 100 j=2,jm1
      do 100 i=2,im1
      dxsq=delx(i)**2
      dysq=dely(j)**2
      dzsq=delz(k)**2
      rdsq1=dxsq*dysq*dzsq/(dxsq+dysq+dzsq)
      rdsq=rdsq1/(3.0*nu+delt*1.0e-10)
      dtvis=amin1(dtvis,rdsq)
  100 continue
c
c      calculate partial areas
c
      call aset
c
c * * calculate beta(ijk) for mesh
c
      call betacal
c
c    ****** skip print **** remove for debug **********************
c
      go to 140
c
c     print betas ac ar at abk for all planes
      kpr=9
      assign 130 to kret
      if (lpr.ne.3) kpr=12
      if (lpr.ne.3) assign 140 to kret
  110 write (kpr,160)
      write (kpr,170)
      do 120 k=1,kmax
      do 120 j=1,jmax
      do 120 i=1,imax
      ijk=nq*(ii5*(k-1)+imax*(j-1)+(i-1))+1
      write (kpr,180) i,j,k,beta(ijk),ac(ijk),ar(ijk),at(ijk),abk(ijk)
  120 continue
      go to kret, (130,140)
  130 assign 140 to kret
      kpr=12
      go to 110
  140 continue
c * * set initial velocity fields into u,v and w arrays
c
      do 150 j=1,jmax
      do 150 i=1,imax
      do 150 k=1,kmax
      call ijkonly
      call setacc
      u(ijk)=0.0
      v(ijk)=0.0
      w(ijk)=0.0
      if (f(ijk).le.emf) go to 150
      if (beta(ipjk).gt.0.0.and.beta(ipjk).ne.1.0) u(ijk)=ui
      if (beta(ijpk).gt.0.0.and.beta(ijpk).ne.1.0) v(ijk)=vi
      if (beta(ijkp).gt.0.0.and.beta(ijkp).ne.1.0) w(ijk)=wi
  150 continue
c
      return
c
  160 format (1h1)
  170 format (4x,1hi,4x,1hj,4x,1hk,7x,4hbeta,12x,2hac,13x,2har,13x,2hat,
     1 13x,3habk)
  180 format (3i5,1p5e15.5)
      end
*dk sola
      subroutine sola
*ca slcom1
c
c * * numerical solution algorithm
c
      include "param.fi"
      msgcode=0
c * * request that msgcode be set to 1
c * * when a message from a controller arrives.
!      call msgflag (msgcode,2)
!      call second (stim)
!      call fadv (1)
      write (59,190)
      trl=stim
      t1=trl
      t2=t1
!      call getjtl (tl)
      tlm=tl-150.0+(1.0-tlimd)*1.0e+10
      if (td.gt.0) go to 160
c
c * * set initial boundary conditions
c
      call bc
c
      go to 40
c
c * * start time cycle
c
   10 continue
      t=t+delt
      ave=1.0
      cycle=cycle+1
      if (nflgc.ge.200.or.nocon.ge.25) t=1.0e+10
      if (mod(cycle,20).eq.0) write (59,180) t,cycle,iter,delt,epsi
      iter=0
      flg=1.0
      fnoc=0
   20 continue
c
c * * explicitly approximate new time level velocities
c
      call tilde
c
c * * set boundary conditions
c
      call bc
c
c * * second order option
c
      if (iorder.eq.1.or.ave.lt.0.75) go to 30
      call secord
      go to 20
c
c * * iteratively adjust cell pressure and velocity
c
   30 if (isor.eq.1) call pressit
      if (isor.eq.0) call prescr
      call bc
c * * update fluid configuration
c
   40 call vfconv
      if (flgc.gt.0.5) go to 110
c
c * * set boundary conditions
c
      call bc
c
c * * determine pressure inperpolation factor and neighbor
c
      call petacal
      call bc
      if (cycle.lt.1.and.isor.ne.2) call prescr
c
c * * plot velocity vectors and free surface on film
c * * lpr=1,  draw only at twplt
c * * lpr=2,  draw and print data on film at twplt
c * * print data on listing only at twprt
c
      if (cycle.eq.0.or.cycle.eq.1) go to 60
      if (t+1.0e-10.ge.twfin) go to 50
      if (t+1.0e-6.lt.twplt) go to 70
   50 twplt=twplt+pltdt
   60 if (lpr.eq.0) go to 70
      if (lpr.eq.1.or.lpr.eq.2) call draw
      if (lpr.eq.2) call lprt
c
c * * print field variable data on paper
c
   70 continue
      if (cycle.eq.0.or.cycle.eq.1) go to 80
      if (t+1.0e-6.lt.twprt) go to 90
      if (twprt.eq.twplt-pltdt.and.lpr.eq.2) go to 90
      twprt=twprt+prtdt
   80 if (lpr.ne.3) go to 90
!      call lprt2
c
c * * set the advance time arrays into the time-n arrays
c
   90 continue
      do 100 k=1,kmax
      do 100 j=1,jmax
      do 100 i=1,imax
      call ijkonly
      un(ijk)=u(ijk)
      vn(ijk)=v(ijk)
      wn(ijk)=w(ijk)
      fn(ijk)=f(ijk)
      pn(ijk)=p(ijk)
  100 continue
c
c * * adjust delt
c
  110 call deltadj
c
      if (msgcode.ne.1) go to 120
c * * get the message from the controller.
c * * this also resets msgcode to zero.
!      call getmsgr (ibuff,nw,1)
c * * check for the message "halt"
      if (ibuff.eq.4hhalt) go to 140
  120 continue
c
c * * test if problem is finished
c
      if (t+1.0e-10.ge.twfin) go to 140
      if (cycle.le.0) go to 10
c
c * * calculate grind time and print out problem data
c
      told=t2
      if (mod(cycle,20).ne.0) go to 130
!      call second (t2)
      xxx=(t2-told)
      xx=xxx*rijk
      if (lpr.ne.0) write (12,170) t,cycle,iter,xx,delt,epsi,vchgt
     1 ,voftot,xxx
      if (lpr.ne.0) write (9,170) t,cycle,iter,xx,delt,epsi,vchgt,voftot
     1 ,xxx
  130 continue
      if (t2-stim.ge.tlm) go to 140
c
c * * write problem dump tape
c
      if (t.ge.twtd) call wrtape
      go to 10
  140 call wrtape
      if (t+1.0e-10.ge.twfin) go to 150
      if (lpr.eq.1.or.lpr.eq.2) call draw
      if (lpr.eq.2) call lprt
  150 continue
!      call gdone
      return
c
c * * read tape dump
c
  160 call rdtape
      go to 10
c
  170 format (3h t=,1pe12.5,7hcycle=,i4,6hiter=,i4,8hgrinds=,e10.3,
     1 6hdelt=,e11.4,6hepsi=,e10.3,7hvchgt=,e10.3,8hvoftot=,e12.5,/5x,
     2"total time for cycle =" ,e11.4)
  180 format (/,3h t=,1pe10.3,7h cycle=,i5,6h iter=,i5,6h delt=,e10.3,6h
     1 epsi=,e10.3)
  190 format (/,25hentering subroutine sola:,/)
      end
*dk surcntr
      subroutine surcntr (i1,i2,j1,j2,k1,k2,mpln,axp,ayp)
*ca slcom1
      include "param.fi"
c
c * * draw free surface on selected plane velocity vector plot
c
      fpl=0.5
c
c * * select indices and coordinates for surface plane to be plotted
c
      go to (10,20,30), mpln
   10 l1=j1
      l2=j2
      m1=i1
      n1=k1
      n2=k2
      go to 40
   20 l1=i1
      l2=i2
      m1=j1
      n1=k1
      n2=k2
      go to 40
   30 l1=i1
      l2=i2
      m1=k1
      n1=j1
      n2=j2
   40 continue
c
c * * multipliers of 1 or 0 to determine plane plotted
c
      axyp=axp*ayp
      axp1=1.0-axp
      ayp1=1.0-ayp
      m=m1
c
c * * draw free surface
c
      do 140 l=l1,l2
      do 140 n=n1,n2
      xlr=1.0
      xll=1.0
      go to (50,60,70), mpln
   50 i=m
      j=l
      k=n
      xlcp=yj(j+1)
      xlc=yj(j)
      xlcm=yj(j-1)
      yncp=zk(k+1)
      ync=zk(k)
      yncm=zk(k-1)
      go to 80
   60 j=m
      i=l
      k=n
      xlcp=xi(i+1)
      xlc=xi(i)
      xlcm=xi(i-1)
      yncp=zk(k+1)
      ync=zk(k)
      yncm=zk(k-1)
      go to 80
   70 k=m
      i=l
      j=n
      xlcp=xi(i+1)
      xlc=xi(i)
      xlcm=xi(i-1)
      yncp=yj(j+1)/rr(i)
      ync=yj(j)/rr(i)
      yncm=yj(j-1)/rr(i)
      if (cyl.lt.0.5.or.mpln.ne.3) go to 80
      yncp=yj(j+1)/x(im1)
      ync=yj(j)/x(im1)
      yncm=yj(j-1)/x(im1)
      xlr=x(i)
      xll=x(i-1)
      if (xll.eq.0.0) xll=em6
   80 continue
      call calcijk
      lmn=ijk
      lpmn=ipjk*axp+ijpk*axp1
      lmnp=ijpk*ayp1+ijkp*ayp
      lpmnp=ijpkp*axp1+ipjkp*axyp+ipjpk*ayp1
      lmnm=ijmk*ayp1+ijkm*ayp
      lmmn=imjk*axp+ijmk*axp1
      lpmnm=ijpkm*axp1+ipjkm*axyp+ipjmk*ayp1
      lmmnp=ijmkp*axp1+imjkp*axyp+imjpk*ayp1
      if (beta(lmn).lt.0.0) go to 140
      fatr=0.25*(f(lmn)+f(lpmn)+f(lmnp)+f(lpmnp))
      fxtr=0.5*(f(lpmnp)+f(lpmn)-f(lmnp)-f(lmn))/(xlcp-xlc)
      fytr=0.5*(f(lmnp)+f(lpmnp)-f(lmn)-f(lpmn))/(yncp-ync)
      fytr=fytr/xlr
      ftrs=fxtr**2+fytr**2
      if (ftrs.eq.0.0) ftrs=1.0e+10
      xtr=0.5*(xlcp+xlc)+(fpl-fatr)*fxtr/ftrs
      xtr=amax1(xtr,xlc)
      xtr=amin1(xtr,xlcp)
      ytr=0.5*(ync+yncp)+(fpl-fatr)*fytr/(ftrs*xlr)
      ytr=amax1(ytr,ync)
      ytr=amin1(ytr,yncp)
      x2ef=xtr-(axp*xbl+axp1*ybf)
      y2ef=ytr-(ayp*zbb+ayp1*ybf)
      if (f(lmn).gt.0.5.and.f(lpmn).gt.0.5) go to 110
      if (f(lmn).lt.0.5.and.f(lpmn).lt.0.5) go to 110
      fabr=0.25*(f(lmn)+f(lpmn)+f(lmnm)+f(lpmnm))
      fxbr=0.5*(f(lpmn)+f(lpmnm)-f(lmn)-f(lmnm))/(xlcp-xlc)
      fybr=0.5*(f(lmn)+f(lpmn)-f(lmnm)-f(lpmnm))/(ync-yncm)
      fybr=fybr/xlr
      fbrs=fxbr**2+fybr**2
      if (fbrs.eq.0.0) fbrs=1.0e+10
      xbmr=0.5*(xlcp+xlc)+(fpl-fabr)*fxbr/fbrs
      xbmr=amax1(xbmr,xlc)
      xbmr=amin1(xbmr,xlcp)
      ybmr=0.5*(ync+yncm)+(fpl-fabr)*fybr/(fbrs*xlr)
      ybmr=amax1(ybmr,yncm)
      ybmr=amin1(ybmr,ync)
      x1ef=xbmr-(axp*xbl+axp1*ybf)
      y1ef=ybmr-(ayp*zbb+ayp1*ybf)
      go to (100,100,90), mpln
   90 if (cyl.lt.0.5) go to 100
      x2ef=xtr*cos(ytr)-xblc
      y2ef=xtr*sin(ytr)-ybfc
      x1ef=xbmr*cos(ybmr)-xblc
      y1ef=xbmr*sin(ybmr)-ybfc
  100 continue
      ix1=fixl(mpln)+x1ef*xconv(mpln)
      iy1=fiyb+y1ef*yconv(mpln)
      ix2=fixl(mpln)+x2ef*xconv(mpln)
      iy2=fiyb+y2ef*yconv(mpln)
!      call drv (ix1,iy1,ix2,iy2)
  110 continue
      if (f(lmn).gt.0.5.and.f(lmnp).gt.0.5) go to 140
      if (f(lmn).lt.0.5.and.f(lmnp).lt.0.5) go to 140
      fatl=0.25*(f(lmn)+f(lmnp)+f(lmmn)+f(lmmnp))
      fxtl=0.5*(f(lmnp)+f(lmn)-f(lmmnp)-f(lmmn))/(xlc-xlcm)
      fytl=0.5*(f(lmmnp)+f(lmnp)-f(lmmn)-f(lmn))/(yncp-ync)
      fytl=fytl/xll
      ftls=fxtl**2+fytl**2
      if (ftls.eq.0.0) ftls=1.0e+10
      xtl=0.5*(xlcm+xlc)+(fpl-fatl)*fxtl/ftls
      xtl=amax1(xtl,xlcm)
      xtl=amin1(xtl,xlc)
      ytl=0.5*(ync+yncp)+(fpl-fatl)*fytl/(ftls*xll)
      ytl=amax1(ytl,ync)
      ytl=amin1(ytl,yncp)
      x1ef=xtl-(axp*xbl+axp1*ybf)
      y1ef=ytl-(ayp*zbb+ayp1*ybf)
      go to (130,130,120), mpln
  120 if (cyl.lt.0.5) go to 130
      x2ef=xtr*cos(ytr)-xblc
      y2ef=xtr*sin(ytr)-ybfc
      x1ef=xtl*cos(ytl)-xblc
      y1ef=xtl*sin(ytl)-ybfc
  130 continue
      ix1=fixl(mpln)+x1ef*xconv(mpln)
      iy1=fiyb+y1ef*yconv(mpln)
      ix2=fixl(mpln)+x2ef*xconv(mpln)
      iy2=fiyb+y2ef*yconv(mpln)
!      call drv (ix1,iy1,ix2,iy2)
  140 continue
c
      return
      end
*dk surfplt
      subroutine surfplt (i1,i2,j1,j2,k1,k2,na,iper)
*ca slcom1
c
c * * draw free surface in perspective
c
      include "param.fi"
c      call fadv (na)
c      call drfp
c      call lincnt (59)
c      write (12,70)
c      write (12,80) name,jnm,dat,clk
c      write (12,90) iter,t,cycle
c      fpl=0.5
c      i1m=i1-1
c      j1m=j1-1
c      k1m=k1-1
c      do 60 i=i1m,i2
c      do 60 j=j1m,j2
c      do 60 k=k1m,k2
c      call calcijk
c      faim=0.0
c      fajm=0.0
c      fakm=0.0
c      crx=cyl/x(im1)
c      if (beta(ijk).lt.0.0.and.beta(ijpk).lt.0.0.and.beta(ijpkp).lt.0.0
c     1 .and.beta(ijkp).lt.0.0) go to 10
c      if (f(ijk).gt.fpl.and.f(ijpk).gt.fpl.and.f(ijpkp).gt.fpl.and.f
c     1 (ijkp).gt.fpl) go to 10
c      if (f(ijk).lt.fpl.and.f(ijpk).lt.fpl.and.f(ijpkp).lt.fpl.and.f
c     1 (ijkp).lt.fpl) go to 10
c      faim=0.125*(f(ijk)+f(ijpk)+f(ijpkp)+f(ijkp)+f(imjk)+f(imjpk)+f
c     1 (imjpkp)+f(imjkp))
c      faix=0.5*(f(ijk)+f(ijpk)+f(ijpkp)+f(ijkp)-4.0*faim)/(xi(i)-xi(i-1)
c     1 )
c      faiy=0.5*(f(imjk)+f(ijk)+f(ijkp)+f(imjkp)-4.0*faim)*rr(i-1)/(yj(j)
c     1 -yj(j+1))
c      faiz=0.5*(f(ijk)+f(ijpk)+f(imjpk)+f(imjk)-4.0*faim)/(zk(k)-zk(k+1)
c     1 )
c      fais=faix**2+faiy**2+faiz**2
c      if (fais.eq.0.0) fais=1.0e+10
c      xaim=0.5*(xi(i-1)+xi(i))+(fpl-faim)*faix/fais
c      xaim=amax1(xaim,xi(i-1))
c      xaim=amin1(xaim,xi(i))
c      yaim=0.5*(yj(j)+yj(j+1))+(fpl-faim)*faiy*rr(i-1)/fais
c      yaim=amax1(yaim,yj(j))
c      yaim=amin1(yaim,yj(j+1))
c      th1=cyl*rx(im1)*yaim
c      yaim=yaim*(1.0-cyl)+xaim*sin(th1)
c      xaim=xaim*cos(th1)
c      zaim=0.5*(zk(k)+zk(k+1))+(fpl-faim)*faiz/fais
c      zaim=amax1(zaim,zk(k))
c      zaim=amin1(zaim,zk(k+1))
c      th1=crx*yaim
c      yaim=xaim*sin(th1)+yaim*(1.0-cyl)
c      xaim=xaim*cos(th1)
c   10 continue
c      if (beta(ijk).lt.0.0.and.beta(ipjk).lt.0.0.and.beta(ipjkp).lt.0.0
c     1 .and.beta(ijkp).lt.0.0) go to 20
c      if (f(ijk).gt.fpl.and.f(ipjk).gt.fpl.and.f(ipjkp).gt.fpl.and.f
c     1 (ijkp).gt.fpl) go to 20
c      if (f(ijk).lt.fpl.and.f(ipjk).lt.fpl.and.f(ipjkp).lt.fpl.and.f
c     1 (ijkp).lt.fpl) go to 20
c      fajm=0.125*(f(ijk)+f(ipjk)+f(ipjkp)+f(ijkp)+f(ijmk)+f(ipjmk)+f
c     1 (ipjmkp)+f(ijmkp))
c      fajx=0.5*(f(ijk)+f(ijkp)+f(ijmkp)+f(ijmk)-4.0*fajm)/(xi(i)-xi(i+1)
c     1 )
c      fajy=0.5*(f(ijk)+f(ipjk)+f(ipjkp)+f(ijkp)-4.0*fajm)*rr(i)/(yj(j)
c     1 -yj(j-1))
c      fajz=0.5*(f(ijk)+f(ijmk)+f(ipjmk)+f(ipjk)-4.0*fajm)/(zk(k)-zk(k+1)
c     1 )
c      fajs=fajx**2+fajy**2+fajz**2
c      if (fajs.eq.0.0) fajs=1.0e+10
c      xajm=0.5*(xi(i)+xi(i+1))+(fpl-fajm)*fajx/fajs
c      xajm=amax1(xajm,xi(i))
c      xajm=amin1(xajm,xi(i+1))
c      yajm=0.5*(yj(j-1)+yj(j))+(fpl-fajm)*fajy*rr(i)/fajs
c      yajm=amax1(yajm,yj(j-1))
c      yajm=amin1(yajm,yj(j))
c      th1=cyl*rx(im1)*yajm
c      yajm=yajm*(1.0-cyl)+xajm*sin(th1)
c      xajm=xajm*cos(th1)
c      zajm=0.5*(zk(k)+zk(k+1))+(fpl-fajm)*fajz/fajs
c      zajm=amax1(zajm,zk(k))
c      zajm=amin1(zajm,zk(k+1))
c      th1=crx*yajm
c      yajm=xajm*sin(th1)+yajm*(1.0-cyl)
c      xajm=xajm*cos(th1)
c   20 continue
c      if (beta(ijk).lt.0.0.and.beta(ipjk).lt.0.0.and.beta(ipjpk).lt.0.0
c     1 .and.beta(ijpk).lt.0.0) go to 30
c      if (f(ijk).gt.fpl.and.f(ipjk).gt.fpl.and.f(ipjpk).gt.fpl.and.f
c     1 (ijpk).gt.fpl) go to 30
c      if (f(ijk).lt.fpl.and.f(ipjk).lt.fpl.and.f(ipjpk).lt.fpl.and.f
c     1 (ijpk).lt.fpl) go to 30
c      fakm=0.125*(f(ijk)+f(ipjk)+f(ipjpk)+f(ijpk)+f(ijkm)+f(ipjkm)+f
c     1 (ipjpkm)+f(ijpkm))
c      fakx=0.5*(f(ijk)+f(ijpk)+f(ijpkm)+f(ijkm)-4.0*fakm)/(xi(i)-xi(i+1)
c     1 )
c      faky=0.5*(f(ijk)+f(ipjk)+f(ipjkm)+f(ijkm)-4.0*fakm)*rr(i)/(yj(j)
c     1 -yj(j+1))
c      fakz=0.5*(f(ijk)+f(ipjk)+f(ipjpk)+f(ijpk)-4.0*fakm)/(zk(k)-zk(k-1)
c     1 )
c      faks=fakx**2+faky**2+fakz**2
c      if (faks.eq.0.0) faks=1.0e+10
c      xakm=0.5*(xi(i)+xi(i+1))+(fpl-fakm)*fakx/faks
c      xakm=amax1(xakm,xi(i))
c      xakm=amin1(xakm,xi(i+1))
c      yakm=0.5*(yj(j)+yj(j+1))+(fpl-fakm)*faky*rr(i)/faks
c      yakm=amax1(yakm,yj(j))
c      yakm=amin1(yakm,yj(j+1))
c      th1=cyl*rx(im1)*yakm
c      yakm=yakm*(1.0-cyl)+xakm*sin(th1)
c      xakm=xakm*cos(th1)
c      zakm=0.5*(zk(k-1)+zk(k))+(fpl-fakm)*fakz/faks
c      zakm=amax1(zakm,zk(k-1))
c      zakm=amin1(zakm,zk(k))
c      th1=crx*yakm
c      yakm=xakm*sin(th1)+yakm*(1.0-cyl)
c      xakm=xakm*cos(th1)
c   30 continue
c      if (faim.eq.0.0.and.fajm.eq.0.0.and.fakm.eq.0.0) go to 60
c      fac=0.125*(f(ijk)+f(ipjk)+f(ipjpk)+f(ijpk)+f(ijkp)+f(ipjkp)+f
c     1 (ipjpkp)+f(ijpkp))
c      facx=0.5*(f(ijk)+f(ijpk)+f(ijpkp)+f(ijkp)-4.0*fac)/(xi(i)-xi(i+1))
c      facy=0.5*(f(ijk)+f(ipjk)+f(ipjkp)+f(ijkp)-4.0*fac)*rr(i)/(yj(j)-yj
c     1 (j+1))
c      facz=0.5*(f(ijk)+f(ipjk)+f(ipjpk)+f(ijpk)-4.0*fac)/(zk(k)-zk(k+1))
c      facs=facx**2+facy**2+facz**2
c      if (facs.eq.0.0) facs=1.0e+10
c      xac=0.5*(xi(i)+xi(i+1))+(fpl-fac)*facx/facs
c      xac=amax1(xac,xi(i))
c      xac=amin1(xac,xi(i+1))
c      yac=0.5*(yj(j)+yj(j+1))+(fpl-fac)*facy*rr(i)/facs
c      yac=amax1(yac,yj(j))
c      yac=amin1(yac,yj(j+1))
c      th1=cyl*rx(im1)*yac
c      yac=yac*(1.0-cyl)+xac*sin(th1)
c      xac=xac*cos(th1)
c      zac=0.5*(zk(k)+zk(k+1))+(fpl-fac)*facz/facs
c      zac=amax1(zac,zk(k))
c      zac=amin1(zac,zk(k+1))
c      th1=crx*yac
c      yac=xac*sin(th1)+yac*(1.0-cyl)
c      xac=xac*cos(th1)
c      call pcnv (ixic,ietac,xac,yac,zac)
c      if (faim.eq.0.0.or.i.eq.i1m) go to 40
c      call pcnv (ixiim,ietaim,xaim,yaim,zaim)
c      call drvec (ixic,ietac,ixiim,ietaim)
c   40 if (fajm.eq.0.0.or.j.eq.j1m) go to 50
c      call pcnv (ixijm,ietajm,xajm,yajm,zajm)
c      call drvec (ixic,ietac,ixijm,ietajm)
c   50 if (fakm.eq.0.0.or.k.eq.k1m) go to 60
c      call pcnv (ixikm,ietakm,xakm,yakm,zakm)
c      call drvec (ixic,ietac,ixikm,ietakm)
c   60 continue
      return
c
c   70 format (21h perspective surface )
c   80 format (1h ,2x,10a8,1x,a8,2(1x,a8))
c   90 format (6x,5hiter=,i5,18x,6htime= ,1pe12.5,12x,7hcycle= ,i4)
      end
*dk surf10n
      subroutine surf10n
c
*ca slcom1
      include "param.fi"
      dimension hx(3,3), hy(3,3), hz(3,3)
      dimension ane(3,3), anw(3,3), aen(3,3), aes(3,3), cee(3,3), cew(3,
     1 3), cnn(3,3), cns(3,3)
      data nowall /1/
      data sang, csang /0.087155742,0.996194698/
      pslim=25.0/amin1(x(im1),y(jm1),z(km1))
      do 1070 k=2,km1
      do 1070 j=2,jm1
      do 1070 i=2,im1
      call calcijk
      if (beta(ijk).lt.0.0) go to 1070
      nfc=nf(ijk)
      if (nfc.eq.0.or.nfc.gt.6) go to 1070
c
c     check if sufficient fluid present
c
      fsum=0.0
      do 10 n=k-1,k+1
      do 10 m=j-1,j+1
      do 10 l=i-1,i+1
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(l-1))+1
      fsum=fsum+f(lmn)
   10 continue
      if (fsum.lt.1.0) go to 1070
c
c     variables for use in six and seven point derivative formulas
c
c     z coordinate pattern
      dzn=0.5*delz(k+1)
      dz=0.5*delz(k)
      dzs=0.5*delz(k-1)
      dzp=dzn+dz
      dzm=dzs+dz
      zp=dzp/dzm
      zm=1.0/zp
      rdzs=1.0/(dzp+dzm)
      dzd=dzp-dzm
      rpdz=1.0/(dzp*dzm)
      frzp=dz/dzp
      frzm=dz/dzm
      ratzp=dzp*rdzs
      ratzm=dzm*rdzs
      dzpa=dzp+dz
      dzma=dzm+dz
      dzdn=dzn-dz
      dzds=dzs-dz
c
c     mirror pattern for theta coordinate (y)
c
      dthe=0.5*dely(j+1)*rx(im1)
      dth=0.5*dely(j)*rx(im1)
      dthw=0.5*dely(j-1)*rx(im1)
      dthp=dthe+dth
      dthm=dth+dthw
      thp=dthp/dthm
      thm=1.0/thp
      rdths=1.0/(dthp+dthm)
      dthd=dthp-dthm
      rpdth=1.0/(dthp*dthm)
      frthp=dth/dthp
      frthm=dth/dthm
      ratthp=dthp*rdths
      ratthm=dthm*rdths
      dthpa=dthp+dth
      dthma=dthm+dth
      dthde=dthe-dth
      dthdw=dthw-dth
c
c     mirror pattern for radial coordinate
c
      dro=0.5*delx(i+1)
      dr=0.5*delx(i)
      dri=0.5*delx(i-1)
      if (i.ne.2) go to 20
      dro=0.5*delx(4)
      dr=0.5*delx(3)
      dri=0.5*delx(2)
   20 continue
      drp=dro+dr
      drm=dr+dri
      rp=drp/drm
      rm=1.0/rp
      rdrs=1.0/(drp+drm)
      drd=drp-drm
      rpdr=1.0/(drp*drm)
      frrp=dr/drp
      frrm=dr/drm
      ratrp=drp*rdrs
      ratrm=drm*rdrs
      drpa=drp+dr
      drma=drm+dr
      drdo=dro-dr
      drdi=dri-dr
c
      dxr=0.5*(delx(i)+delx(i+1))
      dxl=0.5*(delx(i)+delx(i-1))
      if (i.eq.2) dxl=0.5*(delx(2)+delx(3))
      if (i.eq.2) dxr=0.5*(delx(3)+delx(4))
      dzt=0.5*(delz(k)+delz(k+1))
      dzb=0.5*(delz(k)+delz(k-1))
      dybk=0.5*(dely(j)+dely(j+1))/rri(i)
      dybkp=0.5*(dely(j)+dely(j+1))/rri(i+1)
      dybkm=0.5*(dely(j)+dely(j+1))/rri(i-1)
      if (cyl.eq.1.0.and.i.eq.2) dybkm=dybk
      dyf=0.5*(dely(j)+dely(j-1))/rri(i)
      dyfp=0.5*(dely(j)+dely(j-1))/rri(i+1)
      dyfm=0.5*(dely(j)+dely(j-1))/rri(i-1)
      if (cyl.eq.1.0.and.i.eq.2) dyfm=dyf
      rxden=1.0/(dxl*dxr*(dxl+dxr))
      ryden=1.0/(dyf*dybk*(dyf+dybk))
      rydenp=1.0/(dyfp*dybkp*(dyfp+dybkp))
      rydenm=1.0/(dyfm*dybkm*(dyfm+dybkm))
      rzden=1.0/(dzt*dzb*(dzt+dzb))
c
c     assume near cylindrical surface - collapse in radial direction
c
c      iscr, iscv, isch are flags set during collapsing down loops
c      if equal to 0 then subsequent slopes may be accurate
c      if equal to 1 then slopes cannot be correct : at least one row
c      or column of the collapsing arrays does not contain a
c      surface cell
c
      iscr=0
      iscv=0
      isch=0
c
      iscre=0
      iscrf=0
      do 150 kk=1,3
      n=k-2+kk
      do 140 jj=1,3
      m=j-2+jj
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(i-1))+1
      lmnm=lmn-nq
      lmnp=lmn+nq
      if (cyl.lt.0.5) go to 70
      if (m.eq.1) go to 30
      if (m.eq.jmax) go to 50
      go to 70
c
c     m=1
c
   30 mm=2
      if (wf.eq.4) mm=jm1
   40 lmn=nq*(imax*jmax*(n-1)+imax*(mm-1)+(i-1))+1
      go to 60
c
c      m=jmax
c
   50 mm=jm1
      if (wbk.eq.4) mm=2
      go to 40
   60 lmnm=lmn-nq
      lmnp=lmn+nq
   70 continue
      fofm=1.0
      if (ac(lmnm).le.em6.and.f(imjk).lt.emf) fofm=0.0
      fofp=1.0
      if (ac(lmnp).le.em6.and.f(ipjk).lt.emf) fofp=0.0
      if (i.ge.3.and.i.le.im2) go to 100
      if (i.eq.2) go to 80
      fofmm=1.0
      lmnmm=lmnm-nq
      if (ac(lmnmm).lt.em6.and.f(imjk).lt.emf) fofmm=0.0
      delta=2.0*((1.0+ac(lmnm)*(f(lmnm)-1.0))*xi(i-1)*delx(i-1)*fofm+(1.
     1 0+ac(lmn)*(f(lmn)-1.0))*xi(i)*delx(i)+(1.0+ac(lmnmm)*(f(lmnmm)-1.
     2 0))*xi(i-2)*delx(i-2)*fofmm)
      deltb=2.0*(xi(i-2)*delx(i-2)*fofmm+xi(i-1)*delx(i-1)*fofm+xi(i)
     1 *delx(i))
      if (nfc.ne.2) hx(kk,jj)=sqrt(x(i-3)**2+delta)-x(i-3)
      if (nfc.eq.2) hx(kk,jj)=x(i)-sqrt(x(i)**2-delta)
      if (nfc.ne.2) hxcol=sqrt(x(i-3)**2+deltb)-x(i-3)-em6
      if (nfc.eq.2) hxcol=x(i)-sqrt(x(i)**2-deltb)-em6
      if (hx(kk,jj).le.0.0) iscre=iscre+1
      if (hx(kk,jj).ge.hxcol) iscrf=iscrf+1
      go to 90
   80 fofpp=1.0
      lmnpp=lmnp+nq
      if (ac(lmnpp).lt.em6.and.f(ipjk).lt.emf) fofpp=0.0
      delta=2.0*((1.0+ac(lmn)*(f(lmn)-1.0))*xi(i)*delx(i)+(1.0+ac(lmnp)*
     1 (f(lmnp)-1.0))*xi(i+1)*delx(i+1)*fofp+(1.0+ac(lmnpp)*(f(lmnpp)-1.
     2 0))*xi(i+2)*delx(i+2)*fofpp)
      deltb=2.0*(xi(i+2)*delx(i+2)*fofpp+xi(i+1)*delx(i+1)*fofp+xi(i)
     1 *delx(i))
      if (nfc.ne.2) hx(kk,jj)=sqrt(delta)
      if (nfc.eq.2) hx(kk,jj)=x(i+2)-sqrt(x(i+2)**2-delta)
      if (nfc.ne.2) hxcol=sqrt(deltb)-em6
      if (nfc.eq.2) hxcol=x(i+2)-sqrt(x(i+2)**2-deltb)-em6
      if (hx(kk,jj).le.0.0) iscre=iscre+1
      if (hx(kk,jj).ge.hxcol) iscrf=iscrf+1
   90 continue
      go to 140
  100 lmnmm=lmnm-nq
      lmnpp=lmnp+nq
      fofmm=1.0
      fofpp=1.0
      imm=i-3
      ipp=i+2
      if (i.ne.3) go to 110
      lmnmm=lmnm
      fofmm=0.0
      imm=i-2
      if (ac(lmnpp).lt.em6.and.f(ipjk).lt.emf) fofpp=0.0
      go to 130
  110 if (i.ne.im2) go to 120
      lmnpp=lmnp
      fofpp=0.0
      ipp=i+1
      if (ac(lmnmm).lt.em6.and.f(imjk).lt.emf) fofmm=0.0
      go to 130
  120 if (ac(lmnmm).lt.em6.and.f(imjk).lt.emf) fofmm=0.0
      if (ac(lmnpp).lt.em6.and.f(ipjk).lt.emf) fofpp=0.0
  130 delta=2.0*((1.0+ac(lmnm)*(f(lmnm)-1.0))*xi(i-1)*delx(i-1)*fofm+(1.
     1 0+ac(lmn)*(f(lmn)-1.0))*xi(i)*delx(i)+(1.0+ac(lmnp)*(f(lmnp)-1.0)
     2 )*xi(i+1)*delx(i+1)*fofp+(1.0+ac(lmnmm)*(f(lmnmm)-1.0))*xi(i-2)
     3 *delx(i-2)*fofmm+(1.0+ac(lmnpp)*(f(lmnpp)-1.0))*xi(i+2)*delx(i+2)
     4 *fofpp)
      deltb=2.0*(xi(i-2)*delx(i-2)*fofmm+xi(i-1)*delx(i-1)*fofm+xi(i)
     1 *delx(i)+xi(i+1)*delx(i+1)*fofp+xi(i+2)*delx(i+2)*fofpp)
      if (nfc.ne.2) hx(kk,jj)=sqrt(x(imm)**2+delta)-x(imm)
      if (nfc.eq.2) hx(kk,jj)=x(ipp)-sqrt(x(ipp)**2-delta)
      if (nfc.ne.2) hxcol=sqrt(x(imm)**2+deltb)-x(imm)-em6
      if (nfc.eq.2) hxcol=x(ipp)-sqrt(x(ipp)**2-deltb)-em6
      if (hx(kk,jj).le.0.0) iscre=iscre+1
      if (hx(kk,jj).ge.hxcol) iscrf=iscrf+1
  140 continue
  150 continue
c
      if (k.ne.km1) go to 160
      hx(3,1)=hx(2,1)
      hx(3,2)=hx(2,2)
      hx(3,3)=hx(2,3)
      go to 170
  160 if (k.ne.2) go to 170
      hx(1,1)=hx(2,1)
      hx(1,2)=hx(2,2)
      hx(1,3)=hx(2,3)
  170 continue
c
      phxz=rzden*((hx(3,2)-hx(2,2))*dzb**2+(hx(2,2)-hx(1,2))*dzt**2)
      phxy=ryden*((hx(2,3)-hx(2,2))*dyf**2+(hx(2,2)-hx(2,1))*dybk**2)
c
c      assume near vertical surface - collapse in azmuthal direction
c
      if (jbar.eq.1) go to 350
      iscve=0
      iscvf=0
      do 320 kk=1,3
      n=k-2+kk
      do 310 ii=1,3
      l=i-2+ii
      if (i.eq.2) l=1+ii
      lmn=nq*(imax*jmax*(n-1)+imax*(j-1)+(l-1))+1
      lmnm=lmn-ii1
      lmnp=lmn+ii1
      if (cyl.lt.0.5) go to 250
      if (j.eq.2) go to 180
      if (j.eq.jm1) go to 210
      if (i.eq.2) go to 240
      go to 250
c
c     j=2
c
  180 if (i.eq.2) go to 200
  190 ll=l
      m=2
      if (wf.eq.4) m=jm1
      lmnm=nq*(imax*jmax*(n-1)+imax*(m-1)+(ll-1))+1
      go to 250
  200 if (l.gt.1) go to 190
      m=jop(j)
      ll=2
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(ll-1))+1
      lmnp=lmn-ii1
      lmnm=lmn+ii1
      if (wf.eq.4) go to 250
      lmnp=lmn
      lmnm=lmn-ii1
      go to 250
c
c     j=jm1
c
  210 if (i.eq.2) go to 230
  220 ll=l
      m=jm1
      if (wbk.eq.4) m=2
      lmnp=nq*(imax*jmax*(n-1)+imax*(m-1)+(ll-1))+1
      go to 250
  230 if (l.gt.1) go to 220
      m=jop(j)
      ll=2
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(ll-1))+1
      lmnm=lmn+ii1
      lmnp=lmn-ii1
      if (wbk.eq.4) go to 250
      lmnm=lmn
      lmnp=lmn+ii1
      go to 250
c
c      i=2
c
  240 if (l.gt.1) go to 250
      ll=2
      m=jop(j)
      lmn=nq*(imax*jmax*(n-1)+imax*(m-1)+(ll-1))+1
      lmnm=lmn+ii1
      lmnp=lmn-ii1
      if (wf.eq.4) go to 250
      lmnm=lmn-ii1
      lmnp=lmn+ii1
  250 continue
      fofm=1.0
      if (ac(lmnm).le.em6.and.f(ijmk).lt.emf) fofm=0.0
      if (j.eq.2.and.wf.ne.4.and.f(lmn).lt.emf) fofm=0.0
      fofp=1.0
      if (ac(lmnp).le.em6.and.f(ijpk).lt.emf) fofp=0.0
      if (j.eq.jm1.and.wbk.ne.4.and.f(lmn).lt.emf) fofp=0.0
      if (j.ge.3.and.j.le.jm2) go to 270
      if (j.eq.2) go to 260
      fofmm=1.0
      lmnmm=lmnm-ii1
      if (ac(lmnmm).lt.em6.and.f(ijmk).lt.emf) fofmm=0.0
      hy(ii,kk)=(1.0+ac(lmnm)*(f(lmnm)-1.0))*dely(j-1)*fofm+(1.0+ac(lmn)
     1 *(f(lmn)-1.0))*dely(j)+(1.0+ac(lmnmm)*(f(lmnmm)-1.0))*dely(j-2)
     2 *fofmm
      hy(ii,kk)=hy(ii,kk)/x(im1)
      hycol=(dely(j-2)*fofmm+dely(j-1)*fofm+dely(j))/x(im1)-em6
      if (hy(ii,kk).le.0.0) iscve=iscve+1
      if (hy(ii,kk).ge.hycol) iscvf=iscvf+1
      go to 310
  260 fofpp=1.0
      lmnpp=lmnp+ii1
      if (ac(lmnpp).lt.em6.and.f(ijpk).lt.emf) fofpp=0.0
      hy(ii,kk)=(1.0+ac(lmn)*(f(lmn)-1.0))*dely(j)+(1.0+ac(lmnp)*(f(lmnp
     1 )-1.0))*dely(j+1)+(1.0+ac(lmnpp)*(f(lmnpp)-1.0))*dely(j+2)
      hy(ii,kk)=hy(ii,kk)/x(im1)
      hycol=(dely(j+2)*fofpp+dely(j+1)*fofp+dely(j))/x(im1)-em6
      if (hy(ii,kk).le.0.0) iscve=iscve+1
      if (hy(ii,kk).ge.hycol) iscvf=iscvf+1
      go to 310
  270 lmnmm=lmnm-ii1
      lmnpp=lmnp+ii1
      fofmm=1.0
      fofpp=1.0
      if (j.ne.3) go to 280
      lmnmm=lmnm
      fofmm=0.0
      if (ac(lmnpp).lt.em6.and.f(ijpk).lt.emf) fofpp=0.0
      go to 300
  280 if (j.ne.jm2) go to 290
      lmnpp=lmnp
      fofpp=0.0
      if (ac(lmnmm).lt.em6.and.f(ijmk).lt.emf) fofmm=0.0
      go to 300
  290 if (ac(lmnmm).lt.em6.and.f(ijmk).lt.emf) fofmm=0.0
      if (ac(lmnpp).lt.em6.and.f(ijpk).lt.emf) fofpp=0.0
  300 hy(ii,kk)=(1.0+ac(lmnm)*(f(lmnm)-1.0))*dely(j-1)*fofm+(1.0+ac(lmn)
     1 *(f(lmn)-1.0))*dely(j)+(1.0+ac(lmnp)*(f(lmnp)-1.0))*dely(j+1)
     2 *fofp+(1.0+ac(lmnmm)*(f(lmnmm)-1.0))*dely(j-2)*fofmm+(1.0+ac
     3 (lmnpp)*(f(lmnpp)-1.0))*dely(j+2)*fofpp
      hy(ii,kk)=hy(ii,kk)/x(im1)
      hycol=(dely(j-2)*fofmm+dely(j-1)*fofm+dely(j)+dely(j+1)*fofp+dely
     1 (j+2)*fofpp)/x(im1)-em6
      if (hy(ii,kk).le.0.0) iscve=iscve+1
      if (hy(ii,kk).ge.hycol) iscvf=iscvf+1
  310 continue
  320 continue
c
      if (k.ne.km1) go to 330
      hy(1,3)=hy(1,2)
      hy(2,3)=hy(2,2)
      hy(3,3)=hy(3,2)
      go to 340
  330 if (k.ne.2) go to 340
      hy(1,1)=hy(1,2)
      hy(2,1)=hy(2,2)
      hy(3,1)=hy(3,2)
  340 continue
c
      phyx=rxden*((hy(3,2)-hy(2,2))*dxl**2+(hy(2,2)-hy(1,2))*dxr**2)
      phyz=rzden*((hy(2,3)-hy(2,2))*dzb**2+(hy(2,2)-hy(2,1))*dzt**2)
      if (i.eq.2) phyx=rxden*((hy(2,2)-hy(3,2))*(dxl**2)+(hy(2,2)-hy(1,2
     1 ))*dxr*(dxr+2.0*dxl))
      if (i.eq.2) phyz=rzden*((hy(1,3)-hy(1,2))*(dzb**2)+(hy(1,2)-hy(1,1
     1 ))*(dzt**2))
      phyx=xi(i)*phyx
      phyz=xi(i)*phyz
      go to 360
  350 phyx=1.0e10
      phyz=1.0e10
  360 continue
c
c      assume near horizontal surface - collapse in axial direction
c
      ische=0
      do 520 ii=1,3
      l=i-2+ii
      if (i.eq.2) l=1+ii
      do 510 jj=1,3
      m=j-2+jj
      lmn=nq*(imax*jmax*(k-1)+imax*(m-1)+(l-1))+1
      lmnm=lmn-ii2
      lmnp=lmn+ii2
      if (cyl.lt.0.5) go to 450
      if (m.eq.1) go to 370
      if (m.eq.jmax) go to 400
      if (l.eq.1) go to 430
      go to 450
c
c     m=1
c
  370 if (l.eq.1) go to 390
      ll=l
      mm=2
      if (wf.eq.4) mm=jm1
  380 lmn=nq*(imax*jmax*(k-1)+imax*(mm-1)+(ll-1))+1
      go to 440
  390 ll=2
      mm=jm2
      if (wf.eq.4) mm=jop(j)+1
      go to 380
c
c     m=jmax
c
  400 if (l.eq.1) go to 420
      ll=l
      mm=jm1
      if (wbk.eq.4) m=2
  410 lmn=nq*(imax*jmax*(k-1)+imax*(mm-1)+(ll-1))+1
      go to 440
  420 ll=2
      mm=3
      if (wbk.eq.4) mm=jop(jm1)-1
      go to 410
c
c      i=2
c
  430 ll=2
      if (jj.eq.1) mm=jop(j+1)
      if (jj.eq.2) mm=jop(j)
      if (jj.eq.3) mm=jop(j-1)
      go to 380
  440 lmnm=lmn-ii2
      lmnp=lmn+ii2
  450 continue
      fofm=1.0
      if (ac(lmnm).le.em6.and.f(ijkm).lt.emf) fofm=0.0
      fofp=1.0
      if (ac(lmnp).le.em6.and.f(ijkp).lt.emf) fofp=0.0
      if (k.ge.3.and.k.le.km2) go to 470
      if (k.eq.2) go to 460
      fofmm=1.0
      lmnmm=lmnm-ii2
      if (ac(lmnmm).lt.em6.and.f(ijkm).lt.emf) fofmm=0.0
      hz(ii,jj)=(1.0+ac(lmnm)*(f(lmnm)-1.0))*delz(k-1)*fofm+(1.0+ac(lmn)
     1 *(f(lmn)-1.0))*delz(k)+(1.0+ac(lmnmm)*(f(lmnmm)-1.0))*delz(k-2)
     2 *fofmm
      hzcol=delz(k-2)*fofmm+delz(k-1)*fofm+delz(k)-em6
      if (hz(ii,jj).le.0.0) ische=ische+1
      if (hz(ii,jj).ge.hzcol) ischf=ischf+1
      go to 510
  460 fofpp=1.0
      lmnpp=lmnp+ii2
      if (ac(lmnpp).lt.em6.and.f(ijkp).lt.emf) fofpp=0.0
      hz(ii,jj)=(1.0+ac(lmn)*(f(lmn)-1.0))*delz(k)+(1.0+ac(lmnp)*(f(lmnp
     1 )-1.0))*delz(k+1)+(1.0+ac(lmnpp)*(f(lmnpp)-1.0))*delz(k+2)
      hzcol=delz(k+2)*fofpp+delz(k+1)*fofp+delz(k)-em6
      if (hz(ii,jj).le.0.0) ische=ische+1
      if (hz(ii,jj).ge.hzcol) ischf=ischf+1
      go to 510
  470 lmnmm=lmnm-ii2
      lmnpp=lmnp+ii2
      fofmm=1.0
      fofpp=1.0
      if (k.ne.3) go to 480
      lmnmm=lmnm
      fofmm=0.0
      if (ac(lmnpp).lt.em6.and.f(ijkp).lt.emf) fofpp=0.0
      go to 500
  480 if (k.ne.km2) go to 490
      lmnpp=lmnp
      fofpp=0.0
      if (ac(lmnmm).lt.em6.and.f(ijkm).lt.emf) fofmm=0.0
      go to 500
  490 if (ac(lmnmm).lt.em6.and.f(ijkm).lt.emf) fofmm=0.0
      if (ac(lmnpp).lt.em6.and.f(ijkp).lt.emf) fofpp=0.0
  500 hz(ii,jj)=(1.0+ac(lmnm)*(f(lmnm)-1.0))*delz(k-1)*fofm+(1.0+ac(lmn)
     1 *(f(lmn)-1.0))*delz(k)+(1.0+ac(lmnp)*(f(lmnp)-1.0))*delz(k+1)
     2 *fofp+(1.0+ac(lmnmm)*(f(lmnmm)-1.0))*delz(k-2)*fofmm+(1.0+ac
     3 (lmnpp)*(f(lmnpp)-1.0))*delz(k+2)*fofpp
      hzcol=delz(k-2)*fofmm+delz(k-1)*fofm+delz(k)+delz(k+1)*fofp+delz(k
     1 +2)*fofpp-em6
      if (hz(ii,jj).le.0.0) ische=ische+1
      if (hz(ii,jj).ge.hzcol) ischf=ischf+1
  510 continue
  520 continue
c
      phzx=rxden*((hz(3,2)-hz(2,2))*dxl**2+(hz(2,2)-hz(1,2))*dxr**2)
      phzy=ryden*((hz(2,3)-hz(2,2))*dyf**2+(hz(2,2)-hz(2,1))*dybk**2)
      if (i.eq.2) phzx=rxden*((hz(2,2)-hz(3,2))*(dxl**2)+(hz(2,2)-hz(1,2
     1 ))*(dxr+2.0*dxl)*dxr)
      if (i.eq.2) phzy=ryden*((hz(1,3)-hz(1,2))*(dyf**2)+(hz(1,2)-hz(1,1
     1 ))*(dybk**2))
c
      if (iscr+iscv+isch.eq.3) go to 720
      qx=phxy**2+phxz**2
      qy=phyx**2+phyz**2
      qz=phzx**2+phzy**2
      if (qx+qy+qz.gt.1.0e-10) go to 530
      ps(ijk)=0.0
      go to 1070
c
  530 if (qx.lt.qy.and.qx.lt.qz) go to 540
      if (qy.lt.qx.and.qy.lt.qz) go to 600
      if (qz.lt.qx.and.qz.lt.qy) go to 660
c
  540 nfslp=1
      if (f(ipjk).gt.f(imjk)) nfslp=2
      if (nfslp.eq.nf(ijk)) go to 720
      if (nfcal.eq.2) go to 590
      if (nfcal.eq.3) go to 550
      go to 720
c      conflict
  550 if (beta(imjk).lt.0.0.or.beta(ipjk).lt.0.0) go to 590
c
      inb=0
      nfref=0
      do 560 jj=1,3
      do 560 ii=1,3
      n=k-1
      m=j-2+jj
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 560
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  560 continue
      do 570 kk=2,3
      do 570 ii=1,3
      n=k-2+kk
      m=j-1
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 570
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  570 continue
c
      if (nf(imjk).lt.1.or.nf(imjk).gt.6) go to 580
      nfnb=nf(imjk)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  580 continue
c
      if (nfref.eq.0) go to 720
c
      if (inb.eq.1) go to 720
c
      if (nfref.ne.nfslp) go to 720
c
  590 nf(ijk)=nfslp
      go to 720
c
  600 nfslp=3
      if (f(ijpk).gt.f(ijmk)) nfslp=4
      if (nfslp.eq.nf(ijk)) go to 720
      if (nfcal.eq.2) go to 650
      if (nfcal.eq.3) go to 610
      go to 720
c      conflict
  610 if (beta(ijmk).lt.0.0.or.beta(ijpk).lt.0.0) go to 650
c
      inb=0
      nfref=0
      do 620 jj=1,3
      do 620 ii=1,3
      n=k-1
      m=j-2+jj
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 620
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  620 continue
      do 630 kk=2,3
      do 630 ii=1,3
      n=k-2+kk
      m=j-1
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 630
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  630 continue
c
      if (nf(imjk).lt.1.or.nf(imjk).gt.6) go to 640
      nfnb=nf(imjk)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  640 continue
c
      if (nfref.eq.0) go to 720
c
      if (inb.eq.1) go to 720
c
      if (nfref.ne.nfslp) go to 720
c
  650 nf(ijk)=nfslp
      go to 720
c
  660 nfslp=5
      if (f(ijkp).gt.f(ijkm)) nfslp=6
      if (nfslp.eq.nf(ijk)) go to 720
      if (nfcal.eq.2) go to 710
      if (nfcal.eq.3) go to 670
      go to 720
c     conflict
  670 if (beta(ijkm).lt.0.0.or.beta(ijkp).lt.0.0) go to 710
c
      inb=0
      nfref=0
      do 680 jj=1,3
      do 680 ii=1,3
      n=k-1
      m=j-2+jj
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 680
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  680 continue
      do 690 kk=2,3
      do 690 ii=1,3
      n=k-2+kk
      m=j-1
      l=i-2+ii
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (nf(lmn).lt.1.or.nf(lmn).gt.6) go to 690
      nfnb=nf(lmn)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  690 continue
c
      if (nf(imjk).lt.1.or.nf(imjk).gt.6) go to 700
      nfnb=nf(imjk)
      if (nfref.eq.0) nfref=nfnb
      if (nfnb.ne.nfref) inb=1
  700 continue
c
      if (nfref.eq.0) go to 720
c
      if (inb.eq.1) go to 720
c
      if (nfref.ne.nfslp) go to 720
c
  710 nf(ijk)=nfslp
      go to 720
c
  720 continue
c
      nfs(ijk)=nfslp
      nff=nf(ijk)
      cth=1.0
      sth=0.0
      flge=1.0
      flgw=1.0
      flgn=1.0
      flgs=1.0
      go to (730,730,790,790,880,880), nff
  730 continue
      rdew=rri(i)*rdy(j)
      rdns=rdz(k)
      if (cyl.lt.0.5) go to 760
      hxe=(dth*dthma*rdths/dthp)*hx(2,3)+(dthe*dthma*rpdth)*hx(2,2)-
     1 (dthe*dth*rdths/dthm)*hx(2,1)
      hxw=(dth*dthpa*rdths/dthm)*hx(2,1)+(dthw*dthpa*rpdth)*hx(2,2)-
     1 (dthw*dth*rdths/dthp)*hx(2,3)
      th=0.5*dely(j)*rx(im1)
      cth=cos(th)
      sth=sin(th)
      rdew=1.0/(xi(i)*2.0*sth)
      if (nff.eq.2) go to 740
      if (i.gt.3) rsurf=x(i-3)+hx(2,2)
      if (i.gt.3) rsurfe=x(i-3)+hxe
      if (i.gt.3) rsurfw=x(i-3)+hxw
      if (i.eq.3) rsurf=x(i-2)+hx(2,2)
      if (i.eq.3) rsurfe=x(i-2)+hxe
      if (i.eq.3) rsurfw=x(i-2)+hxw
      if (i.eq.2) rsurf=x(i-1)+hx(2,2)
      if (i.eq.2) rsurf=amax1(rsurf,0.1*delx(2))
      if (i.eq.2) rsurfe=x(i-1)+hxe
      if (i.eq.2) rsurfe=amax1(rsurfe,0.1*delx(2))
      if (i.eq.2) rsurfw=x(i-1)+hxw
      if (i.eq.2) rsurfw=amax1(rsurfw,0.1*delx(2))
      go to 750
  740 if (i.lt.im2) rsurf=x(i+2)-hx(2,2)
      if (i.lt.im2) rsurfe=x(i+2)-hxe
      if (i.lt.im2) rsurfw=x(i+2)-hxw
      if (i.eq.im2) rsurf=x(i+1)-hx(2,2)
      if (i.eq.im2) rsurfe=x(i+1)-hxe
      if (i.eq.im2) rsurfw=x(i+1)-hxw
      if (i.eq.im1) rsurf=x(i)-hx(2,2)
      if (i.eq.im1) rsurfe=x(i)-hxe
      if (i.eq.im1) rsurfw=x(i)-hxw
  750 rdew=1.0/(rsurf*2.0*sth)
      if (nff.eq.2) sth=-sth
  760 continue
      afe=abk(ijk)
      if (cyl.eq.1.0.and.afe.lt.em6) afe=1.0
      afw=abk(ijmk)
      if (cyl.eq.1.0.and.afw.lt.em6) afw=1.0
      afn=at(ijk)
      afs=at(ijkm)
      afloor=ar(imjk)
      if (nff.eq.2) afloor=ar(ijk)
      if (nowall.eq.1) afloor=1.0
c  **   use flg- to indicate wall adhesion at cell bottom
      if (afe.gt.em6.and.f(ijpk).lt.emf) flge=afloor
      if (afw.gt.em6.and.f(ijmk).lt.emf) flgw=afloor
      if (afn.gt.em6.and.f(ijkp).lt.emf) flgn=afloor
      if (afs.gt.em6.and.f(ijkm).lt.emf) flgs=afloor
      phxzb=(hx(2,2)-hx(1,2))/dzb
      phxzt=(hx(3,2)-hx(2,2))/dzt
      if (cyl.lt.0.5) go to 990
      dhna=0.5*(phxzt+phxzb)
      fcyl=0.5*delz(k)*dhna*rxi(i)
      if (nff.eq.2) fcyl=-fcyl
      flgn=flgn*(1.0+fcyl)
      flgs=flgs*(1.0-fcyl)
c
c      building four sets of six point coefficients
c
c      coordinates are such that (e,w) goes with (j+1,j-1)
c      coordinates are such that (n,s) goes with (k+1,k-1)
c
c      z derivative on e face (for dhne)
c
      ane(3,3)=(1.0+3.0*frthp*ratzm)/(3.0*dzp)
      ane(2,3)=-(1.0+zp)*ane(3,3)+(frthp/dzm)
      ane(1,3)=zp*ane(3,3)-(frthp/dzm)
      ane(3,2)=-ane(3,3)+zm*rdzs
      ane(2,2)=(1.0+zp)*ane(3,3)-(frthp/dzm)+dzd*rpdz
      ane(1,2)=-zp*ane(3,3)+(frthp/dzm)-zp*rdzs
      ane(3,1)=0.0
      ane(2,1)=0.0
      ane(1,1)=0.0
c
c     z derivative on w face (for dhnw)
c
      anw(1,1)=-(1.0+3.0*frthm*ratzp)/(3.0*dzm)
      anw(2,1)=-(1.0+zm)*anw(1,1)-(frthm/dzp)
      anw(3,1)=zm*anw(1,1)+(frthm/dzp)
      anw(1,2)=-anw(1,1)-zp*rdzs
      anw(2,2)=(1.0+zm)*anw(1,1)+(frthm/dzp)+dzd*rpdz
      anw(3,2)=-zm*anw(1,1)-(frthm/dzp)+zm*rdzs
      anw(1,3)=0.0
      anw(2,3)=0.0
      anw(3,3)=0.0
c
c      theta derivative on n face  (for dhen)
c
      aen(3,3)=(1.0+3.0*frzp*ratthm)/(3.0*dthp)
      aen(3,2)=-(1.0+thp)*aen(3,3)+(frzp/dthm)
      aen(3,1)=thp*aen(3,3)-(frzp/dthm)
      aen(2,3)=-aen(3,3)+thm*rdths
      aen(2,2)=(1.0+thp)*aen(3,3)-(frzp/dthm)+dthd*rpdth
      aen(2,1)=-thp*aen(3,3)+(frzp/dthm)-thp*rdths
      aen(1,3)=0.0
      aen(1,2)=0.0
      aen(1,1)=0.0
c
c      theta derivative on s face (for dhes)
c
      aes(1,1)=-(1.0+3.0*frzm*ratthp)/(3.0*dthm)
      aes(1,2)=-(1.0+thm)*aes(1,1)-(frzm/dthp)
      aes(1,3)=thm*aes(1,1)+(frzm/dthp)
      aes(2,1)=-aes(1,1)-thp*rdths
      aes(2,2)=(1.0+thm)*aes(1,1)+(frzm/dthp)+dthd*rpdth
      aes(2,3)=-thm*aes(1,1)-(frzm/dthp)+thm*rdths
      aes(3,1)=0.0
      aes(3,2)=0.0
      aes(3,3)=0.0
c
c      building four sets of seven point coefficients
c
c      theta derivative on e face (for dhee)
c
      cee(3,3)=rpdth*rdths*ratzm*(dthde*(dthma**2)+(dthma+dth)*(dthe**2)
     1 )/(3.0*dthp)
      cee(2,3)=-(1.0+zp)*cee(3,3)+(dthma+dth)*rdths/dthp
      cee(1,3)=zp*cee(3,3)
      cee(3,2)=-cee(3,3)
      cee(2,2)=(1.0+zp)*cee(3,3)+(dthde-dthm)*rpdth
      cee(1,2)=-zp*cee(3,3)
      cee(3,1)=0.0
      cee(2,1)=-dthde*rdths/dthm
      cee(1,1)=0.0
c
c     theta derivative on w face  (for dhew)
c
      cew(1,1)=-rpdth*rdths*ratzp*(dthdw*(dthpa**2)+(dthpa+dth)*(dthw**2
     1 ))/(3.0*dthm)
      cew(2,1)=-(1.0+zm)*cew(1,1)-(dthpa+dth)*rdths/dthm
      cew(3,1)=zm*cew(1,1)
      cew(1,2)=-cew(1,1)
      cew(2,2)=(1.0+zm)*cew(1,1)-(dthdw-dthp)*rpdth
      cew(3,2)=-zm*cew(1,1)
      cew(1,3)=0.0
      cew(2,3)=dthdw*rdths/dthp
      cew(3,3)=0.0
c
c     z derivative on n face  (for dhnn)
c
      cnn(3,3)=rpdz*rdzs*ratthm*(dzdn*(dzma**2)+(dzma+dz)*(dzn**2))/(3.0
     1 *dzp)
      cnn(3,2)=-(1.0+thp)*cnn(3,3)+(dzma+dz)*rdzs/dzp
      cnn(3,1)=thp*cnn(3,3)
      cnn(2,3)=-cnn(3,3)
      cnn(2,2)=(1.0+thp)*cnn(3,3)+(dzdn-dzm)*rpdz
      cnn(2,1)=-thp*cnn(3,3)
      cnn(1,3)=0.0
      cnn(1,2)=-dzdn*rdzs/dzm
      cnn(1,1)=0.0
c
c      z derivative on s face  (for dhns)
c
      cns(1,1)=-rpdz*rdzs*ratthp*(dzds*(dzpa**2)+(dzpa+dz)*(dzs**2))/(3.
     1 0*dzm)
      cns(1,2)=-(1.0+thm)*cns(1,1)-(dzpa+dz)*rdzs/dzm
      cns(1,3)=thm*cns(1,1)
      cns(2,1)=-cns(1,1)
      cns(2,2)=(1.0+thm)*cns(1,1)-(dzds-dzp)*rpdz
      cns(2,3)=-thm*cns(1,1)
      cns(3,1)=0.0
      cns(3,2)=dzds*rdzs/dzp
      cns(3,3)=0.0
c
c     obtain radial coordinates for converting theta
c     derivatives to gradients
c
      rrsurf=1.0/rsurf
      rrsurfe=1.0/rsurfe
      rrsurfw=1.0/rsurfw
c
c      initialize gradients for use in do loops
c
      dhee=0.0
      dhew=0.0
      dhen=0.0
      dhes=0.0
      dhne=0.0
      dhnw=0.0
      dhnn=0.0
      dhns=0.0
c
c      nested do loops for construction of gradients
c
c      six or seven point formulas for gradients imposed by
c      nine point summations (3 or 2 terms have been set to zero)
c
      do 780 kk=1,3
      do 770 jj=1,3
      dhne=dhne+ane(kk,jj)*hx(kk,jj)
      dhnw=dhnw+anw(kk,jj)*hx(kk,jj)
      dhen=dhen+aen(kk,jj)*hx(kk,jj)*rrsurf
      dhes=dhes+aes(kk,jj)*hx(kk,jj)*rrsurf
      dhee=dhee+cee(kk,jj)*hx(kk,jj)*rrsurfe
      dhew=dhew+cew(kk,jj)*hx(kk,jj)*rrsurfw
      dhnn=dhnn+cnn(kk,jj)*hx(kk,jj)
      dhns=dhns+cns(kk,jj)*hx(kk,jj)
  770 continue
  780 continue
      go to 990
  790 continue
      rdew=rdx(i)
      rdns=rdz(k)
      phinw=0.25*(hy(1,3)+hy(2,3)+hy(1,2)+hy(2,2))
      phine=0.25*(hy(2,3)+hy(3,3)+hy(2,2)+hy(3,2))
      phisw=0.25*(hy(1,2)+hy(2,2)+hy(1,1)+hy(2,1))
      phise=0.25*(hy(2,2)+hy(3,2)+hy(2,1)+hy(3,1))
      phiave=0.5*(phine+phise)
      phiavw=0.5*(phinw+phisw)
      phiavn=0.5*(phinw+phine)
      phiavs=0.5*(phisw+phise)
      if (j.eq.2) phim=(yj(2)+dely(1))/x(im1)
      if (j.ne.2) phim=(yj(j)-y(j-2))/x(im1)
c ******  y=2 special case ********
      dele=phiave-phim
      delw=phiavw-phim
      if (i.eq.2) delw=0.0
      deln=phiavn-phim
      dels=phiavs-phim
      acose=cos(dele)
      asine=sin(dele)
      acosw=cos(delw)
      asinw=sin(delw)
      acosn=cos(deln)
      asinn=sin(deln)
      acoss=cos(dels)
      asins=sin(dels)
      afe=ar(ijk)
      afw=ar(imjk)
      if (cyl.gt.0.5.and.i.eq.2) afw=1.0
      afn=at(ijk)
      afs=at(ijkm)
      afloor=abk(ijmk)
      if (nff.eq.4) afloor=abk(ijk)
      if (nowall.eq.1) afloor=1.0
      if (cyl.eq.1.0.and.afloor.lt.em6) afloor=1.0
      if (afe.gt.em6.and.f(ipjk).lt.emf) flge=afloor
      if (afw.gt.em6.and.f(imjk).lt.emf) flgw=afloor
      if (afn.gt.em6.and.f(ijkp).lt.emf) flgn=afloor
      if (afs.gt.em6.and.f(ijkm).lt.emf) flgs=afloor
c
c      building four sets of six point coefficients
c
c      coordinates are such that (e,w) goes with (i+1,i-1)
c      coordinates are such that (n,s) goes with (k+1,k-1)
c
c      z derivative on e face (for dhne)
c
      ane(3,3)=(1.0+3.0*frrp*ratzm)/(3.0*dzp)
      ane(2,3)=-(1.0+zp)*ane(3,3)+(frrp/dzm)
      ane(1,3)=zp*ane(3,3)-(frrp/dzm)
      ane(3,2)=-ane(3,3)+zm*rdzs
      ane(2,2)=(1.0+zp)*ane(3,3)-(frrp/dzm)+dzd*rpdz
      ane(1,2)=-zp*ane(3,3)+(frrp/dzm)-zp*rdzs
      ane(3,1)=0.0
      ane(2,1)=0.0
      ane(1,1)=0.0
      if (i.ne.2) go to 800
      ane(3,1)=ane(3,2)
      ane(2,1)=ane(2,2)
      ane(1,1)=ane(1,2)
      ane(3,2)=ane(3,3)
      ane(2,2)=ane(2,3)
      ane(1,2)=ane(1,3)
      ane(3,3)=0.0
      ane(2,3)=0.0
      ane(1,3)=0.0
  800 continue
c
c     z derivative on w face (for dhnw)
c
      anw(1,1)=-(1.0+3.0*frrm*ratzp)/(3.0*dzm)
      anw(2,1)=-(1.0+zm)*anw(1,1)-(frrm/dzp)
      anw(3,1)=zm*anw(1,1)+(frrm/dzp)
      anw(1,2)=-anw(1,1)-zp*rdzs
      anw(2,2)=(1.0+zm)*anw(1,1)+(frrm/dzp)+dzd*rpdz
      anw(3,2)=-zm*anw(1,1)-(frrm/dzp)+zm*rdzs
      anw(1,3)=0.0
      anw(2,3)=0.0
      anw(3,3)=0.0
c
c      r derivative on n face  (for dhen)
c
      aen(3,3)=(1.0+3.0*frzp*ratrm)/(3.0*drp)
      if (i.eq.2) aen(3,3)=-(1.0+rm+3.0*frzp*ratrm)/(3.0*drp)
      aen(3,2)=-(1.0+rp)*aen(3,3)+(frzp/drm)
      aen(3,1)=rp*aen(3,3)-(frzp/drm)
      aen(2,3)=-aen(3,3)+rm*rdrs
      aen(2,2)=(1.0+rp)*aen(3,3)-(frzp/drm)+drd*rpdr
      aen(2,1)=-rp*aen(3,3)+(frzp/drm)-rp*rdrs
      aen(1,3)=0.0
      aen(1,2)=0.0
      aen(1,1)=0.0
      if (i.ne.2) go to 810
      aen(2,3)=aen(2,3)-2.0*rm*rdrs
      aen(2,2)=aen(2,2)+(2.0/drp)
      aen(2,1)=aen(2,1)-(2.0*rdrs)
  810 continue
c
c     r derivative on s face (for dhes)
c
      aes(1,1)=-(1.0+3.0*frzm*ratrp)/(3.0*drm)
      aes(1,2)=-(1.0+rm)*aes(1,1)-(frzm/drp)
      aes(1,3)=rm*aes(1,1)+(frzm/drp)
      aes(2,1)=-aes(1,1)-rp*rdrs
      aes(2,2)=(1.0+rm)*aes(1,1)+(frzm/drp)+drd*rpdr
      aes(2,3)=-rm*aes(1,1)-(frzm/drp)+rm*rdrs
      aes(3,1)=0.0
      aes(3,2)=0.0
      aes(3,3)=0.0
      if (i.ne.2) go to 820
      aes(1,3)=(1.0+rm-3.0*frzm*ratrm)/(3.0*drp)
      aes(1,2)=-(1.0+rp)*aes(1,3)+(frzm/drm)
      aes(1,1)=rp*aes(1,3)-(frzm/drm)
      aes(2,3)=-aes(1,3)-(rm*rdrs)
      aes(2,2)=(1.0+rp)*aes(1,3)-(frzm/drm)+(1.0/drp)+(1.0/drm)
      aes(2,1)=-rp*aes(1,3)+(frzm/drm)-(1.0/drm)-(rdrs)
  820 continue
c
c      building four sets of seven point coefficients
c
c      r derivative on e face (for dhee)
c
      cee(3,3)=rpdr*rdrs*ratzm*(drdo*(drma**2)+(drma+dr)*(dro**2))/(3.0
     1 *drp)
      cee(2,3)=-(1.0+zp)*cee(3,3)+(drma+dr)*rdrs/drp
      cee(1,3)=zp*cee(3,3)
      cee(3,2)=-cee(3,3)
      cee(2,2)=(1.0+zp)*cee(3,3)+(drdo-drm)*rpdr
      cee(1,2)=-zp*cee(3,3)
      cee(3,1)=0.0
      cee(2,1)=-drdo*rdrs/drm
      cee(1,1)=0.0
      if (i.ne.2) go to 830
      cee(1,1)=-zp*cee(3,3)
      cee(2,1)=-(1.0+zm)*cee(1,1)-(drpa+dr)*rdrs/drm
      cee(3,1)=zm*cee(1,1)
      cee(1,2)=-cee(1,1)
      cee(2,2)=(1.0+zm)*cee(1,1)-(drdi-drp)*rpdr
      cee(3,2)=-zm*cee(1,1)
      cee(1,3)=0.0
      cee(2,3)=drdi*rdrs/drp
      cee(3,3)=0.0
  830 continue
c
c     r derivative on w face  (for dhew)
c
      cew(1,1)=-rpdr*rdrs*ratzp*(drdi*(drpa**2)+(drpa+dr)*(dri**2))/(3.0
     1 *drm)
      cew(2,1)=-(1.0+zm)*cew(1,1)-(drpa+dr)*rdrs/drm
      cew(3,1)=zm*cew(1,1)
      cew(1,2)=-cew(1,1)
      cew(2,2)=(1.0+zm)*cew(1,1)-(drdi-drp)*rpdr
      cew(3,2)=-zm*cew(1,1)
      cew(1,3)=0.0
      cew(2,3)=drdi*rdrs/drp
      cew(3,3)=0.0
c
c     z derivative on n face  (for dhnn)
c
      cnn(3,3)=rpdz*rdzs*ratrm*(dzdn*(dzma**2)+(dzma+dz)*(dzn**2))/(3.0
     1 *dzp)
      cnn(3,2)=-(1.0+rp)*cnn(3,3)+(dzma+dz)*rdzs/dzp
      cnn(3,1)=rp*cnn(3,3)
      cnn(2,3)=-cnn(3,3)
      cnn(2,2)=(1.0+rp)*cnn(3,3)+(dzdn-dzm)*rpdz
      cnn(2,1)=-rp*cnn(3,3)
      cnn(1,3)=0.0
      cnn(1,2)=-dzdn*rdzs/dzm
      cnn(1,1)=0.0
      if (i.ne.2) go to 840
      cnn(3,2)=-(1.0+rp)*cnn(3,3)
      cnn(3,1)=rp*cnn(3,3)+(1.0/dzp)-dzdn*rdzs/dzp
      cnn(2,3)=-cnn(3,3)
      cnn(2,2)=(1.0+rp)*cnn(3,3)
      cnn(2,1)=-rp*cnn(3,3)+(dzdn-dzm)*rpdz
      cnn(1,2)=0.0
      cnn(1,1)=-dzdn*rdzs/dzm
  840 continue
c
c      z derivative on s face  (for dhns)
c
      cns(1,1)=-rpdz*rdzs*ratrp*(dzds*(dzpa**2)+(dzpa+dz)*(dzs**2))/(3.0
     1 *dzm)
      cns(1,2)=-(1.0+rm)*cns(1,1)-(dzpa+dz)*rdzs/dzm
      cns(1,3)=rm*cns(1,1)
      cns(2,1)=-cns(1,1)
      cns(2,2)=(1.0+rm)*cns(1,1)-(dzds-dzp)*rpdz
      cns(2,3)=-rm*cns(1,1)
      cns(3,1)=0.0
      cns(3,2)=dzds*rdzs/dzp
      cns(3,3)=0.0
      if (i.ne.2) go to 850
      cns(1,3)=rm*cns(1,1)
      cns(1,2)=-(1.0+rp)*cns(1,3)
      cns(1,1)=rp*cns(1,3)-(1.0/dzm)+dzds*rdzs/dzm
      cns(2,3)=-cns(1,3)
      cns(2,2)=(1.0+rp)*cns(1,3)
      cns(2,1)=-rp*cns(1,3)-(dzds-dzp)*rpdz
      cns(3,2)=0.0
      cns(3,1)=dzds*rdzs/dzp
  850 continue
c
c      initialize gradients for use in do loops
c
      dhee=0.0
      dhew=0.0
      dhen=0.0
      dhes=0.0
      dhne=0.0
      dhnw=0.0
      dhnn=0.0
      dhns=0.0
c
c      nested do loops for construction of gradients
c
c      multiply gradients by appropiate radius on each face
c      to obtain variables with dimensions of angles
c
c      six or seven point formulas for gradients imposed by
c      nine point summations (3 or 2 terms have been set to zero)
c
      do 870 kk=1,3
      do 860 ii=1,3
      dhne=dhne+ane(kk,ii)*hy(ii,kk)*x(i)
      dhnw=dhnw+anw(kk,ii)*hy(ii,kk)*x(i-1)
      dhen=dhen+aen(kk,ii)*hy(ii,kk)*xi(i)
      dhes=dhes+aes(kk,ii)*hy(ii,kk)*xi(i)
      dhee=dhee+cee(kk,ii)*hy(ii,kk)*x(i)
      dhew=dhew+cew(kk,ii)*hy(ii,kk)*x(i-1)
      dhnn=dhnn+cnn(kk,ii)*hy(ii,kk)*xi(i)
      dhns=dhns+cns(kk,ii)*hy(ii,kk)*xi(i)
  860 continue
  870 continue
      go to 990
  880 continue
      rdew=rdx(i)
      rdns=rdy(j)*rri(i)
      afe=ar(ijk)
      afw=ar(imjk)
      if (i.eq.2.and.cyl.eq.1.0.and.x(1).eq.0.0) afw=1.0
      afn=abk(ijk)
      if (j.eq.jm1.and.cyl.eq.1.0.and.jc2pi.eq.0) afn=1.0
      afs=abk(ijmk)
      if (j.eq.2.and.cyl.eq.1.0.and.jc2pi.eq.0) afs=1.0
      afloor=at(ijkm)
      if (nff.eq.6) afloor=at(ijk)
      if (nowall.eq.1) afloor=1.0
      if (afe.gt.em6.and.f(ipjk).lt.emf) flge=afloor
      if (afw.gt.em6.and.f(imjk).lt.emf) flgw=afloor
      if (afn.gt.em6.and.f(ijpk).lt.emf) flgn=afloor
      if (afs.gt.em6.and.f(ijmk).lt.emf) flgs=afloor
      if (cyl.gt.0.5) flge=flge*x(i)*rxi(i)
      if (cyl.gt.0.5) flgw=flgw*x(i-1)*rxi(i)
c
c      building four sets of six point coefficients
c
c      coordinates are such that (e,w) goes with (i+1,i-1)
c      coordinates are such that (n,s) goes with (j+1,j-1)
c
c      theta derivative on e face (for dhne)
c
      ane(3,3)=(1.0+3.0*frrp*ratthm)/(3.0*dthp)
      ane(2,3)=-(1.0+thp)*ane(3,3)+(frrp/dthm)
      ane(1,3)=thp*ane(3,3)-(frrp/dthm)
      ane(3,2)=-ane(3,3)+thm*rdths
      ane(2,2)=(1.0+thp)*ane(3,3)-(frrp/dthm)+dthd*rpdth
      ane(1,2)=-thp*ane(3,3)+(frrp/dthm)-thp*rdths
      ane(3,1)=0.0
      ane(2,1)=0.0
      ane(1,1)=0.0
      if (i.ne.2) go to 890
      ane(3,1)=ane(3,2)
      ane(2,1)=ane(2,2)
      ane(1,1)=ane(1,2)
      ane(3,2)=ane(3,3)
      ane(2,2)=ane(2,3)
      ane(1,2)=ane(1,3)
      ane(3,3)=0.0
      ane(2,3)=0.0
      ane(1,3)=0.0
  890 continue
c
c     theta derivative on w face (for dhnw)
c
      anw(1,1)=-(1.0+3.0*frrm*ratthp)/(3.0*dthm)
      anw(2,1)=-(1.0+thm)*anw(1,1)-(frrm/dthp)
      anw(3,1)=thm*anw(1,1)+(frrm/dthp)
      anw(1,2)=-anw(1,1)-thp*rdths
      anw(2,2)=(1.0+thm)*anw(1,1)+(frrm/dthp)+dthd*rpdth
      anw(3,2)=-thm*anw(1,1)-(frrm/dthp)+thm*rdths
      anw(1,3)=0.0
      anw(2,3)=0.0
      anw(3,3)=0.0
c
c      r derivative on n face  (for dhen)
c
      aen(3,3)=(1.0+3.0*frthp*ratrm)/(3.0*drp)
      if (i.eq.2) aen(3,3)=-(1.0+rm+3.0*frthp*ratrm)/(3.0*drp)
      aen(3,2)=-(1.0+rp)*aen(3,3)+(frthp/drm)
      aen(3,1)=rp*aen(3,3)-(frthp/drm)
      aen(2,3)=-aen(3,3)+rm*rdrs
      aen(2,2)=(1.0+rp)*aen(3,3)-(frthp/drm)+drd*rpdr
      aen(2,1)=-rp*aen(3,3)+(frthp/drm)-rp*rdrs
      aen(1,3)=0.0
      aen(1,2)=0.0
      aen(1,1)=0.0
      if (i.ne.2) go to 900
      aen(2,3)=aen(2,3)-2.0*rm*rdrs
      aen(2,2)=aen(2,2)+(2.0/drp)
      aen(2,1)=aen(2,1)-(2.0*rdrs)
  900 continue
c
c     r derivative on s face (for dhes)
c
      aes(1,1)=-(1.0+3.0*frthm*ratrp)/(3.0*drm)
      aes(1,2)=-(1.0+rm)*aes(1,1)-(frthm/drp)
      aes(1,3)=rm*aes(1,1)+(frthm/drp)
      aes(2,1)=-aes(1,1)-rp*rdrs
      aes(2,2)=(1.0+rm)*aes(1,1)+(frthm/drp)+drd*rpdr
      aes(2,3)=-rm*aes(1,1)-(frthm/drp)+rm*rdrs
      aes(3,1)=0.0
      aes(3,2)=0.0
      aes(3,3)=0.0
      if (i.ne.2) go to 910
      aes(1,3)=(1.0+rm-3.0*frthm*ratrm)/(3.0*drp)
      aes(1,2)=-(1.0+rp)*aes(1,3)+(frthm/drm)
      aes(1,1)=rp*aes(1,3)-(frthm/drm)
      aes(2,3)=-aes(1,3)-(rm*rdrs)
      aes(2,2)=(1.0+rp)*aes(1,3)-(frthm/drm)+(1.0/drp)+(1.0/drm)
      aes(2,1)=-rp*aes(1,3)+(frthm/drm)-(1.0/drm)-(rdrs)
  910 continue
c
c      building four sets of seven point coefficients
c
c      r derivative on e face (for dhee)
c
      cee(3,3)=rpdr*rdrs*ratthm*(drdo*(drma**2)+(drma+dr)*(dro**2))/(3.0
     1 *drp)
      cee(2,3)=-(1.0+thp)*cee(3,3)+(drma+dr)*rdrs/drp
      cee(1,3)=thp*cee(3,3)
      cee(3,2)=-cee(3,3)
      cee(2,2)=(1.0+thp)*cee(3,3)+(drdo-drm)*rpdr
      cee(1,2)=-thp*cee(3,3)
      cee(3,1)=0.0
      cee(2,1)=-drdo*rdrs/drm
      cee(1,1)=0.0
      if (i.ne.2) go to 920
      cee(1,1)=-thp*cee(3,3)
      cee(2,1)=-(1.0+thm)*cee(1,1)-(drpa+dr)*rdrs/drm
      cee(3,1)=thm*cee(1,1)
      cee(1,2)=-cee(1,1)
      cee(2,2)=(1.0+thm)*cee(1,1)-(drdi-drp)*rpdr
      cee(3,2)=-thm*cee(1,1)
      cee(1,3)=0.0
      cee(2,3)=drdi*rdrs/drp
      cee(3,3)=0.0
  920 continue
c
c     r derivative on w face  (for dhew)
c
      cew(1,1)=-rpdr*rdrs*ratthp*(drdi*(drpa**2)+(drpa+dr)*(dri**2))/(3.
     1 0*drm)
      cew(2,1)=-(1.0+thm)*cew(1,1)-(drpa+dr)*rdrs/drm
      cew(3,1)=thm*cew(1,1)
      cew(1,2)=-cew(1,1)
      cew(2,2)=(1.0+thm)*cew(1,1)-(drdi-drp)*rpdr
      cew(3,2)=-thm*cew(1,1)
      cew(1,3)=0.0
      cew(2,3)=drdi*rdrs/drp
      cew(3,3)=0.0
c
c     theta derivative on n face  (for dhnn)
c
      cnn(3,3)=rpdth*rdths*ratrm*(dthde*(dthma**2)+(dthma+dth)*(dthe**2)
     1 )/(3.0*dthp)
      cnn(3,2)=-(1.0+rp)*cnn(3,3)+(dthma+dth)*rdths/dthp
      cnn(3,1)=rp*cnn(3,3)
      cnn(2,3)=-cnn(3,3)
      cnn(2,2)=(1.0+rp)*cnn(3,3)+(dthde-dthm)*rpdth
      cnn(2,1)=-rp*cnn(3,3)
      cnn(1,3)=0.0
      cnn(1,2)=-dthde*rdths/dthm
      cnn(1,1)=0.0
      if (i.ne.2) go to 930
      cnn(3,2)=-(1.0+rp)*cnn(3,3)
      cnn(3,1)=rp*cnn(3,3)+(1.0/dthp)-dthde*rdths/dthp
      cnn(2,3)=-cnn(3,3)
      cnn(2,2)=(1.0+rp)*cnn(3,3)
      cnn(2,1)=-rp*cnn(3,3)+(dthde-dthm)*rpdth
      cnn(1,2)=0.0
      cnn(1,1)=-dthde*rdths/dthm
  930 continue
c
c      theta derivative on s face  (for dhns)
c
      cns(1,1)=-rpdth*rdths*ratrp*(dthdw*(dthpa**2)+(dthpa+dth)*(dthw**2
     1 ))/(3.0*dthm)
      cns(1,2)=-(1.0+rm)*cns(1,1)-(dthpa+dth)*rdths/dthm
      cns(1,3)=rm*cns(1,1)
      cns(2,1)=-cns(1,1)
      cns(2,2)=(1.0+rm)*cns(1,1)-(dthdw-dthp)*rpdth
      cns(2,3)=-rm*cns(1,1)
      cns(3,1)=0.0
      cns(3,2)=dthdw*rdths/dthp
      cns(3,3)=0.0
      if (i.ne.2) go to 940
      cns(1,3)=rm*cns(1,1)
      cns(1,2)=-(1.0+rp)*cns(1,3)
      cns(1,1)=rp*cns(1,3)-(1.0/dthm)+dthdw*rdths/dthm
      cns(2,3)=-cns(1,3)
      cns(2,2)=(1.0+rp)*cns(1,3)
      cns(2,1)=-rp*cns(1,3)-(dthdw-dthp)*rpdth
      cns(3,2)=0.0
      cns(3,1)=dthdw*rdths/dthp
  940 continue
c
c      initialize gradients for use in do loops
c
      dhee=0.0
      dhew=0.0
      dhen=0.0
      dhes=0.0
      dhne=0.0
      dhnw=0.0
      dhnn=0.0
      dhns=0.0
c
c      nested do loops for construction of gradients
c
c      six or seven point formulas for gradients imposed by
c      nine point summations (3 or 2 terms have been set to zero)
c
      do 980 jj=1,3
      do 970 ii=1,3
      dhne=dhne+ane(jj,ii)*hz(ii,jj)*rx(i)
      if (i.eq.2) go to 950
      dhnw=dhnw+anw(jj,ii)*hz(ii,jj)*rx(i-1)
      go to 960
  950 continue
      dhnw=1.0
  960 continue
      dhen=dhen+aen(jj,ii)*hz(ii,jj)
      dhes=dhes+aes(jj,ii)*hz(ii,jj)
      dhee=dhee+cee(jj,ii)*hz(ii,jj)
      dhew=dhew+cew(jj,ii)*hz(ii,jj)
      dhnn=dhnn+cnn(jj,ii)*hz(ii,jj)*rxi(i)
      dhns=dhns+cns(jj,ii)*hz(ii,jj)*rxi(i)
  970 continue
  980 continue
  990 continue
      if (nowall.eq.0) go to 1000
      afe=1.0
      afw=1.0
      afn=1.0
      afs=1.0
 1000 continue
      term1=1.0+dhne*dhne
      term2=1.0+dhnw*dhnw
      rhde=sqrt(term1+dhee*dhee)
      rhdw=sqrt(term2+dhew*dhew)
      if (nff.eq.3.or.nff.eq.4) go to 1020
      few=flge*(afe*dhee/rhde+(1.0-afe)*csang)-flgw*(afw*dhew/rhdw-(1.0
     1 -afw)*csang)
c      add cyl terms for nf=1 or 2
      if (sth.ne.0.0) go to 1010
      go to 1030
 1010 continue
      few=few*cth-sth*(flge*(afe*term1/rhde+(1.0-afe)*sang)+flgw*(afw
     1 *term2/rhdw+(1.0-afw)*sang))
      go to 1030
 1020 continue
      t2re=(1.0+dhne**2)/rhde
      t2rw=(1.0+dhnw**2)/rhdw
      few=flge*(afe*(dhee/rhde*acose+t2re*asine)+(1.0-afe)*csang)-flgw*
     1 (afw*(dhew/rhdw*acosw+t2rw*asinw)-(1.0-afw)*csang)
c
 1030 continue
      term1=1.0+dhen*dhen
      term2=1.0+dhes*dhes
      rhdn=sqrt(term1+dhnn*dhnn)
      rhds=sqrt(term2+dhns*dhns)
      if (nff.eq.3.or.nff.eq.4) go to 1040
      fns=flgn*(afn*dhnn/rhdn+(1.0-afn)*csang)-flgs*(afs*dhns/rhds-(1.0
     1 -afs)*csang)
      go to 1050
 1040 continue
      t2rn=-dhen*dhnn/rhdn
      t2rs=-dhes*dhns/rhds
      fns=flgn*(afn*(dhnn/rhdn*acosn+t2rn*asinn)+(1.0-afn)*csang)-flgs*
     1 (afs*(dhns/rhds*acoss+t2rs*asins)-(1.0-afs)*csang)
 1050 continue
c    add wall adhesion at bottom of cell
      if (nowall.eq.1) go to 1060
      few=few-(2.0-flge-flgw)*sang
      fns=fns-(2.0-flgn-flgs)*sang
c    limit curvature by cell size
 1060 afew=abs(few)
      few=rdew*amin1(afew,2.0)*sign(1.0,few)
      afns=abs(fns)
      fns=rdns*amin1(afns,2.0)*sign(1.0,fns)
c
      fewfns=few+fns
      ps(ijk)=-sigma*amin1(abs(fewfns),pslim)*sign(1.0,fewfns)
c
 1070 continue
c
c      fix floor and ceiling for all j and k i=3,4......im1
c
      pslimi=pslim*sigma
      do 1150 k=2,km1
      do 1150 j=2,jm1
      do 1150 i=2,im1
      call ijkajct
      pn(ijk)=0.0
      if (beta(ijk).le.0.0) go to 1150
      if (nf(ijk).lt.1.or.nf(ijk).gt.6) go to 1150
c
c     treat interior cells of i,j,k mesh by general rule,
c     exception for edge cells
c
      psav=0.0
      riav=0.0
      do 1140 mm=1,3
      do 1140 nn=1,3
      do 1140 ll=1,3
      l=i-2+ll
      m=j-2+mm
      n=k-2+nn
      lmn=nq*(ii5*(n-1)+imax*(m-1)+(l-1))+1
      if (ll.eq.2.and.mm.eq.2.and.nn.eq.2) go to 1140
      if (l.ne.1) go to 1080
      lpmn=lmn+nq
      ps(lmn)=ps(lpmn)
 1080 if (l.ne.imax) go to 1090
      lmmn=lmn-nq
      ps(lmn)=ps(lmmn)
 1090 if (m.ne.1) go to 1100
      lmpn=lmn+nq*imax
      ps(lmn)=ps(lmpn)
 1100 if (m.ne.jmax) go to 1110
      lmmn=lmn-nq*imax
      ps(lmn)=ps(lmmn)
 1110 if (n.ne.1) go to 1120
      lmnp=lmn+nq*ii5
      ps(lmn)=ps(lmnp)
 1120 if (n.ne.kmax) go to 1130
      lmnm=lmn-nq*ii5
      ps(lmn)=ps(lmnm)
 1130 continue
      if (abs(ps(lmn)).lt.1.0) go to 1140
      riav=riav+1.0
      psav=psav+ps(lmn)
 1140 continue
      if (riav.eq.0.0) go to 1150
      psavi=psav/riav
      pn(ijk)=ps(ijk)
      if (abs(ps(ijk)).gt.1.5*abs(psavi).or.abs(ps(ijk)).gt.pslimi) pn
     1 (ijk)=amin1(pslimi,1.5*abs(psavi))*sign(1.0,psavi)
      if (abs(ps(ijk)).lt.0.5*abs(psavi)) pn(ijk)=0.5*abs(psavi)*sign(1.
     1 0,psavi)
 1150 continue
c
      do 1160 k=1,kmax
      do 1160 j=1,jmax
      do 1160 i=1,imax
      call ijkonly
      ps(ijk)=pn(ijk)
      if (i.eq.1.or.i.eq.imax.or.j.eq.1.or.j.eq.jmax.or.k.eq.1.or.k.eq
     1 .kmax) ps(ijk)=0.0
 1160 continue
c
c      fix for i=2 case
c
      i=2
      do 1180 k=2,km1
      do 1180 j=2,jm1
      call ijkajct
      if (beta(ijk).le.0.0) go to 1180
      if (nf(ijk).lt.1.or.nf(ijk).gt.6) go to 1180
c
c     treat interior cells of j,k mesh by general rule,
c     exception for edge cells
c
      psav=0.0
      riav=0.0
      do 1170 mm=1,3
      do 1170 nn=1,3
      m=j-2+mm
      n=k-2+nn
      lmn=nq*(ii5*(n-1)+imax*(m-1)+2)+1
      lmpn=lmn+nq*imax
      lmmn=lmn-nq*imax
      lmnp=lmn+nq*ii5
      lmnm=lmn-nq*ii5
      if (m.eq.1) ps(lmn)=ps(lmpn)
      if (n.eq.1) ps(lmn)=ps(lmnp)
      if (m.eq.jmax) ps(lmn)=ps(lmmn)
      if (n.eq.kmax) ps(lmn)=ps(lmnm)
      if (abs(ps(lmn)).lt.1.0) go to 1170
      riav=riav+1.0
      psav=psav+ps(lmn)
 1170 continue
      if (riav.eq.0.0) go to 1180
      psavi=psav/riav
      if (abs(ps(ijk)).gt.1.5*abs(psavi).or.abs(ps(ijk)).gt.pslimi) ps
     1 (ijk)=amin1(pslimi,1.5*abs(psavi))*sign(1.0,psavi)
      if (abs(ps(ijk)).lt.0.5*abs(psavi)) ps(ijk)=0.5*abs(psavi)*sign(1.
     1 0,psavi)
 1180 continue
c
      i=2
      do 1190 k=1,kmax,km1
      do 1190 j=1,jmax,jm1
      call ijkonly
      ps(ijk)=0.0
 1190 continue
c
c  **********
      return
      end
*dk surf10o
*dk tilde
      subroutine tilde
*ca slcom1
c
c * * explicitly approximate new time level velocities
c
      include "param.fi"
      do 50 j=2,jm1
      do 50 i=2,im1
      do 40 k=2,km1
      call calcijk
      udum=0.0
      vdum=0.0
      wdum=0.0
      if (beta(ijk).lt.0.0) go to 30
c
      ij=i+(j-1)*imax
      gx=gxa(ij)
      gy=gya(ij)
c * * u-velocity
c
      if (beta(ipjk).lt.0.0) go to 10
      if (f(ijk)+f(ipjk).lt.emf) go to 10
      sgu=sign(1.0,un(ijk))
      dudr=(un(ipjk)-un(ijk))*rdx(i+1)
      dudl=(un(ijk)-un(imjk))*rdx(i)
      if (ar(ipjk).lt.em6) dudr=0.0
      if (ar(imjk).lt.em6) dudl=0.0
      dxav=delx(i)+delx(i+1)+alpha*sgu*(delx(i+1)-delx(i))
      fux=un(ijk)*(delx(i)*dudr+delx(i+1)*dudl+alpha*sgu*(delx(i+1)*dudl
     1 -delx(i)*dudr))/dxav
      vbdyt=(delx(i)*vn(ipjk)+delx(i+1)*vn(ijk))/(delx(i)+delx(i+1))
      vbdyb=(delx(i)*vn(ipjmk)+delx(i+1)*vn(ijmk))/(delx(i)+delx(i+1))
      dvdy=vbdyt-vbdyb
      if (abk(ijk).lt.em6.or.abk(ipjk).lt.em6.or.abk(ijmk).lt.em6.or.abk
     1 (ipjmk).lt.em6) dvdy=0.0
      vav=0.5*(vbdyt+vbdyb)
      dyf=0.5*(dely(j)+dely(j+1))
      dya=0.5*(dely(j-1)+dely(j))
      dudf=rr(i)*(un(ijpk)-un(ijk))/dyf
      duda=rr(i)*(un(ijk)-un(ijmk))/dya
      if (ar(ijpk).lt.em6) dudf=0.0
      if (ar(ijmk).lt.em6) duda=0.0
      sgv=sign(1.0,vav)
      dyav=dyf+dya+alpha*sgv*(dyf-dya)
      fuy=(dya*dudf+dyf*duda+alpha*sgv*(dyf*duda-dya*dudf))*(vav/dyav)
      wbdyt=(delx(i+1)*wn(ijk)+delx(i)*wn(ipjk))/(delx(i)+delx(i+1))
      wbdyb=(delx(i+1)*wn(ijkm)+delx(i)*wn(ipjkm))/(delx(i)+delx(i+1))
      wav=0.5*(wbdyt+wbdyb)
      dzt=0.5*(delz(k)+delz(k+1))
      dzb=0.5*(delz(k-1)+delz(k))
      dudt=(un(ijkp)-un(ijk))/dzt
      dudb=(un(ijk)-un(ijkm))/dzb
      if (ar(ijkp).lt.em6) dudt=0.0
      if (ar(ijkm).lt.em6) dudb=0.0
      sgw=sign(1.0,wav)
      dzav=dzt+dzb+alpha*sgw*(dzt-dzb)
      fuz=(dzb*dudt+dzt*dudb+alpha*sgw*(dzt*dudb-dzb*dudt))*(wav/dzav)
      dudxsq=2.0*(dudr-dudl)/(delx(i)+delx(i+1))
      dudysq=rdy(j)*(dudf-duda)*rr(i)
      dudzsq=rdz(k)*(dudt-dudb)
      dudcyl=1.0*rx(i)*(delx(i)*dudr+delx(i+1)*dudl)/(delx(i)+delx(i+1))
     1 -2.0*rr(i)*rx(i)*rdy(j)*dvdy-un(ijk)*rx(i)**2
      visx=nu*(dudxsq+dudysq+dudzsq+cyl*dudcyl)
      fuc=cyl*(-0.25*(vbdyt+vbdyb)**2)*rx(i)
      fru=-radps*(vbdyt+vbdyb)-x(i)*radps**2
      udum=un(ijk)+delt*((p(ijk)-p(ipjk))*2.0/(delx(i)+delx(i+1))+gx-fux
     1 -fuy-fuz+visx-fuc-fru)
   10 continue
c
c * * v-velocity
c
      if (beta(ijpk).lt.0.0) go to 20
      if (f(ijk)+f(ijpk).lt.emf) go to 20
      vbdyr=(delx(i+1)*vn(ijk)+delx(i)*vn(ipjk))/(delx(i)+delx(i+1))
      vbdyl=(delx(i)*vn(imjk)+delx(i-1)*vn(ijk))/(delx(i-1)+delx(i))
      ubdyr=(dely(j)*un(ijpk)+dely(j+1)*un(ijk))/(dely(j)+dely(j+1))
      ubdyl=(dely(j)*un(imjpk)+dely(j+1)*un(imjk))/(dely(j)+dely(j+1))
      uav=0.5*(ubdyr+ubdyl)
      dxr=0.5*(delx(i)+delx(i+1))
      dxl=0.5*(delx(i-1)+delx(i))
      dvdr=(vn(ipjk)-vn(ijk))/dxr
      dvdl=(vn(ijk)-vn(imjk))/dxl
      if (abk(ipjk).lt.em6) dvdr=0.0
      if (abk(imjk).lt.em6) dvdl=0.0
      sgu=sign(1.0,uav)
      dxav=dxr+dxl+alpha*sgu*(dxr-dxl)
      fvx=(dxr*dvdl+dxl*dvdr+alpha*sgu*(dxr*dvdl-dxl*dvdr))*(uav/dxav)
      fvc=cyl*(vbdyr*ubdyr+vbdyl*ubdyl)*0.5*rxi(i)
      sgv=sign(1.0,vn(ijk))
      dvdf=(vn(ijpk)-vn(ijk))*rdy(j+1)*rri(i)
      dvda=(vn(ijk)-vn(ijmk))*rdy(j)*rri(i)
      if (abk(ijpk).lt.em6) dvdf=0.0
      if (abk(ijmk).lt.em6) dvda=0.0
      dyav=dely(j)+dely(j+1)+alpha*sgv*(dely(j+1)-dely(j))
      fvy=vn(ijk)*(dely(j)*dvdf+dely(j+1)*dvda+alpha*sgv*(dely(j+1)*dvda
     1 -dely(j)*dvdf))/dyav
      wbdyt=(dely(j+1)*wn(ijk)+dely(j)*wn(ijpk))/(dely(j)+dely(j+1))
      wbdyb=(dely(j+1)*wn(ijkm)+dely(j)*wn(ijpkm))/(dely(j)+dely(j+1))
      wav=0.5*(wbdyt+wbdyb)
      dzt=0.5*(delz(k)+delz(k+1))
      dzb=0.5*(delz(k-1)+delz(k))
      dvdt=(vn(ijkp)-vn(ijk))/dzt
      dvdb=(vn(ijk)-vn(ijkm))/dzb
      if (abk(ijkp).lt.em6) dvdt=0.0
      if (abk(ijkm).lt.em6) dvdb=0.0
      sgw=sign(1.0,wav)
      dzav=dzt+dzb+alpha*sgw*(dzt-dzb)
      fvz=(dzt*dvdb+dzb*dvdt+alpha*sgw*(dzt*dvdb-dzb*dvdt))*(wav/dzav)
      dvdxsq=rdx(i)*(dvdr-dvdl)
      dvdysq=2.0*rri(i)*(dvdf-dvda)/(dely(j)+dely(j+1))
      dvdzsq=rdz(k)*(dvdt-dvdb)
      dudy=un(ijpk)+un(imjpk)-un(ijk)-un(imjk)
      if (ar(ijpk).lt.em6.or.ar(imjpk).lt.em6.or.ar(ijk).lt.em6.or.ar
     1 (imjk).lt.em6) dudy=0.0
      dvdcyl=rxi(i)*0.5*(dvdr+dvdl)+2.0*rxi(i)*rri(i)*dudy/(dely(j)+dely
     1 (j+1))-vn(ijk)*rxi(i)**2
      visy=nu*(dvdxsq+dvdysq+dvdzsq+cyl*dvdcyl)
      frv=radps*(ubdyr+ubdyl)
      vdum=vn(ijk)+delt*((p(ijk)-p(ijpk))*2.0/(dely(j)+dely(j+1))*rri(i)
     1 +gy-fvx-fvy-fvz+visy-fvc-frv)
   20 continue
c
c * * w-velocity
c
      if (beta(ijkp).lt.0.0) go to 30
      if (f(ijk)+f(ijkp).lt.emf) go to 30
      ubdyr=(delz(k+1)*un(ijk)+delz(k)*un(ijkp))/(delz(k)+delz(k+1))
      ubdyl=(delz(k+1)*un(imjk)+delz(k)*un(imjkp))/(delz(k)+delz(k+1))
      uav=0.5*(ubdyr+ubdyl)
      dxr=0.5*(delx(i)+delx(i+1))
      dxl=0.5*(delx(i-1)+delx(i))
      dwdr=(wn(ipjk)-wn(ijk))/dxr
      dwdl=(wn(ijk)-wn(imjk))/dxl
      if (at(ipjk).lt.em6) dwdr=0.0
      if (at(imjk).lt.em6) dwdl=0.0
      sgu=sign(1.0,uav)
      dxav=dxr+dxl+alpha*sgu*(dxr-dxl)
      fwx=(dxr*dwdl+dxl*dwdr+alpha*sgu*(dxr*dwdl-dxl*dwdr))*(uav/dxav)
      vbdyf=(delz(k+1)*vn(ijk)+delz(k)*vn(ijkp))/(delz(k)+delz(k+1))
      vbdya=(delz(k+1)*vn(ijmk)+delz(k)*vn(ijmkp))/(delz(k)+delz(k+1))
      vav=0.5*(vbdyf+vbdya)
      dyf=0.5*(dely(j)+dely(j+1))
      dya=0.5*(dely(j-1)+dely(j))
      dwdf=rri(i)*(wn(ijpk)-wn(ijk))/dyf
      dwda=rri(i)*(wn(ijk)-wn(ijmk))/dya
      if (at(ijpk).lt.em6) dwdf=0.0
      if (at(ijmk).lt.em6) dwda=0.0
      sgv=sign(1.0,vav)
      dyav=dyf+dya+alpha*sgv*(dyf-dya)
      fwy=(dyf*dwda+dya*dwdf+alpha*sgv*(dyf*dwda-dya*dwdf))*(vav/dyav)
      sgw=sign(1.0,wn(ijk))
      dwdt=(wn(ijkp)-wn(ijk))*rdz(k+1)
      dwdb=(wn(ijk)-wn(ijkm))*rdz(k)
      if (at(ijkp).lt.0.0) dwdt=0.0
      if (at(ijkm).lt.em6) dwdb=0.0
      dzav=delz(k)+delz(k+1)+alpha*sgw*(delz(k+1)-delz(k))
      fwz=(delz(k)*dwdt+delz(k+1)*dwdb+alpha*sgw*(delz(k+1)*dwdb-delz(k)
     1 *dwdt))*(wn(ijk)/dzav)
      dwdxsq=rdx(i)*(dwdr-dwdl)
      dwdysq=rri(i)*rdy(j)*(dwdf-dwda)
      dwdzsq=2.0*(dwdt-dwdb)/(delz(k+1)+delz(k))
      dwdcyl=rxi(i)*0.5*(dwdr+dwdl)
      visz=nu*(dwdxsq+dwdysq+dwdzsq+cyl*dwdcyl)
      wdum=wn(ijk)+delt*((p(ijk)-p(ijkp))*2.0/(delz(k)+delz(k+1))+gz-fwx
     1 -fwy-fwz+visz)
   30 continue
c
c * * compute velocities considering second order
c     accurate option parameter ave
c
      u(ijk)=ave*udum+(1.0-ave)*u(ijk)
      v(ijk)=ave*vdum+(1.0-ave)*v(ijk)
      w(ijk)=ave*wdum+(1.0-ave)*w(ijk)
   40 continue
   50 continue
c
c  **********
      return
      end
*dk velv
      subroutine velv (i1,i2,j1,j2,kk1,kk2,na,iper)
*ca slcom1
c
c * * draw velocity vectors in plane or perspective plots
c
c * * determine plane to be plotted
c
      include "param.fi"
      ifm=0
      jfm=0
      kfm=0
      if (i1.ne.i2) go to 10
      ifm=i1
      mpln=1
      icon=j2-j1+1
      jcon=kk2-kk1+1
   10 if (j1.ne.j2) go to 20
      jfm=j1
      mpln=2
      icon=i2-i1+1
      jcon=kk2-kk1+1
   20 if (kk1.ne.kk2) go to 30
      kfm=kk1
      mpln=3
      icon=i2-i1+1
      jcon=j2-j1+1
   30 continue
c
c * * draw plane or perspective plot frame
c
!      call fadv (na)
      if (iper.eq.2) go to 40
      axp=1.0
      if (mpln.eq.1) axp=0.0
      ayp=1.0
      if (mpln.eq.3) ayp=0.0
      if (na.gt.0) call drf (ifm,jfm,kfm)
      go to 50
   40 if (na.gt.0) call drfp
c
c * * set up to draw velocity vectors
c
   50 ij=1
      kntij=1
      vmax=0.0
      vmaxi=0.0
      vmaxj=0.0
      vmaxk=0.0
      do 110 k=2,kmax
      do 110 j=2,jmax
      do 110 i=2,imax
      call calcijk
      if (f(ijk).lt.0.5) go to 110
      if (ac(ijk).lt.0.5) go to 110
      go to (60,70,80), mpln
   60 if (i.ne.i1) go to 110
      if (j.lt.j1) go to 110
      if (j.gt.j2) go to 110
      if (k.lt.kk1) go to 110
      if (k.gt.kk2) go to 110
      go to 90
   70 if (j.ne.j1) go to 110
      if (k.lt.kk1) go to 110
      if (k.gt.kk2) go to 110
      if (i.lt.i1) go to 110
      if (i.gt.i2) go to 110
      go to 90
   80 if (k.ne.kk1) go to 110
      if (j.lt.j1) go to 110
      if (j.gt.j2) go to 110
      if (i.lt.i1) go to 110
      if (i.gt.i2) go to 110
   90 continue
      sth1=0.0
      cth1=1.0
      cylm=1.0
      if (mpln.ne.3.and.iper.eq.1) go to 100
      th1=cyl*rx(im1)*yj(j)
      sth1=sin(th1)
      cth1=cos(th1)
      cylm=1.0-cyl
  100 continue
c
c * * determine maximum vector length
c
      viner=xi(i)*radps
      viner=0.0
      uave=0.5*(u(ijk)+u(imjk))
      vave=0.5*(v(ijk)+v(ijmk))+viner
      uvv(ij)=uave*cth1-vave*sth1
      vvv(ij)=uave*sth1+vave*cth1
      wvv(ij)=0.5*(w(ijk)+w(ijkm))
      xpc(ij)=xi(i)*cth1
      ypc(ij)=yj(j)*cylm+xi(i)*sth1
      zpc(ij)=zk(k)
      vmax=amax1(vmax,abs(uvv(ij)),abs(vvv(ij)),abs(wvv(ij)))
      vmaxi=amax1(vmaxi,abs(vvv(ij)),abs(wvv(ij)))
      vmaxj=amax1(vmaxj,abs(uvv(ij)),abs(wvv(ij)))
      vmaxk=amax1(vmaxk,abs(uvv(ij)),abs(vvv(ij)))
      if (beta(ijk).lt.0.0.or.f(ijk).lt.0.5) uvv(ij)=1.0e+10
      kt1=km1
      if (k.gt.kt1) uvv(ij)=+1.0e+10
      kntij=kntij+1
      ij=ij+nq2
      if (ij.gt.nvec) go to 310
  110 continue
      kntij=kntij-1
c
c * * scale vector length
c
      go to (120,130,140), mpln
  120 drou=velmx*delmn/(vmaxi+1.0e-10)
      go to 150
  130 drou=velmx*delmn/(vmaxj+1.0e-10)
      go to 150
  140 drou=velmx*delmn/(vmaxk+1.0e-10)
  150 continue
      ij=1
      if (mpln.ne.3) go to 160
      xbl=xblc
      ybf=ybfc
  160 continue
      do 210 l=1,kntij
      x1p=xpc(ij)
      y1p=ypc(ij)
      z1p=zpc(ij)
      u1p=uvv(ij)
      v1p=vvv(ij)
      w1p=wvv(ij)
      if (u1p.gt.1.0e+06) go to 210
      go to (170,190), iper
c
c * * plane vector plot
c
  170 continue
      x1ef=axp*(x1p-xbl)+(1.0-axp)*(y1p-ybf)
      y1ef=ayp*(z1p-zbb)+(1.0-ayp)*(y1p-ybf)
      x2ef=axp*(x1p+u1p*drou-xbl)+(1.0-axp)*(y1p+v1p*drou-ybf)
      y2ef=ayp*(z1p+w1p*drou-zbb)+(1.0-ayp)*(y1p+v1p*drou-ybf)
      ix1=fixl(mpln)+x1ef*xconv(mpln)
      iy1=fiyb+y1ef*yconv(mpln)
      ix2=fixl(mpln)+x2ef*xconv(mpln)
      iy2=fiyb+y2ef*yconv(mpln)
!      call drv (ix1,iy1,ix2,iy2)
c    draw arrowheads
      fix1=ix1
      fix2=ix2
      fiy1=iy1
      fiy2=iy2
      x34=fix1+0.85*(fix2-fix1)
      y34=fiy1+0.85*(fiy2-fiy1)
      dx234=fix2-x34
      dy234=fiy2-y34
      d234=sqrt(dx234**2+dy234**2)
      da234=d234*0.57735027
c      tan(15 deg)=0.57735027
      cosb=da234
      sinb=0.0
      if (dx234.eq.0.0) go to 180
      beta234=1.57079633-atan(dy234/dx234)
      cosb=cos(beta234)*da234
      sinb=sin(beta234)*da234
  180 ix3=x34-cosb
      ix4=x34+cosb
      iy3=y34+sinb
      iy4=y34-sinb
!      call drv (ix2,iy2,ix3,iy3)
!      call drv (ix2,iy2,ix4,iy4)
      go to 210
c
c * * perspective vector plot
c
  190 call pcnv (ixi1,ieta1,x1p,y1p,z1p)
      x2p=x1p+u1p*drou
      y2p=y1p+v1p*drou
      z2p=z1p+w1p*drou
      call pcnv (ixi2,ieta2,x2p,y2p,z2p)
      call drvec (ixi1,ieta1,ixi2,ieta2)
c    draw arrowheads
      fix1=ixi1
      fix2=ixi2
      fiy1=ieta1
      fiy2=ieta2
      x34=fix1+0.85*(fix2-fix1)
      y34=fiy1+0.85*(fiy2-fiy1)
      dx234=fix2-x34
      dy234=fiy2-y34
      d234=sqrt(dx234**2+dy234**2)
      da234=d234*0.57735027
c      tan(15 deg)=0.57735027
      cosb=da234
      sinb=0.0
      if (dx234.eq.0.0) go to 200
      beta234=1.57079633-atan(dy234/dx234)
      cosb=cos(beta234)*da234
      sinb=sin(beta234)*da234
  200 ix3=x34-cosb
      ix4=x34+cosb
      iy3=y34+sinb
      iy4=y34-sinb
!      call drv (ixi2,ieta2,ix3,iy3)
!      call drv (ixi2,ieta2,ix4,iy4)
  210 ij=ij+nq2
c
c * * draw constant plane free surface
c
      if (iper.ne.1) go to 220
      call surcntr (i1,i2,j1,j2,kk1,kk2,mpln,axp,ayp)
c
c * * write titles and info on plot
c
  220 continue
      xbl=x(1)
      ybf=y(1)
c      call lincnt (59)
      go to (230,240,250), mpln
  230 write (12,350) vmax,vmaxi
      go to 260
  240 write (12,360) vmax,vmaxj
      go to 260
  250 write (12,370) vmax,vmaxk
  260 continue
      write (12,380) cycle,t,(name(i),i=1,5),jnm
      if (iper.eq.2) go to 300
      go to (270,280,290), mpln
  270 write (12,320) ifm,j1,j2,kk1,kk2
      go to 300
  280 write (12,330) jfm,i1,i2,kk1,kk2
      go to 300
  290 write (12,340) kfm,i1,i2,j1,j2
  300 return
c
c * * error - stop and print
c
  310 continue
      write (9,390)
      write (59,390)
      write (12,390)
!      call exita (3020)
	stop !main
c
  320 format (13h  constant i=,i4,17h   surface     j=,i3,3h to,i3,7h
     1  k=,i3,3h to,i3)
  330 format (13h  constant j=,i4,17h   surface     i=,i3,3h to,i3,7h
     1  k=,i3,3h to,i3)
  340 format (13h  constant k=,i4,17h   surface     i=,i3,3h to,i3,7h
     1  j=,i3,3h to,i3)
  350 format (18h  velocity vectors,4x,7h  vmax=,1pe12.5,4x,6hvmaxi=,e12
     1 .5)
  360 format (18h  velocity vectors,4x,7h  vmax=,1pe12.5,4x,6hvmaxj=,e12
     1 .5)
  370 format (18h  velocity vectors,4x,7h  vmax=,1pe12.5,4x,6hvmaxk=,e12
     1 .5)
  380 format (8h  cycle=,i5,4h  t=,1pe12.5,5x,5a8,1x,a8)
  390 format (51h0 - - - - error - - twx array for plotting exceeded)
      end
*dk vchgcal
      subroutine vchgcal
*ca slcom1
c
c * * calculate total volume change
      include "param.fi"
      voftot=0.0
c
      do 40 k=2,km1
      do 40 j=2,jm1
      do 40 i=2,im1
      call calcijk
      if (beta(ijk).lt.0.0) go to 40
      voftot=voftot+f(ijk)*delx(i)*dely(j)/rri(i)*delz(k)*ac(ijk)
      vchg=0.0
      if (f(ijk).gt.emf.and.f(ijk).lt.emf1) go to 20
      if (f(ijk).ge.emf1) go to 10
      vchg=f(ijk)
      f(ijk)=0.0
      go to 20
   10 continue
      vchg=-(1.0-f(ijk))
      f(ijk)=1.0
   20 continue
      vchgt=vchgt+vchg*delx(i)*dely(j)/rri(i)*delz(k)*ac(ijk)
      if (f(ijk).lt.emf1) go to 40
      if (f(ipjk).lt.emf) go to 30
      if (f(imjk).lt.emf) go to 30
      if (f(ijpk).lt.emf) go to 30
      if (f(ijmk).lt.emf) go to 30
      if (f(ijkp).lt.emf) go to 30
      if (f(ijkm).lt.emf) go to 30
      go to 40
   30 f(ijk)=f(ijk)-1.1*emf
      vchg=1.1*emf
      vchgt=vchgt+vchg*delx(i)*dely(j)/rri(i)*delz(k)*ac(ijk)
   40 continue
c
      return
      end
*dk vfconv
      subroutine vfconv
*ca slcom1
c
c * * convect volume of fluid function f
c
      include "param.fi"
      data flgcs, kount /0.0,0/
      if (cycle.lt.1) go to 100
      kount=0
      flgcs=flgcs+flgc
      flgc=0.0
      do 90 k=1,km1
      do 90 j=1,jm1
      do 90 i=1,im1
      call calcijk
      vx=u(ijk)*delt
      vy=v(ijk)*delt
      vz=w(ijk)*delt
      abvx=abs(vx)
      abvy=abs(vy)
      abvz=abs(vz)
      if (nf(ijk).gt.7) go to 20
      if (beta(ijk).lt.0.0) go to 20
      if (abvx.lt.0.9*delx(i).and.abvy.lt.0.9*dely(j)/rri(i).and.abvz.lt
     1 .0.9*delz(k)) go to 20
      flgc=1.0
      write (9,130) i,j,k,u(ijk),v(ijk),w(ijk),t,delt,cycle,iter
      write (12,130) i,j,k,u(ijk),v(ijk),w(ijk),t,delt,cycle,iter
      if (flgcs.lt.100.0) go to 10
      write (12,110) flgcs
      write (9,110) flgcs
      call exit
   10 continue
      kount=kount+1
      if (kount.lt.20) go to 20
      write (12,120) kount
      write (9,120) kount
      call exit
   20 if (beta(ijk).le.0.0) go to 90
      if (beta(ipjk).le.0.0) go to 40
c
c * * convect f in x-direction
c
      ia=i+1
      id=i
      iajk=ipjk
      idjk=ijk
      idm=max0(i-1,1)
      rb=ar(ijk)*(1.0-cyl+cyl*x(i))
      ra=ac(ipjk)*(1.0-cyl+cyl*xi(i+1))
      rd=ac(ijk)*(1.0-cyl+cyl*xi(i))
      if (vx.ge.0.0) go to 30
      ia=i
      id=i+1
      iajk=ijk
      idjk=ipjk
      idm=min0(i+2,imax)
      ra=ac(ijk)*(1.0-cyl+cyl*xi(i))
      rd=ac(ipjk)*(1.0-cyl+cyl*xi(i+1))
   30 continue
      iad=ia
      idmjk=nq*(ii5*(k-1)+imax*(j-1)+(idm-1))+1
      if (nf(idjk).eq.3.or.nf(idjk).eq.4) iad=id
      if (nf(idjk).eq.5.or.nf(idjk).eq.6) iad=id
      if (fn(iajk).lt.emf.or.fn(idmjk).lt.emf) iad=ia
      iadjk=nq*(ii5*(k-1)+imax*(j-1)+(iad-1))+1
      fdm=amax1(fn(idmjk),fn(idjk))
      if (f(idmjk).lt.emf.and.f(iajk).lt.emf) fdm=amax1(fdm,0.10)
      if (beta(idmjk).lt.0.0) fdm=1.0
      fx1=fn(iadjk)*abvx+amax1((fdm-fn(iadjk))*abvx-(fdm-fn(idjk))*delx
     1 (id),0.0)
      fx=amin1(fx1,fn(idjk)*delx(id))
      f(idjk)=f(idjk)-fx*rdx(id)*(rb/rd)
      f(iajk)=f(iajk)+fx*rdx(ia)*(rb/ra)
   40 if (beta(ijpk).le.0.0) go to 60
c
c * * convect f in y-direction
c
      ja=j+1
      jd=j
      ijak=ijpk
      ijdk=ijk
      jdm=max0(j-1,1)
      rb=abk(ijk)
      ra=ac(ijpk)
      rd=ac(ijk)
      if (vy.ge.0.0) go to 50
      ja=j
      jd=j+1
      ijak=ijk
      ijdk=ijpk
      jdm=min0(j+2,jmax)
      ra=ac(ijk)
      rd=ac(ijpk)
   50 continue
      jad=ja
      ijdmk=nq*(ii5*(k-1)+imax*(jdm-1)+(i-1))+1
      if (nf(ijdk).eq.5.or.nf(ijdk).eq.6) jad=jd
      if (nf(ijdk).eq.1.or.nf(ijdk).eq.2) jad=jd
      if (fn(ijak).lt.emf.or.fn(ijdmk).lt.emf) jad=ja
      ijadk=nq*(ii5*(k-1)+imax*(jad-1)+(i-1))+1
      fdm=amax1(fn(ijdmk),fn(ijdk))
      if (f(ijdmk).lt.emf.and.f(ijak).lt.emf) fdm=amax1(fdm,0.10)
      if (beta(ijdmk).lt.0.0) fdm=1.0
      fy1=fn(ijadk)*abvy+amax1((fdm-fn(ijadk))*abvy-(fdm-fn(ijdk))*dely
     1 (jd)/rri(i),0.0)
      fy=amin1(fy1,fn(ijdk)*dely(jd)/rri(i))
      f(ijdk)=f(ijdk)-fy*rdy(jd)*rri(i)*(rb/rd)
      f(ijak)=f(ijak)+fy*rdy(ja)*rri(i)*(rb/ra)
c
c * * convect f in z-direction
c
   60 if (beta(ijkp).le.0.0) go to 80
      ka=k+1
      kd=k
      ijka=ijkp
      ijkd=ijk
      kdm=max0(k-1,1)
      rb=at(ijk)
      ra=ac(ijkp)
      rd=ac(ijk)
      if (vz.ge.0.0) go to 70
      ka=k
      kd=k+1
      ijka=ijk
      ijkd=ijkp
      kdm=min0(k+2,kmax)
      ra=ac(ijk)
      rd=ac(ijkp)
   70 continue
      kad=ka
      ijkdm=nq*(ii5*(kdm-1)+imax*(j-1)+(i-1))+1
      if (nf(ijkd).eq.1.or.nf(ijkd).eq.2) kad=kd
      if (nf(ijkd).eq.3.or.nf(ijkd).eq.4) kad=kd
      if (fn(ijka).lt.emf.or.fn(ijkdm).lt.emf) kad=ka
      ijkad=nq*(ii5*(kad-1)+imax*(j-1)+(i-1))+1
      fdm=amax1(fn(ijkdm),fn(ijkd))
      if (f(ijkdm).lt.emf.and.f(ijka).lt.emf) fdm=amax1(fdm,0.10)
      if (beta(ijkdm).lt.0.0) fdm=1.0
      fz1=fn(ijkad)*abvz+amax1((fdm-fn(ijkad))*abvz-(fdm-fn(ijkd))*delz
     1 (kd),0.0)
      fz=amin1(fz1,fn(ijkd)*delz(kd))
      f(ijkd)=f(ijkd)-fz*rdz(kd)*(rb/rd)
      f(ijka)=f(ijka)+fz*rdz(ka)*(rb/ra)
   80 if (nf(ijk).ne.0) go to 90
      if (idefm.gt.0) go to 90
      if (beta(ijk).eq.1.0) go to 90
      f(ijk)=f(ijk)+delt*fn(ijk)*d(ijk)
   90 continue
c
c * * calculate total volume change
c
  100 call vchgcal
c
      return
c
  110 format ('vfconv - too many large velocities  flgcs=',f10.5)
  120 format ('vfconv - too many large velocities  kount=',f10.5)
  130 format ('vfconv - velocities too large'/3x,3i3,5f12.5,2i5)
      end
*dk wrtape
      subroutine wrtape
*ca slcom1
c
c * * write magmetic tape
c
      include "param.fi"
      write (9,10) numtd,t,cycle
      if (lpr.gt.0) write (12,10) numtd,t,cycle
      write (59,10) numtd,t,cycle
      twtd=t+tddt
      write (8) (aa(n),n=1,ncr1)
      write (8) (basc(n),n=1,ncr2)
      write (8) (basc1(n),n=1,ncr2)
      write (8) (basc2(n),n=1,ncr2)
      write (8) (basc3(n),n=1,ncr2)
      numtd=numtd+1
c
      return
c
   10 format (12h0 tape dump=,i3,7h  at t=,1pe12.5,7h cycle=,i5)
      end
*dk xindf
      subroutine xindf
*ca slcom1
c
c * * set arrays to indefinite
c
      include "param.fi"
      ilim=ii5*kmax*nq
      xindef=1e10            !data xindef /17770000000000777777b/
      iup=ibasc
      do 10 i=1,iup
      basc(i)=xindef
      basc2(i)=xindef
   10 basc3(i)=xindef
      ncr1=1	!locf(zlast)-locf(aa)+1
      ncr2=1	!locf(basc(ilim))-locf(basc)+1
c
      return
      end
*dk drawq
      subroutine drawq
*ca slcom1
c
      include "param.fi"
      parameter (ibar2q=2*ibar2-2)
      common betaq(ibar2q,kbar2), fq(ibar2q,kbar2), acq(ibar2q,kbar2),
     1 arq(ibar2q,kbar2), atq(ibar2q,kbar2), uq(ibar2q,kbar2), vq(ibar2q
     2 ,kbar2), xq(ibar2q), xiq(ibar2q), yq(kbar2), yjq(kbar2), im1q,
     3 jm1q, jbar2q, jmaxq, imaxq, sf, xshft, yshft, delxq(ibar2q),
     4 delyq(kbar2), xminq, xmaxq, yminq, ymaxq
c
c      move quanities from 3d mesh to 2d plot mesh
c
c      plotting j=2 and j=jm1 planes for 180 deg problem
c      plotting j=2 and j=jop(2) planes for 360 deg problem
c
c      im1q=ibar2q-1
c      jm1q=km1
c      jbar2q=kbar2
c      jmaxq=kmax
c      imaxq=ibar2q
c
c      j=2
c      do 10 l=imax,ibar2q
c      do 10 k=1,kmax
c      i=l-imax+2
c      call ijkonly
c      betaq(l,k)=beta(ijk)
c      fq(l,k)=f(ijk)
c      arq(l,k)=ar(ijk)
c      atq(l,k)=at(ijk)
c      acq(l,k)=ac(ijk)
c      uq(l,k)=u(ijk)
c      vq(l,k)=w(ijk)
c   10 continue
c
c      j=jop(2)
c
c      do 20 l=1,im1
c      do 20 k=1,kmax
c      i=imax-l+1
c      call ijkajct
c      betaq(l,k)=beta(ijk)
c      fq(l,k)=f(ijk)
c      arq(l,k)=ar(imjk)
c      if (l.eq.im1) arq(l,k)=1.0
c      atq(l,k)=at(ijk)
c      acq(l,k)=ac(ijk)
c      uq(l,k)=-u(imjk)
c      if (l.eq.im1) uq(l,k)=0.0
c      vq(l,k)=w(ijk)
c   20 continue
c
c      do 30 l=imax,ibar2q
c      i=l-im1+1
c      xq(l)=x(i)+x(im1)
c      delxq(l)=delx(i)
c   30 continue
c
c      do 40 l=1,im1
c      i=imax-l
c      xq(l)=x(im1)-x(i)
c      delxq(l)=delx(i)
c   40 continue
c
c      do 50 l=2,ibar2q
c      xiq(l)=(xq(l)+xq(l-1))*0.5
c   50 continue
c      xiq(1)=-xiq(2)
c
c      do 60 k=1,kmax
c      yq(k)=z(k)
c      yjq(k)=zk(k)
c      delyq(k)=delz(k)
c   60 continue
c
c +++ set constant terms for plotting
c
c      xminq=xq(1)
c      xmaxq=xq(im1q)
c      yminq=yq(1)
c      ymaxq=yq(jm1q)
c      d1=xmaxq-xminq
c      d2=ymaxq-yminq
c      d3=amax1(d1,d2)
c      sf=1.0/d3
c      xshft=0.5*(1.0-d1*sf)
c      yshft=0.5*(1.0-d2*sf)
c
c * * determine maximum vector length
c
c      vmax=0.0
c      vmaxj=0.0
c      do 70 k=2,km1
c      do 70 j=2,jm1
c      do 70 i=2,im1
c      call ijkajct
c      if (f(ijk).lt.0.5) go to 70
c      if (ac(ijk).lt.0.5) go to 70
c      if (j.ne.2.and.j.ne.jop(2)) go to 70
c      uave=0.5*(u(ijk)+u(imjk))
c      vave=0.5*(v(ijk)+v(ijmk))
c      wave=0.5*(w(ijk)+w(ijkm))
c      vmax=amax1(vmax,abs(uave),abs(vave),abs(wave))
c      vmaxj=amax1(vmaxj,abs(uave),abs(wave))
c   70 continue
c
c * * scale vector length
c
c      velmx1=velmx*delmn/(vmaxj+1.0e-10)
c
c +++
c +++ draw velocity vector plot
c +++
c      call adv (1)
c      call lincnt (1)
c     write(12,200) t,cycle,name
c      call frameq (xminq,xmaxq,ymaxq,yminq)
c      call drwobsq
c      do 80 i=2,im1q
c     acelh=0.5*cyl+1.0-cyl
c      do 80 j=2,jm1q
c      if (fq(i,j).lt.0.5) go to 80
c      if (acq(i,j).lt.acelh) go to 80
c      xxcc=xiq(i)
c      yycc=0.5*(yq(j)+yq(j-1))
c      umpl=0.0
c      vmpl=0.0
c      if (arq(i,j).gt.em6) umpl=1.0
c      if (arq(i-1,j).gt.em6) umpl=1.0/(1.0+umpl)
c      if (atq(i,j).gt.em6) vmpl=1.0
c      if (atq(i,j-1).gt.em6) vmpl=1.0/(1.0+vmpl)
c      uplt=umpl*(uq(i-1,j)+uq(i,j))
c      vplt=vmpl*(vq(i,j-1)+vq(i,j))
c      uvec=uplt*velmx1+xxcc
c      vvec=vplt*velmx1+yycc
c      call drwveca (xxcc,yycc,uvec,vvec,1)
c      call pltpt (xxcc,yycc,53b,1)
c   80 continue
c
c +++ draw free surface
c
c      fpl=0.5
c      do 100 i=2,im1q
c      do 100 j=2,jm1q
c      if (betaq(i,j).lt.0.0) go to 100
c      fatr=0.25*(fq(i,j)+fq(i+1,j)+fq(i,j+1)+fq(i+1,j+1))
c      fxtr=0.5*(fq(i+1,j+1)+fq(i+1,j)-fq(i,j+1)-fq(i,j))/(xiq(i+1)-xiq(i
c     1 ))
c      fytr=0.5*(fq(i,j+1)+fq(i+1,j+1)-fq(i,j)-fq(i+1,j))/(yjq(j+1)-yjq(j
c     1 ))
c      ftrs=fxtr**2+fytr**2
c      if (ftrs.eq.0.0) ftrs=1.0e+10
c      xtr=0.5*(xiq(i+1)+xiq(i))+(fpl-fatr)*fxtr/ftrs
c      xtr=amax1(xtr,xiq(i))
c      xtr=amin1(xtr,xiq(i+1))
c      ytr=0.5*(yjq(j)+yjq(j+1))+(fpl-fatr)*fytr/ftrs
c      ytr=amax1(ytr,yjq(j))
c      ytr=amin1(ytr,yjq(j+1))
c      if (fq(i,j).gt.0.5.and.fq(i+1,j).gt.0.5) go to 90
c      if (fq(i,j).lt.0.5.and.fq(i+1,j).lt.0.5) go to 90
c      fabr=0.25*(fq(i,j)+fq(i+1,j)+fq(i,j-1)+fq(i+1,j-1))
c      fxbr=0.5*(fq(i+1,j)+fq(i+1,j-1)-fq(i,j)-fq(i,j-1))/(xiq(i+1)-xiq(i
c     1 ))
c      fybr=0.5*(fq(i,j)+fq(i+1,j)-fq(i,j-1)-fq(i+1,j-1))/(yjq(j)-yjq(j-1
c     1 ))
c      fbrs=fxbr**2+fybr**2
c      if (fbrs.eq.0.0) fbrs=1.0e+10
c      xxbr=0.5*(xiq(i+1)+xiq(i))+(fpl-fabr)*fxbr/fbrs
c      xxbr=amax1(xxbr,xiq(i))
c      xxbr=amin1(xxbr,xiq(i+1))
c      ybr=0.5*(yjq(j)+yjq(j-1))+(fpl-fabr)*fybr/fbrs
c      ybr=amax1(ybr,yjq(j-1))
c      ybr=amin1(ybr,yjq(j))
c      call drwvec (xxbr,ybr,xtr,ytr,1)
c   90 continue
c      if (fq(i,j).gt.0.5.and.fq(i,j+1).gt.0.5) go to 100
c      if (fq(i,j).lt.0.5.and.fq(i,j+1).lt.0.5) go to 100
c      fatl=0.25*(fq(i,j)+fq(i,j+1)+fq(i-1,j)+fq(i-1,j+1))
c      fxtl=0.5*(fq(i,j+1)+fq(i,j)-fq(i-1,j+1)-fq(i-1,j))/(xiq(i)-xiq(i-1
c     1 ))
c      fytl=0.5*(fq(i-1,j+1)+fq(i,j+1)-fq(i-1,j)-fq(i,j))/(yjq(j+1)-yjq(j
c     1 ))
c      ftls=fxtl**2+fytl**2
c      if (ftls.eq.0.0) ftls=1.0e+10
c      xtl=0.5*(xiq(i-1)+xiq(i))+(fpl-fatl)*fxtl/ftls
c      xtl=amax1(xtl,xiq(i-1))
c      xtl=amin1(xtl,xiq(i))
c      ytl=0.5*(yjq(j)+yjq(j+1))+(fpl-fatl)*fytl/ftls
c      ytl=amax1(ytl,yjq(j))
c      ytl=amin1(ytl,yjq(j+1))
c      call drwvec (xtl,ytl,xtr,ytr,1)
c  100 continue
c
c      call lincnt (59)
c      write (12,120) vmax,vmaxj
c      write (12,130) cycle,t,(name(i),i=1,5),jnm
c      write (12,110) 2,im1,2,km1
      return
c
c  110 format (21h constant j=2 and jm1,17h   surface     i=,i3,3h to,i3,
c     1 7h     k=,i3,3h to,i3)
c  120 format (18h  velocity vectors,4x,7h  vmax=,1pe12.5,4x,6hvmaxj=,e12
c     1 .5)
c  130 format (8h  cycle=,i5,4h  t=,1pe12.5,5x,5a8,1x,a8)
      end
*dk drwobsq
      subroutine drwobsq
*ca slcom1
c
      include "param.fi"
      parameter (ibar2q=2*ibar2-2)
      common betaq(ibar2q,kbar2), fq(ibar2q,kbar2), acq(ibar2q,kbar2),
     1 arq(ibar2q,kbar2), atq(ibar2q,kbar2), uq(ibar2q,kbar2), vq(ibar2q
     2 ,kbar2), xq(ibar2q), xiq(ibar2q), yq(kbar2), yjq(kbar2), im1q,
     3 jm1q, jbar2q, jmaxq, imaxq, sf, xshft, yshft, delxq(ibar2q),
     4 delyq(kbar2), xminq, xmaxq, yminq, ymaxq
c
c +++ draw around all obstacles
c
c      do 170 i=2,im1q
c      atr=1.0-em6
c      atl=1.0-em6
c      if (i.eq.2.) atl=-em6
c      atc=1.0-em6
c      do 170 j=2,jm1q
c      if (acq(i,j).lt.em6) go to 170
c      afr=1.0
c      aft=1.0
c      afl=1.0
c      afb=1.0
c      if (arq(i,j).lt.atr) afr=arq(i,j)/atr
c      if (atq(i,j).lt.atc) aft=atq(i,j)/atc
c      if (arq(i-1,j).lt.atl) afl=arq(i-1,j)/atl
c      if (atq(i,j-1).lt.atc) afb=atq(i,j-1)/atc
c      if (acq(i,j).ge.atc) go to 120
c      if (i.eq.2) afl=afr-em6
c      if (i.eq.im1q) afr=afl-em6
c      if (j.eq.2) afb=aft-em6
c      if (j.eq.jm1q) aft=afb-em6
c      if ((aft+afb).lt.em6.or.(afl+afr).lt.em6) go to 170
c      m=1
c      amn=afb+afr
c      if ((afr+aft).gt.amn) go to 10
c      m=2
c      amn=afr+aft
c   10 if ((aft+afl).gt.amn) go to 20
c      m=3
c      amn=aft+afl
c   20 if ((afl+afb).gt.amn) go to 30
c      m=4
c   30 go to (40,60,80,100), m
c   40 x1=xq(i-1)+aft*delxq(i)
c      y1=yq(j)
c      if (aft.lt.1.0) go to 50
c      y1=y1-afr*delyq(j)
c   50 x2=xq(i-1)
c      zy2=yq(j)-afl*delyq(j)
c      if (afl.lt.1.0) go to 160
c      x2=x2+afb*delxq(i)
c      go to 160
c   60 x1=xq(i-1)
c      y1=yq(j-1)+afl*delyq(j)
c      if (afl.lt.1.0) go to 70
c      x1=x1+aft*delxq(i)
c   70 x2=xq(i-1)+afb*delxq(i)
c      zy2=yq(j-1)
c      if (afb.lt.1.0) go to 160
c      zy2=zy2+afr*delyq(j)
c      go to 160
c   80 x1=xq(i)-afb*delxq(i)
c      y1=yq(j-1)
c      if (afb.lt.1.0) go to 90
c      y1=y1+afl*delyq(j)
c   90 x2=xq(i)
c      zy2=yq(j-1)+afr*delyq(j)
c      if (afr.lt.1.0) go to 160
c      x2=x2-aft*delxq(i)
c      go to 160
c  100 x1=xq(i)
c      y1=yq(j)-afr*delyq(j)
c      if (afr.lt.1.0) go to 110
c      x1=x1-afb*delxq(i)
c  110 x2=xq(i)-aft*delxq(i)
c      zy2=yq(j)
c      if (aft.lt.1.0) go to 160
c      zy2=zy2-afl*delyq(j)
c      go to 160
c  120 if (afr.gt.em6) go to 130
c      x1=xq(i)
c      y1=yq(j-1)
c      x2=x1
c      zy2=yq(j)
c      call drwvec (x1,y1,x2,zy2,1)
c  130 if (aft.gt.em6) go to 140
c      x1=xq(i-1)
c      y1=yq(j)
c      x2=xq(i)
c      zy2=y1
c      call drwvec (x1,y1,x2,zy2,1)
c  140 if (afl.gt.em6) go to 150
c      x1=xq(i-1)
c      y1=yq(j)
c      x2=x1
c      zy2=yq(j-1)
c      call drwvec (x1,y1,x2,zy2,1)
c  150 if (afb.gt.em6) go to 170
c      x1=xq(i-1)
c      y1=yq(j-1)
c      x2=xq(i)
c      zy2=y1
c      call drwvec (x1,y1,x2,zy2,1)
c  160 call drwvec (x1,y1,x2,zy2,1)
c  170 continue
c      return
      end
*dk drwvec
      subroutine drwvec (xone,yone,xtwo,ytwo,isym)
*ca slcom1
      include "param.fi"
      parameter (ibar2q=2*ibar2-2)
      common betaq(ibar2q,kbar2), fq(ibar2q,kbar2), acq(ibar2q,kbar2),
     1 arq(ibar2q,kbar2), atq(ibar2q,kbar2), uq(ibar2q,kbar2), vq(ibar2q
     2 ,kbar2), xq(ibar2q), xiq(ibar2q), yq(kbar2), yjq(kbar2), im1q,
     3 jm1q, jbar2q, jmaxq, imaxq, sf, xshft, yshft, delxq(ibar2q),
     4 delyq(kbar2), xminq, xmaxq, yminq, ymaxq
c
c +++ draw a vector
c +++ provides a system dependant call
c
c      ic=0
c      x1=xone
c      y1=yone
c      x2=xtwo
c      zy2=ytwo
c     x01=(x1-xminq)*sf+xshft
c      y01=(y1-yminq)*sf+yshft
c      x02=(x2-xminq)*sf+xshft
c      y02=(zy2-yminq)*sf+yshft
c     ix1=16.+900.0*x01
c      ix2=16.+900.0*x02
c      iy1=16.+900.0*(1.0-y01)
c      iy2=16.+900.0*(1.0-y02)
c      call drv (ix1,iy1,ix2,iy2)
      return
c
      entry drwveca(xone,yone,xtwo,ytwo,isym)
c
c +++ draw a vector with arrowheads
c
c      ic=0
c      x1=xone
c      y1=yone
c      x2=xtwo
c      zy2=ytwo
c      x01=(x1-xminq)*sf+xshft
c      y01=(y1-yminq)*sf+yshft
c      x02=(x2-xminq)*sf+xshft
c      y02=(zy2-yminq)*sf+yshft
c      ix1=16.+900.0*x01
c      ix2=16.+900.0*x02
c      iy1=16.+900.0*(1.0-y01)
c      iy2=16.+900.0*(1.0-y02)
c      call drv (ix1,iy1,ix2,iy2)
c    draw arrowheads
c      fix1=ix1
c      fix2=ix2
c      fiy1=iy1
c      fiy2=iy2
c      x34=fix1+0.85*(fix2-fix1)
c      y34=fiy1+0.85*(fiy2-fiy1)
c      dx234=fix2-x34
c      dy234=fiy2-y34
c      d234=sqrt(dx234**2+dy234**2)
c      da234=d234*0.57735027
c      tan(15 deg)=0.57735027
c      cosb=da234
c      sinb=0.0
c      if (dx234.eq.0.0) go to 10
c      beta234=1.57079633-atan(dy234/dx234)
c      cosb=cos(beta234)*da234
c      sinb=sin(beta234)*da234
c   10 ix3=x34-cosb
c      ix4=x34+cosb
c      iy3=y34+sinb
c      iy4=y34-sinb
c      call drv (ix2,iy2,ix3,iy3)
c      call drv (ix2,iy2,ix4,iy4)
      return
      end
*dk frame
      subroutine frameq (xxl,xxr,yyt,yyb)
c
c +++ draw a frame around the plot
c
c      call drwvec (xxl,yyt,xxr,yyt,0)
c      call drwvec (xxl,yyt,xxl,yyb,0)
c      call drwvec (xxl,yyb,xxr,yyb,0)
c      call drwvec (xxr,yyb,xxr,yyt,0)
      return
      end


	SUBROUTINE SETPAR
	INCLUDE 'PARAM.FI'	
      
      DELT=0.01
      FLHT=0.0
      GY=-9.81
      PLTDT=0.1
      PRTDT=0.1
      TWFIN=1.0
      UI=0.0
      VELMX=5.0
      WL=2
      WR=2
      WT=2
      WB=2
      AUTOT=1.0
	
	
	
C	*** MESH SETUP
	NKX=1
	XL(1)=0.0
	XL(2)=1.0
	XC(1)=0.5
	
	NXL=8
	NXR=8
	DXMN(1)=1.0/16.0

	NKY=1
	YL(1)=0.0
	YL(2)=0.1
	YC(1)=0.5
	NYL=1
	NYR=1
	DYMN(1)=0.05
	
	NKZ=1
	ZL(1)=0.0
	ZL(2)=1.0
	ZC(1)=0.5
	NZL=8
	NZR=8
	DZMN(1)=1.0/16.0

C     PHYSICS SETUP
      NU=1E-6
      RHOF=1000.0
      
	END

	SUBROUTINE VTK(IOPST,FNAME)
      INCLUDE 'PARAM.FI'	
      
      INTEGER IOPST,IOS
      CHARACTER(200) FNAME
      REAL UTEMP(IBAR2,JBAR2),VTEMP(IBAR2,JBAR2)
      INTEGER I,J
      

	open(IOpst,file=fname,status='new',iostat=ios)
	write(IOpst,'(A)')   '# vtk DataFile Version 2.0'
	write(IOpst,'(A,A,A)')   'vof result output'
	write(IOpst,'(A)')      'ASCII'

	!grid data

	
      END

