      program platem
      implicit double precision (a-h,o-z)
      include 'dftpol.id.inc'      
      dimension c(maxmon,0:maxel),
     *etheta(0:maxel),ch(maxmon,0:maxel),fdgr(0:maxel),
     *Pn(1001),hh(1001),gp(1001)
      write(*,*) 'Monomer 1 constrained to be within'        
      write(*,*) 'one sigma from z=closew, i.e the left wall'  
      write(*,*) '(flexible grafting)'    
      ifc = 38
      ins = 49
      isd = 61
      pi = acos(-1.d0)
      bk = 1.38066D-23
      avno = 6.02214D23
      fourpi = 4.d0*pi
      twopi = 2.d0*pi
      fourpi = 4.d0*pi
      rtwelve = 1.d0/12.d0
      rthree = 1.d0/3.d0
      volfact = fourpi/3.d0
      rnine = 1.d0/9.d0
      rthree = 1.d0/3.d0
      rphi = fourpi*0.2d0 
      aphi = fourpi*0.5d0
      pis = pi/6.d0
      pit = pi/3.d0
      pif = pi/4.d0
      ddtol = 0.0000001d0
      open (ifc,file='densitydistribution',form='formatted')    
      open (ins,file='input',form='formatted')
      rewind ifc
      rewind ins
      T = 1.d0
      dz = 0.05d0
      closew = 0.d0
      read(ins,*) 
      read(ins,*) nmon
      read(ins,*) 
      read(ins,*) surfdens
      read(ins,*) 
      read(ins,*) hstart,dh,nh
      sclosew = closew
      surfdens = dabs(surfdens)
      write(*,*) 'grafting density: ',surfdens
      rrT = 1.d0/T
      rdz = 1.d0/dz

c      twopidz = twopi*dz
      twopidz = 0.5d0*dz

      irdz = int(rdz+0.001d0)
      istart = int(closew/dz+0.01d0)
      istp1 = istart+1 
      ists = int(sclosew/dz+0.01d0)
      istp1s = ists+1 
      ism = int(1.d0/dz+0.01d0)
      inw = ism+int(closew/dz+0.01d0)
      pie = pi/8.d0
      dzpie = pie*dz
      rnmon = real(nmon)
      rrnmon = 1.d0/rnmon      
      rrcmon = 1.d0/(rnmon-2.d0)
      checknm = abs(rnmon*0.5d0-real(int(rnmon*0.5d0+0.0001)))
      write(*,*) 'no. of monomers/polymer: ',nmon

      h = hstart-dh
      do ih = 1,nh
      h = h+dh
      nfack = int(h/dz+0.01d0)
      islut = int((h-closew)/dz+0.01d0)      
      isluts = int((h-sclosew)/dz+0.01d0)  
      imitt = nfack/2    
      do 36 iz = istp1,islut
      emtrams = 0.d0
      cmtrams = 0.d0
      ebelam(iz) = 1.d0
 36   ehbclam(iz) = 1.d0
      do i = istp1,istp1+ism-1
      etheta(i) = ebelam(i)
      enddo
      do  245 iz = istp1,inw
      sume = 0.d0
      rsume = 0.d0
      do 345 jz = istp1,iz+ism-1
      rsume = etheta(jz)+rsume
 345  sume = ebelam(jz)+sume
      trams = ehbclam(iz)*twopidz
      ch(nmon-1,iz) = (0.5d0*ebelam(iz+ism)+sume)*trams
 245  c(nmon-1,iz) = (0.5d0*etheta(iz+ism)+rsume)*trams
      do 445 iz = inw+1,islut-ism
      sume = 0.5d0*ebelam(iz-ism)
      rsume = 0.5d0*etheta(iz-ism)
      do 545 jz = iz-ism+1,iz+ism-1
      rsume = etheta(jz)+rsume
 545  sume = ebelam(jz)+sume
      trams = ehbclam(iz)*twopidz
      ch(nmon-1,iz) = (0.5d0*ebelam(iz+ism)+sume)*trams
 445  c(nmon-1,iz) = (0.5d0*etheta(iz+ism)+rsume)*trams
      do iz = islut-ism+1,islut
      sume = 0.5d0*ebelam(iz-ism)
      rsume = 0.5d0*etheta(iz-ism)
      do jz = iz-ism+1,islut
      rsume = etheta(jz)+rsume
      sume = ebelam(jz)+sume
      enddo
      ch(nmon-1,iz) = ehbclam(iz)*sume*twopidz
      c(nmon-1,iz) = etheta(iz)*rsume*twopidz
      enddo

      k = nmon-1
      do 745 mmm = 2,nmon-2 
      k = k-1
      do  845 iz = istp1,inw
      sume = 0.d0
      do 945 jz = istp1,iz+ism-1
 945  sume = ehbclam(jz)*c(k+1,jz)+sume
      c(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*c(k+1,iz+ism)+
     *sume)*twopidz
      sume = 0.d0
      do jz = istp1,iz+ism-1
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
 845  ch(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(k+1,iz+ism)+
     *sume)*twopidz
      do 1045 iz = inw+1,islut-ism
      sume = 0.5d0*ehbclam(iz-ism)*c(k+1,iz-ism)
      do 1145 jz = iz-ism+1,iz+ism-1
 1145 sume = ehbclam(jz)*c(k+1,jz)+sume
      c(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*c(k+1,iz+ism)+
     *sume)*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(k+1,iz-ism)
      do jz = iz-ism+1,iz+ism-1
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
 1045 ch(k,iz) = ehbclam(iz)*(0.5d0*ehbclam(iz+ism)*ch(k+1,iz+ism)+
     *sume)*twopidz
      do iz = islut-ism+1,islut
      sume = 0.5d0*ehbclam(iz-ism)*c(k+1,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*c(k+1,jz)+sume
      enddo
      c(k,iz) = ehbclam(iz)*sume*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(k+1,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*ch(k+1,jz)+sume
      enddo
      ch(k,iz) = ehbclam(iz)*sume*twopidz
      enddo
 745  continue

      do  1245 iz = istp1,inw
      sume = 0.d0
      do 1345 jz = istp1,iz+ism-1
 1345 sume = ehbclam(jz)*c(2,jz)+sume
      c(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*c(2,iz+ism)+sume)*
     *twopidz
      sume = 0.d0
      do 8645 jz = istp1,iz+ism-1
 8645 sume = ehbclam(jz)*ch(2,jz)+sume
 1245 ch(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
     *twopidz
      do 1445 iz = inw+1,islut-ism
      sume = 0.5d0*ehbclam(iz-ism)*c(2,iz-ism)
      do 1545 jz = iz-ism+1,iz+ism-1
 1545 sume = ehbclam(jz)*c(2,jz)+sume
      c(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*c(2,iz+ism)+sume)*
     *twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(2,iz-ism)
      do 8745 jz = iz-ism+1,iz+ism-1
 8745 sume = ehbclam(jz)*ch(2,jz)+sume
 1445 ch(1,iz) = ebelam(iz)*(0.5d0*ehbclam(iz+ism)*ch(2,iz+ism)+sume)*
     *twopidz
      do iz = islut-ism+1,islut
      sume = 0.5d0*ehbclam(iz-ism)*c(2,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*c(2,jz)+sume
      enddo
      c(1,iz) = ebelam(iz)*sume*twopidz
      sume = 0.5d0*ehbclam(iz-ism)*ch(2,iz-ism)
      do jz = iz-ism+1,islut
      sume = ehbclam(jz)*ch(2,jz)+sume
      enddo
      ch(1,iz) = ebelam(iz)*sume*twopidz
      enddo
      sumgrmon = 0.d0
      do i = istp1s,istp1s+ism-1
      sumgrmon = ch(1,i)+sumgrmon
      enddo
      sumgrmon = sumgrmon*dz
      fnorm = surfdens/sumgrmon 
      do i = istp1s,istp1s+ism-1
      fdgr(i) = ch(1,i)*fnorm
      enddo          
      ddmax = 0.d0
      z = -0.5d0*dz
      do 9 i = istp1s,islut
      z = z+dz
      if (z.lt.rnmon) then
      dumsum = 0.d0 
      do 10 k = 2,nmon-1
 10   dumsum = c(k,i)*ch(nmon+1-k,i)+dumsum
      ttfem(i) = c(1,i)*fnorm
      dumsum = dumsum*fnorm
      if (i.lt.(istp1+ism)) then
      tfdmon(i) = dumsum+ttfem(i)+fdgr(i)
      else
      tfdmon(i) = dumsum+ttfem(i)
      endif
      else
      tfdmon(i) = 0.d0
      ttfem(i) = 0.d0
      endif
 9    continue
      z = -0.5d0*dz
      k = islut+1
      do i = istp1,imitt
      k = k-1
      tfdm = tfdmon(i)+tfdmon(k)
      tfem = ttfem(i)+ttfem(k)
      z = z+dz
      if (z.lt.rnmon) then      
      ddiff = abs(tfdm-fdmon(i))/tfdm
      if (ddiff.gt.ddmax) ddmax = ddiff
      fem(i) = tfem
      fdmon(i) = tfdm
      else
      fdmon(i) = 0.d0
      fem(i) = 0.d0
      endif
      enddo
      jz = imitt+1
      do iz = imitt+1,isluts
      jz = jz-1
      fdgr(iz) = fdgr(jz)
      fdmon(iz) = fdmon(jz)
      fem(iz) = fem(jz)
      enddo
      niz = 0
      sumfdm = 0.d0
      Freen = 0.d0
      z = 0.d0
      do 500 iz = istp1,islut
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      niz = niz+1
      fdm = fdmon(iz)      
      fde = fem(iz)
      fdc = fdm-fde-fdgr(iz)
      sumfdm = fdm+sumfdm
 500  Freen = Freen+fdc*bclamb+fde*belamb+fdm*(Vexm(iz)-rrnmon)+
     *fdgr(iz)*(belamb+dlog(fnorm))
      Freen = Freen*dz
       avfdm = sumfdm/real(islut-istart)
      gammam = 0.5d0*(2.*sumfdm*dz-bdm*h)

      Fs2m = 0.d0
      Fslw = 0.d0
      Pmacross = 0.d0
      avfdm = 0.d0
      avfem = 0.d0
      z = sclosew-0.5d0*dz
      do 61 iz = istp1s,isluts
      z = z+dz
      avfdm = fdmon(iz)+avfdm
 61   avfem = fem(iz)+avfem

      write(*,*)
      write(*,*) 'grand potential: '
      write(*,*) Freen
      write(*,*)
      fp1S = fdmon(isluts)
      fn1S = fdmon(isluts-2)
      c0Skv = fdmon(isluts-1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwcmq.lt.0.d0) fwcmq = 0.d0
      fwcml = fp1S+0.5d0*(fp1S-c0Skv)
      if (fwcml.lt.0.d0) fwcmq = 0.d0
      Fs2m = fwcmq-Fs2m*dz
      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwcmq.lt.0.d0) fwcmq = 0.d0
      fwcml = fp1S+0.5d0*(fp1S-c0Skv)
      if (fwcml.lt.0.d0) fwcmq = 0.d0
      Fslw = fwcmq-Fslw*dz

      fp1S = fdmon(istp1+ism)
      fn1S = fdmon(istp1+ism+2)
      c0Skv = fdmon(istp1+ism+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fmdp = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nm(ism+delta) - quad. extr. = ',fmdp
      fp1S = fdmon(istp1+ism-1)
      fn1S = fdmon(istp1+ism-3)
      c0Skv = fdmon(istp1+ism-2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fmdm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nm(izslab-delta) - quad. extr. = ',fmdm

      fp1S = fdgr(istp1+ism-1)
      fn1S = fdgr(istp1+ism-3)
      c0Skv = fdgr(istp1+ism-2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
      fgrbr = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'nmgraft(izslab-delta) - quad. extr. = ',fgrbr
      Pdonn = fmdp-fmdm
      write(*,*)  'Pdonn = ',Pdonn

      Pnlw = Fslw+Pdonn
      Pnet = Fs2m+Pdonn
      write(*,*) 
      write(*,*)'pressure (right wall): '
      write(*,*)  Fs2m-fgrbr
      write(*,*)  Pdonn+fwcmq
      write(*,*)  Fs2m-Pdonn-fgrbr
      write(*,*)  Fs2m+Pdonn+fgrbr
      write(*,*)
      write(*,*)'pressure (left wall): '
      write(*,*)  Fslw-fgrbr
      write(*,*)  Pnlw
      write(*,*)

      Pmacross = -2.d0*Pmacross*dz
      write(*,*) 'monomer contact density (left wall): '   
      write(*,*)  fwcmq
      avfdm = avfdm/real((islut+1-istp1))
      avfem = avfem/real((islut+1-istp1))
      write(*,*) 'avfdm = ',avfdm
      fp1S = fdmon(imitt+1)
      fn1S = fdmon(imitt+3)
      c0Skv = fdmon(imitt+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      write(*,*) 'monomer mid plane density:  '   
      write(*,*)  fmm
      hh(ih) = h
      Pn(ih) = Pdonn+fwcmq
      gp(ih) = Freen
      enddo

      rewind ifc
      z = sclosew-0.5d0*dz
      do iz = istp1s,isluts
      z = z+dz
      write(ifc,*) z,fdmon(iz),fem(iz),fdgr(iz)
      enddo      
      iFc = 83
      open (iFc,file='Fh',form='formatted')
      rewind iFc
      idPc = 85
      open (idPc,file='dFh',form='formatted')
      rewind idPc
      iPc = 87
      open (iPc,file='Ph',form='formatted')
      rewind iPc
      minig = 1
      maxig = nh
      do ig = minig,maxig
      dP = Pn(ig)-Pn(maxig)
      PP = Pn(ig)
      if (dabs(dP).lt.1E-16) dP = 0.d0
      if (dabs(PP).lt.1E-16) PP = 0.d0
      dgp = gp(ig)-gp(maxig)
      if (dabs(dgp).lt.1E-16) dgp = 0.d0
      write(iPc,*) hh(ig),dP,PP
      write(iFc,*) hh(ig),dgp,gp(ig)
      enddo
      if (maxig.gt.minig) then
      ig = maxig-1
      dgp = gp(ig+1)-gp(ig)
      if (dabs(dgp).gt.1.E-16) then
      Pdnm = -(dgp)/(hh(ig+1)-hh(ig))
      else
      Pdnm = 0.d0
      endif
      do ig = minig,maxig-1
      dgp = gp(ig+1)-gp(ig)
      s = hh(ig)+0.5d0*(hh(ig+1)-hh(ig))
      if (dabs(dgp).gt.1.E-16) then
      Pdnet = -(dgp)/(hh(ig+1)-hh(ig))
      dPd = Pdnet-Pdnm
      if (dabs(dPd).lt.1E-16) dPd = 0.d0      
      if (dabs(Pdnet).lt.1E-16) Pdnet = 0.d0  
      write(idPc,*) s,dPd,Pdnet
      else
      write(idPc,*) s,0.d0,0.d0
      endif
      enddo
      endif
      STOP
      END

