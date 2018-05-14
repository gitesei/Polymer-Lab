      program platem
      implicit double precision (a-h,o-z)
      include 'dftpol.id.inc'
      dimension c(maxmon,0:maxel),fps(0:maxel),
     *chA(0:maxel),chB(0:maxel),
     *Pn(10001),hh(10001),gp(10001),aPn(10001),agp(10001)
      ifc = 38
      ins = 49
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
      ddtol = 0.000002d0
      bdm = 0.1
      dz = 0.05
      closew = 0.d0
      T = 1.d0
      write(*,*) 'IDEAL CHAINS, NON-ADSORBING SURFACES'
      open (ifc,file='densitydistribution',form='formatted')    
      open (ins,file='input',form='formatted')
      rewind ifc
      rewind ins
      read(ins,*) 
      read(ins,*) nmon
      read(ins,*) 
      read(ins,*)  hstart,dh,nh
      sclosew = closew
      rrT = 1.d0/T
      rdz = 1.d0/dz
      twopidz = 0.5d0*dz
      irdz = int(rdz+0.001d0)

      istart = int(closew/dz+0.01d0)
      istp1 = istart+1 
      ists = int(sclosew/dz+0.01d0)
      istp1s = ists+1 
      ism = int(1.d0/dz+0.01d0)
      isms = int(q1/dz+0.01d0)
      inw = ism+int(closew/dz+0.01d0)
      inws = isms+int(sclosew/dz+0.01d0)
      pie = pi/8.d0
      dzpie = pie*dz
      rnmon = real(nmon)
      rrnmon = 1.d0/rnmon      
      bdpol = bdm/rnmon
      Pb = bdpol
      chempp = dlog(bdpol)
      hwdist = 0.d0
      write(*,*) 'hwdist = ',hwdist
c      cmpar = dfloat(ncmpar)
      scalem = chempp/(2.d0*rnmon)
      emscale = 2.d0*scalem
      write(*,*) 'monomer density (bdm) =',bdm
      write(*,*) 'polymer density (bdpol) =',bdpol
      write(*,*) 'no. of monomers/polymer = ',nmon
      write(*,*) 'polymer chemical pot. (betamu) = ',chempp
      write(*,*) 'total bulk pressure = ',Pb
      write(*,*) 'grid spacing (dz) = ',dz

      h = hstart-dh
      do ih = 1,nh         
      h = h+dh
      write(*,*) 
      write(*,*) 'SEPARATION: ',h
      islut = int((h-closew)/dz+0.01d0)  
      isluts = int((h-sclosew)/dz+0.01d0)          
      nfack = int(h/dz+0.01d0)
      imitt = nfack/2

      do 36 iz = istp1,imitt
      ehbclam(iz) = dexp(-0.5d0*Vexm(iz)+scalem)
 36   ebelam(iz) = dexp(-Vexs(iz)+emscale)
      jz = imitt+1
      do 502 iz = imitt+1,islut
      jz = jz-1
      ebelam(iz) = ebelam(jz)
 502  ehbclam(iz) = ehbclam(jz)      

      nAB = 1
      do  245 iz = istp1,inw
      sume = 0.d0
      do 345 jz = istp1,iz+ism-1
 345  sume = ebelam(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*ebelam(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 245  c(nmon-1,iz) = tuu
  
      do 445 iz = inw+1,imitt
      sume = 0.5d0*ebelam(iz-ism)
      do 545 jz = iz-ism+1,iz+ism-1
 545  sume = ebelam(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*ebelam(iz+ism)+sume)*twopidz
c      chA(iz) = tuu
      chA(iz) = tuu*ehbclam(iz)
 445  c(nmon-1,iz) = tuu
      jz = imitt+1
      do 4445 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
 4445 c(nmon-1,iz) = tuu

      k = nmon-1
      do 745 mmm = 2,nmon-2
      k = k-1
      nAB = nAB+1
      if (mod(nAB,2).eq.0) then

      do  845 iz = istp1,inw
      sume = 0.d0
      do 945 jz = istp1,iz+ism-1
 945  sume = chA(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
 845  c(k,iz) = tuu
      do 1045 iz = inw+1,imitt
      sume = 0.5d0*chA(iz-ism)
      do 1145 jz = iz-ism+1,iz+ism-1
 1145 sume = chA(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      chB(iz) = tuu*ehbclam(iz)
 1045 c(k,iz) = tuu 
      jz = imitt+1
      do 7045 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chB(jz)
      chB(iz) = tuu
 7045 c(k,iz) = tuu
      else

      do  895 iz = istp1,inw
      sume = 0.d0
      do 995 jz = istp1,iz+ism-1
 995  sume = chB(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 895  c(k,iz) = tuu
      do 9045 iz = inw+1,imitt
      sume = 0.5d0*chB(iz-ism)
      do 9145 jz = iz-ism+1,iz+ism-1
 9145 sume = chB(jz)+sume
      tuu = ehbclam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      chA(iz) = tuu*ehbclam(iz)
 9045 c(k,iz) = tuu 
      jz = imitt+1
      do 9345 iz = imitt+1,imitt+ism
      jz = jz-1
      tuu = chA(jz)
      chA(iz) = tuu
 9345 c(k,iz) = tuu
      endif

 745  continue

      nAB = nAB+1
      if (mod(nAB,2).eq.0) then

      do  1245 iz = istp1,inw
      sume = 0.d0
      do 1345 jz = istp1,iz+ism-1
 1345 sume = chA(jz)+sume
 1245 c(1,iz) = ebelam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      do 1445 iz = inw+1,imitt
      sume = 0.5d0*chA(iz-ism)
      do 1545 jz = iz-ism+1,iz+ism-1
 1545 sume = chA(jz)+sume
 1445 c(1,iz) = ebelam(iz)*(0.5d0*chA(iz+ism)+sume)*twopidz
      else

      do  6245 iz = istp1,inw
      sume = 0.d0
      do 6345 jz = istp1,iz+ism-1
 6345 sume = chB(jz)+sume
 6245 c(1,iz) = ebelam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      do 6445 iz = inw+1,imitt
      sume = 0.5d0*chB(iz-ism)
      do 6545 jz = iz-ism+1,iz+ism-1
 6545 sume = chB(jz)+sume
 6445 c(1,iz) = ebelam(iz)*(0.5d0*chB(iz+ism)+sume)*twopidz
      endif

      do 9 i = istp1,imitt
      dumsum = 0.d0 
      do 10 k = 2,nmon-1
 10   dumsum = c(k,i)*c(nmon+1-k,i)+dumsum
      fem(i) = 2.d0*c(1,i)
 9    fdmon(i) = dumsum+fem(i)

      jz = imitt+1
      do 11 iz = imitt+1,isluts
      jz = jz-1
      do k = 1,nmon-1
      c(k,iz) = c(k,jz)
      enddo
      ehbclam(iz) = ehbclam(jz)
      ebelam(iz) = ebelam(jz)
      fdmon(iz) = fdmon(jz)
 11   fem(iz) = fem(jz)
      do 4747 iz = istp1,islut-ism
      sfps = 0.d0
      do 5757 k = 2,nmon-2
 5757 sfps = ehbclam(iz)*c(nmon+1-k,iz)*
     *ehbclam(iz+ism)*c(k+1,iz+ism)+sfps
      sfps = ehbclam(iz)*c(2,iz)*ebelam(iz+ism)+sfps
 4747 fps(iz) = ebelam(iz)*ehbclam(iz+ism)*c(2,iz+ism)+sfps

      sumfdm = 0.d0
      Freen = 0.d0
      chFreen = 0.d0
      Fpol = 0.d0
      do 500 iz = istp1,imitt
      bclamb = 2.d0*(dlog(ehbclam(iz))-scalem)
      belamb = dlog(ebelam(iz))-emscale
      fdm = fdmon(iz)      
      fde = fem(iz)      
      fdc = (fdmon(iz)-fde)
      sumfdm = fdm+sumfdm
      chFreen = chFreen-fdm*rrnmon
 500  Freen = Freen+fdc*bclamb+fde*belamb-fdm*rrnmon
      Freen = 2.d0*Freen*dz
      chFreen = 2.d0*chFreen*dz
c      write(*,*) 'Freen = ',Freen
c      write(*,*) 'chFreen = ',chFreen
      avfdm = sumfdm/real(islut-istart)
      gammam = 0.5d0*(2.*sumfdm*dz-bdm*h)

      Fs2m = 0.d0
      Pmacross = 0.d0
      Psacross = 0.d0
      avfdm = 0.d0
      write(*,*)
      write(*,*) 'grand potential '
      write(*,'(1e25.14)') Freen+Pb*h
      write(*,*)
      write(*,*) 'grand potential (check): '
      write(*,'(1e25.14)') chFreen+Pb*h
      Pmacross = -2.d0*Pmacross*dz
      fp1S = fdmon(istp1)
      fn1S = fdmon(istp1+2)
      c0Skv = fdmon(istp1+1)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fwcmq = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
      if (fwcmq.lt.0.d0) fwcmq = 0.d0
C      ** EXTRAPOLATION BY USE OF LINEAR EXPRESSION **
      fwcml = fp1S+0.5d0*(fp1S-c0Skv)
      if (fwcml.lt.0.d0) fwcmq = 0.d0
      Fs2m = fwcmq-Fs2m*dz
      Pnet = Fs2m-bdpol
      write(*,*) 
      write(*,*) 'monomer contact density - quad. extr.:',fwcmq
      write(*,*) 
      write(*,*)'net pressure, at wall: '
      write(*,'(1e25.14)')  Pnet
      write(*,*) 

      fp1S = fdmon(imitt+1)
      fn1S = fdmon(imitt+3)
      c0Skv = fdmon(imitt+2)
      c1Skv = (fp1S-fn1S)*0.5d0
      c2Skv = (fp1S+fn1S-2.d0*c0Skv)*0.5d0
C      ** EXTRAPOLATION BY USE OF QUADRATIC EXPRESSION **
      fmm = (c0Skv+c1Skv*1.5d0+c2Skv*2.25d0)
c      write(*,*) 'monomer mid plane density - quad. extr.:',fmm
      Pidmid = fmm
      write(*,*) 'ideal part of mid plane pressure:'
      write(*,'(1e25.14)') Pidmid
      write(*,*) 
      write(*,*) 'Pmacross:'
      write(*,'(1e25.14)') Pmacross
      sumfps = 0.d0
      do 688 i = imitt-ism+1,imitt
 688  sumfps = fps(i)+sumfps
      sumfps = sumfps*twopidz
      Pbr = -2.d0*sumfps
      write(*,*) 
      write(*,*) 'total bridging pressure: '
      write(*,'(1e25.14)') Pbr
      write(*,*) 
      write(*,*) 'net pressure, at mid plane: '
      write(*,'(1e25.14)') Pbr+Pidmid+Pmacross-bdpol
c      avfdm = avfdm/dfloat((islut+1-istp1))
c      write(*,*) 
c      write(*,*)'average monomer density, <n(z)> = '
c      write(*,'(1e25.14)')  avfdm
c      write(*,*) 
c      write(*,*)'0.5*<n(z)>*h '
c      write(*,'(1e25.14)')  0.5d0*h*avfdm
      write(*,*) 
      hh(ih) = h
      Pn(ih) = Pnet
      gp(ih) = Freen+Pb*h
      agp(ih) = chFreen+Pb*h
      aPn(ih) = Pbr+Pidmid+Pmacross-bdpol
      enddo

      rewind ifc
      z = sclosew-0.5d0*dz
      do iz = istp1s,isluts
      z = z+dz
      write(ifc,*) z,fdmon(iz),fem(iz)
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
      write(iPc,*) hh(ig),Pn(ig),aPn(ig)
      write(iFc,*) hh(ig),gp(ig)-gp(maxig),gp(ig)
      enddo
      if (maxig.gt.minig) then
      ig = maxig-1
      Pdnm = -(gp(ig+1)-gp(ig))/(hh(ig+1)-hh(ig))
      aPdnm = -(agp(ig+1)-agp(ig))/(hh(ig+1)-hh(ig))
      do 2 ig = minig,maxig-1
      Pdnet = -(gp(ig+1)-gp(ig))/(hh(ig+1)-hh(ig))
      aPdnet = -(agp(ig+1)-agp(ig))/(hh(ig+1)-hh(ig))
      s = hh(ig)+0.5d0*(hh(ig+1)-hh(ig))
2     write(idPc,*) s,Pdnet,aPdnet
      endif

      STOP
      END
