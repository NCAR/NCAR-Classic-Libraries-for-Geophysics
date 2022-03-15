mkdir -p ./lib
mkdir -p ./objs
( cd ./src; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f ../lib/libspherepack.a ../objs/advec.o ../objs/alf.o ../objs/divec.o ../objs/dives.o ../objs/divgc.o ../objs/divgs.o ../objs/gaqd.o ../objs/geo2math.o ../objs/gradec.o ../objs/grades.o ../objs/gradgc.o ../objs/gradgs.o ../objs/helmsph.o ../objs/hrfft.o ../objs/idivec.o ../objs/idives.o ../objs/idivgc.o ../objs/idivgs.o ../objs/idvtec.o ../objs/idvtes.o ../objs/idvtgc.o ../objs/idvtgs.o ../objs/igradec.o ../objs/igrades.o ../objs/igradgc.o ../objs/igradgs.o ../objs/ihgeod.o ../objs/isfvpec.o ../objs/isfvpes.o ../objs/isfvpgc.o ../objs/isfvpgs.o ../objs/islapec.o ../objs/islapes.o ../objs/islapgc.o ../objs/islapgs.o ../objs/ivlapec.o ../objs/ivlapes.o ../objs/ivlapgc.o ../objs/ivlapgs.o ../objs/ivrtec.o ../objs/ivrtes.o ../objs/ivrtgc.o ../objs/ivrtgs.o ../objs/lfim.o ../objs/lfin.o ../objs/lfp.o ../objs/lfpt.o ../objs/sfvpec.o ../objs/sfvpes.o ../objs/sfvpgc.o ../objs/sfvpgs.o ../objs/shaec.o ../objs/shaes.o ../objs/shagc.o ../objs/shags.o ../objs/shallow.o ../objs/shigc.o ../objs/shigs.o ../objs/shpe.o ../objs/shpg.o ../objs/shsec.o ../objs/shses.o ../objs/shsgc.o ../objs/shsgs.o ../objs/slapec.o ../objs/slapes.o ../objs/slapgc.o ../objs/slapgs.o ../objs/sphcom.o ../objs/sshifte.o ../objs/trssph.o ../objs/trvsph.o ../objs/vhaec.o ../objs/vhaes.o ../objs/vhagc.o ../objs/vhags.o ../objs/vhsec.o ../objs/vhses.o ../objs/vhsgc.o ../objs/vhsgs.o ../objs/visequ.o ../objs/visgau.o ../objs/visgeo.o ../objs/vlapec.o ../objs/vlapes.o ../objs/vlapgc.o ../objs/vlapgs.o ../objs/vrtec.o ../objs/vrtes.o ../objs/vrtgc.o ../objs/vrtgs.o ../objs/vshifte.o ../objs/vsurf.o ../objs/vtsec.o ../objs/vtses.o ../objs/vtsgc.o ../objs/vtsgs.o 
gfortran  -c advec.f -o ../objs/advec.o
gfortran  -c alf.f -o ../objs/alf.o
gfortran  -c divec.f -o ../objs/divec.o
gfortran  -c dives.f -o ../objs/dives.o
gfortran  -c divgc.f -o ../objs/divgc.o
gfortran  -c divgs.f -o ../objs/divgs.o
gfortran  -c gaqd.f -o ../objs/gaqd.o
gfortran  -c geo2math.f -o ../objs/geo2math.o
gfortran  -c gradec.f -o ../objs/gradec.o
gfortran  -c grades.f -o ../objs/grades.o
gfortran  -c gradgc.f -o ../objs/gradgc.o
gfortran  -c gradgs.f -o ../objs/gradgs.o
gfortran  -c helmsph.f -o ../objs/helmsph.o
gfortran  -c hrfft.f -o ../objs/hrfft.o
gfortran  -c idivec.f -o ../objs/idivec.o
gfortran  -c idives.f -o ../objs/idives.o
gfortran  -c idivgc.f -o ../objs/idivgc.o
gfortran  -c idivgs.f -o ../objs/idivgs.o
gfortran  -c idvtec.f -o ../objs/idvtec.o
gfortran  -c idvtes.f -o ../objs/idvtes.o
gfortran  -c idvtgc.f -o ../objs/idvtgc.o
gfortran  -c idvtgs.f -o ../objs/idvtgs.o
gfortran  -c igradec.f -o ../objs/igradec.o
gfortran  -c igrades.f -o ../objs/igrades.o
gfortran  -c igradgc.f -o ../objs/igradgc.o
gfortran  -c igradgs.f -o ../objs/igradgs.o
gfortran  -c ihgeod.f -o ../objs/ihgeod.o
gfortran  -c isfvpec.f -o ../objs/isfvpec.o
gfortran  -c isfvpes.f -o ../objs/isfvpes.o
gfortran  -c isfvpgc.f -o ../objs/isfvpgc.o
gfortran  -c isfvpgs.f -o ../objs/isfvpgs.o
gfortran  -c islapec.f -o ../objs/islapec.o
gfortran  -c islapes.f -o ../objs/islapes.o
gfortran  -c islapgc.f -o ../objs/islapgc.o
gfortran  -c islapgs.f -o ../objs/islapgs.o
gfortran  -c ivlapec.f -o ../objs/ivlapec.o
gfortran  -c ivlapes.f -o ../objs/ivlapes.o
gfortran  -c ivlapgc.f -o ../objs/ivlapgc.o
gfortran  -c ivlapgs.f -o ../objs/ivlapgs.o
gfortran  -c ivrtec.f -o ../objs/ivrtec.o
gfortran  -c ivrtes.f -o ../objs/ivrtes.o
gfortran  -c ivrtgc.f -o ../objs/ivrtgc.o
gfortran  -c ivrtgs.f -o ../objs/ivrtgs.o
gfortran  -c lfim.f -o ../objs/lfim.o
gfortran  -c lfin.f -o ../objs/lfin.o
gfortran  -c lfp.f -o ../objs/lfp.o
gfortran  -c lfpt.f -o ../objs/lfpt.o
gfortran  -c sfvpec.f -o ../objs/sfvpec.o
gfortran  -c sfvpes.f -o ../objs/sfvpes.o
gfortran  -c sfvpgc.f -o ../objs/sfvpgc.o
gfortran  -c sfvpgs.f -o ../objs/sfvpgs.o
gfortran  -c shaec.f -o ../objs/shaec.o
gfortran  -c shaes.f -o ../objs/shaes.o
gfortran  -c shagc.f -o ../objs/shagc.o
gfortran  -c shags.f -o ../objs/shags.o
gfortran  -c shallow.f -o ../objs/shallow.o
gfortran  -c shigc.f -o ../objs/shigc.o
gfortran  -c shigs.f -o ../objs/shigs.o
gfortran  -c shpe.f -o ../objs/shpe.o
gfortran  -c shpg.f -o ../objs/shpg.o
gfortran  -c shsec.f -o ../objs/shsec.o
gfortran  -c shses.f -o ../objs/shses.o
gfortran  -c shsgc.f -o ../objs/shsgc.o
gfortran  -c shsgs.f -o ../objs/shsgs.o
gfortran  -c slapec.f -o ../objs/slapec.o
gfortran  -c slapes.f -o ../objs/slapes.o
gfortran  -c slapgc.f -o ../objs/slapgc.o
gfortran  -c slapgs.f -o ../objs/slapgs.o
gfortran  -c sphcom.f -o ../objs/sphcom.o
gfortran  -c sshifte.f -o ../objs/sshifte.o
gfortran  -c trssph.f -o ../objs/trssph.o
gfortran  -c trvsph.f -o ../objs/trvsph.o
gfortran  -c vhaec.f -o ../objs/vhaec.o
gfortran  -c vhaes.f -o ../objs/vhaes.o
gfortran  -c vhagc.f -o ../objs/vhagc.o
gfortran  -c vhags.f -o ../objs/vhags.o
gfortran  -c vhsec.f -o ../objs/vhsec.o
gfortran  -c vhses.f -o ../objs/vhses.o
gfortran  -c vhsgc.f -o ../objs/vhsgc.o
gfortran  -c vhsgs.f -o ../objs/vhsgs.o
gfortran  -c visequ.f -o ../objs/visequ.o
gfortran  -c visgau.f -o ../objs/visgau.o
gfortran  -c visgeo.f -o ../objs/visgeo.o
gfortran  -c vlapec.f -o ../objs/vlapec.o
gfortran  -c vlapes.f -o ../objs/vlapes.o
gfortran  -c vlapgc.f -o ../objs/vlapgc.o
gfortran  -c vlapgs.f -o ../objs/vlapgs.o
gfortran  -c vrtec.f -o ../objs/vrtec.o
gfortran  -c vrtes.f -o ../objs/vrtes.o
gfortran  -c vrtgc.f -o ../objs/vrtgc.o
gfortran  -c vrtgs.f -o ../objs/vrtgs.o
gfortran  -c vshifte.f -o ../objs/vshifte.o
gfortran  -c vsurf.f -o ../objs/vsurf.o
gfortran  -c vtsec.f -o ../objs/vtsec.o
gfortran  -c vtses.f -o ../objs/vtses.o
gfortran  -c vtsgc.f -o ../objs/vtsgc.o
gfortran  -c vtsgs.f -o ../objs/vtsgs.o
/usr/bin/ar -rv ../lib/libspherepack.a ../objs/advec.o ../objs/alf.o ../objs/divec.o ../objs/dives.o ../objs/divgc.o ../objs/divgs.o ../objs/gaqd.o ../objs/geo2math.o ../objs/gradec.o ../objs/grades.o ../objs/gradgc.o ../objs/gradgs.o ../objs/helmsph.o ../objs/hrfft.o ../objs/idivec.o ../objs/idives.o ../objs/idivgc.o ../objs/idivgs.o ../objs/idvtec.o ../objs/idvtes.o ../objs/idvtgc.o ../objs/idvtgs.o ../objs/igradec.o ../objs/igrades.o ../objs/igradgc.o ../objs/igradgs.o ../objs/ihgeod.o ../objs/isfvpec.o ../objs/isfvpes.o ../objs/isfvpgc.o ../objs/isfvpgs.o ../objs/islapec.o ../objs/islapes.o ../objs/islapgc.o ../objs/islapgs.o ../objs/ivlapec.o ../objs/ivlapes.o ../objs/ivlapgc.o ../objs/ivlapgs.o ../objs/ivrtec.o ../objs/ivrtes.o ../objs/ivrtgc.o ../objs/ivrtgs.o ../objs/lfim.o ../objs/lfin.o ../objs/lfp.o ../objs/lfpt.o ../objs/sfvpec.o ../objs/sfvpes.o ../objs/sfvpgc.o ../objs/sfvpgs.o ../objs/shaec.o ../objs/shaes.o ../objs/shagc.o ../objs/shags.o ../objs/shallow.o ../objs/shigc.o ../objs/shigs.o ../objs/shpe.o ../objs/shpg.o ../objs/shsec.o ../objs/shses.o ../objs/shsgc.o ../objs/shsgs.o ../objs/slapec.o ../objs/slapes.o ../objs/slapgc.o ../objs/slapgs.o ../objs/sphcom.o ../objs/sshifte.o ../objs/trssph.o ../objs/trvsph.o ../objs/vhaec.o ../objs/vhaes.o ../objs/vhagc.o ../objs/vhags.o ../objs/vhsec.o ../objs/vhses.o ../objs/vhsgc.o ../objs/vhsgs.o ../objs/visequ.o ../objs/visgau.o ../objs/visgeo.o ../objs/vlapec.o ../objs/vlapes.o ../objs/vlapgc.o ../objs/vlapgs.o ../objs/vrtec.o ../objs/vrtes.o ../objs/vrtgc.o ../objs/vrtgs.o ../objs/vshifte.o ../objs/vsurf.o ../objs/vtsec.o ../objs/vtses.o ../objs/vtsgc.o ../objs/vtsgs.o 
ar: creating archive ../lib/libspherepack.a
a - ../objs/advec.o
a - ../objs/alf.o
a - ../objs/divec.o
a - ../objs/dives.o
a - ../objs/divgc.o
a - ../objs/divgs.o
a - ../objs/gaqd.o
a - ../objs/geo2math.o
a - ../objs/gradec.o
a - ../objs/grades.o
a - ../objs/gradgc.o
a - ../objs/gradgs.o
a - ../objs/helmsph.o
a - ../objs/hrfft.o
a - ../objs/idivec.o
a - ../objs/idives.o
a - ../objs/idivgc.o
a - ../objs/idivgs.o
a - ../objs/idvtec.o
a - ../objs/idvtes.o
a - ../objs/idvtgc.o
a - ../objs/idvtgs.o
a - ../objs/igradec.o
a - ../objs/igrades.o
a - ../objs/igradgc.o
a - ../objs/igradgs.o
a - ../objs/ihgeod.o
a - ../objs/isfvpec.o
a - ../objs/isfvpes.o
a - ../objs/isfvpgc.o
a - ../objs/isfvpgs.o
a - ../objs/islapec.o
a - ../objs/islapes.o
a - ../objs/islapgc.o
a - ../objs/islapgs.o
a - ../objs/ivlapec.o
a - ../objs/ivlapes.o
a - ../objs/ivlapgc.o
a - ../objs/ivlapgs.o
a - ../objs/ivrtec.o
a - ../objs/ivrtes.o
a - ../objs/ivrtgc.o
a - ../objs/ivrtgs.o
a - ../objs/lfim.o
a - ../objs/lfin.o
a - ../objs/lfp.o
a - ../objs/lfpt.o
a - ../objs/sfvpec.o
a - ../objs/sfvpes.o
a - ../objs/sfvpgc.o
a - ../objs/sfvpgs.o
a - ../objs/shaec.o
a - ../objs/shaes.o
a - ../objs/shagc.o
a - ../objs/shags.o
a - ../objs/shallow.o
a - ../objs/shigc.o
a - ../objs/shigs.o
a - ../objs/shpe.o
a - ../objs/shpg.o
a - ../objs/shsec.o
a - ../objs/shses.o
a - ../objs/shsgc.o
a - ../objs/shsgs.o
a - ../objs/slapec.o
a - ../objs/slapes.o
a - ../objs/slapgc.o
a - ../objs/slapgs.o
a - ../objs/sphcom.o
a - ../objs/sshifte.o
a - ../objs/trssph.o
a - ../objs/trvsph.o
a - ../objs/vhaec.o
a - ../objs/vhaes.o
a - ../objs/vhagc.o
a - ../objs/vhags.o
a - ../objs/vhsec.o
a - ../objs/vhses.o
a - ../objs/vhsgc.o
a - ../objs/vhsgs.o
a - ../objs/visequ.o
a - ../objs/visgau.o
a - ../objs/visgeo.o
a - ../objs/vlapec.o
a - ../objs/vlapes.o
a - ../objs/vlapgc.o
a - ../objs/vlapgs.o
a - ../objs/vrtec.o
a - ../objs/vrtes.o
a - ../objs/vrtgc.o
a - ../objs/vrtgs.o
a - ../objs/vshifte.o
a - ../objs/vsurf.o
a - ../objs/vtsec.o
a - ../objs/vtses.o
a - ../objs/vtsgc.o
a - ../objs/vtsgs.o
( cd ./test; /usr/bin/gmake clean; /usr/bin/gmake )
rm -f  testrvsph.exe tdiv.exe testrssph.exe testsshifte.exe testvshifte.exe testvtsgs.exe tgaqd.exe tgrad.exe tidvt.exe tsha.exe tshpe.exe tshpg.exe tslap.exe tvha.exe tvlap.exe tvrt.exe tvts.exe
rm -f testrvsph.exe
gfortran  testrvsph.f -o testrvsph.exe -L../lib -l spherepack
./testrvsph.exe


 EQUALLY SPACED TO GAUSSIAN GRID TRANSFER 
 trvsph input arguments: 
 intl =  0
 igride(1) = -1   igride(2) =  1
 nlone =  36   nlate =  19
 ive =  1
 igridg(1) =  2   igridg(2) =  0
 nlong = 128   nlatg =  64
 ivg =  0
 lsave =   21814   lwork =   76000   ldwork =  8321


 trvsph output: 
 ier =        0  lsvmin =   21786  lwkmin =   75097
 least squares error in u =  0.330E-06
 least squares error in v =  0.263E-06


 GAUSSIAN TO EQUALLY SPACED GRID TRANSFER 
 trvsph input arguments: 
 intl =  0
 igridg(1) =  2   igridg(2) =  0
 nlong = 128   nlatg =  64
 ivg =  0
 igride(1) = -1   igride(2) =  1
 nlone =  36   nlate =  19
 ive =  1
 lsave =   21814   lwork =   76000   ldwork =  8321


 trvsph output: 
 ier =        0  lsvmin =   21808  lwkmin =   75097
 least squares error in u =  0.438E-06
 least squares error in v =  0.308E-06
rm -f tdiv.exe
gfortran  tdiv.f -o tdiv.exe -L../lib -l spherepack
./tdiv.exe
 nlat =       15
 nlon =        9
   nt =        2
     gaqd
  ier =        0
 thtg
  0.1551E+00 0.3561E+00 0.5582E+00 0.7606E+00 0.9631E+00 0.1166E+01 0.1368E+01 0.1571E+01
  0.1773E+01 0.1976E+01 0.2178E+01 0.2381E+01 0.2583E+01 0.2786E+01 0.2986E+01
     ****
     ****
 icas =        1
     **ec
     shs 
 ierr =        0
     vha 
 ierr =        0
     shsi
 ierr =        0
     div 
 ierr =        0
 nlat =       15
 nlon =        9
 err2 =  0.68834E-06
     **ec
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    30375
 lwrk =    10125
     idvi
 ierr =        0
     idiv
 ierr =        0
 prtb
  0.3940E-07 0.4126E-08
 errv =  0.18176E-06
 errw =  0.15609E-06
     ****
     ****
 icas =        2
     **es
     shsi
 ierr =        0
     div 
 ierr =        0
 err2 =  0.68834E-06
     **es
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    30375
 lwrk =    10125
     idvi
 ierr =        0
     idiv
 ierr =        0
 prtb
  0.3940E-07 0.4126E-08
 errv =  0.18176E-06
 errw =  0.15609E-06
     ****
     ****
 icas =        3
     **gc
     shsi
 ierr =        0
     div 
 ierr =        0
 err2 =  0.74114E-06
     **gc
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    30375
 lwrk =    10125
     idvi
 ierr =        0
     idiv
 ierr =        0
 prtb
  0.2678E-07-0.7305E-08
 errv =  0.19173E-06
 errw =  0.17670E-06
     ****
     ****
 icas =        4
     **gs
     shsi
 ierr =        0
     div 
 ierr =        0
 err2 =  0.74114E-06
     **gs
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    30375
 lwrk =    10125
     idvi
 ierr =        0
     idiv
 ierr =        0
 prtb
  0.2678E-07-0.7305E-08
 errv =  0.19069E-06
 errw =  0.17670E-06
rm -f testrssph.exe
gfortran  testrssph.f -o testrssph.exe -L../lib -l spherepack
./testrssph.exe


 EQUALLY SPACED TO GAUSSIAN GRID TRANSFER 
 trssph input arguments: 
 intl =  0
 igride(1) = -1   igride(2) =  1
 nlone =  36   nlate =  19
 igridg(1) =  2   igridg(2) =  0
 nlong = 194   nlatg =  92
 lsave =   22213   lwork =   53347   ldwork =  8832


 trssph output: 
 ier =  0  lsvmin =   22213  lwkmin =   53347
 least squares error =  0.268E-06

 GAUSSIAN TO EQUALLY SPACED GRID TRANSFER 
 trssph input arguments: 
 intl =  0
 igridg(1) =  2   igridg(2) =  0
 nlong = 194   nlatg =  92
 igride(1) = -1   igride(2) =  1
 nlone =  36   nlate =  19
 lsave =   22213  lwork =   53347  ldwork =    8832


 trssph output: 
 ier =  0  lsvmin =   22213  lwkmin =   53347
 least squares error =  0.109E-05
rm -f testsshifte.exe
gfortran  testsshifte.f -o testsshifte.exe -L../lib -l spherepack
./testsshifte.exe
 sshifte arguments
 ioff =  0 nlon = 144 nlat =  72
 lsave =   608 lwork = 21024
 ier =  0
 least squares error =  0.503E-06
 sshifte arguments
 ioff =  1 nlon = 144 nlat =  72
 lsave =   608 lwork = 21024
 ier =  0
 least squares error =  0.734E-06
rm -f testvshifte.exe
gfortran  testvshifte.f -o testvshifte.exe -L../lib -l spherepack
./testvshifte.exe
 vshifte arguments
 ioff =  0 nlon = 144 nlat =  72
 lsave =   608 lwork = 21024
 ier =  0
 least squares error 
 err2u =  0.408E-06 err2v =  0.359E-06
 vshifte arguments
 ioff =  1 nlon = 144 nlat =  72
 lsave =   608 lwork = 21024
 ier =  0
 least squares error 
 err2u =  0.641E-06 err2v =  0.545E-06
rm -f testvtsgs.exe
gfortran  testvtsgs.f -o testvtsgs.exe -L../lib -l spherepack
./testvtsgs.exe
 testvtsgs: error in vt    1.756446E+02 error in wt    1.920956E+02
rm -f tgaqd.exe
gfortran  tgaqd.f -o tgaqd.exe -L../lib -l spherepack
./tgaqd.exe
 check1 on dble cn    2.00000000000D+00
 check2 on dble cn    1.27826212744D+00
 check on dble cn    1.13060255061D+00
 check1 on dble cn    2.00000000000D+00
 check2 on dble cn    1.27834264789D+00
 check on dble cn    1.13063815957D+00
 check1 on dble cn    2.00000000000D+00
 check2 on dble cn    1.27826212744D+00
 check on dble cn    1.13060255061D+00
 tdoub    2.090000E-04
 tdoub gaqd   5.099998E-05
 nlat    63  points    7.067899D-17 weights    6.592125D-15
  sumw    2.000000000000000000000000000000D+00
  sums    2.000000238418579E+00
 weights: nlat     63 irele      2 dmax    1.781561E-06 rmax    2.471997E-06
 points:  nlat     63 irele      3 dmax    3.602842E-08 rmax    7.414911E-08
 max difference in mu   2.306705E-07
 nlat    63 tsing   2.799998E-05 tdoub   5.099998E-05
rm -f tgrad.exe
gfortran  tgrad.f -o tgrad.exe -L../lib -l spherepack
./tgrad.exe
 nlat =       33
 nlon =       18
   nt =        4
     gaqd
  ier =        0
 thtg
  0.7178E-01 0.1648E+00 0.2583E+00 0.3520E+00 0.4457E+00 0.5394E+00 0.6332E+00 0.7269E+00
  0.8207E+00 0.9144E+00 0.1008E+01 0.1102E+01 0.1196E+01 0.1289E+01 0.1383E+01 0.1477E+01
  0.1571E+01 0.1665E+01 0.1758E+01 0.1852E+01 0.1946E+01 0.2040E+01 0.2133E+01 0.2227E+01
  0.2321E+01 0.2415E+01 0.2508E+01 0.2602E+01 0.2696E+01 0.2790E+01 0.2883E+01 0.2977E+01
  0.3070E+01
     ****
     ****
 icas =        1
     **ec
     shai
 ierr =        0
     sha 
 ierr =        0
     vhci
 ierr =        0
     grad
 ierr =        0
 errv =  0.63749E-06
 errw =  0.49424E-06
     vhai
 ierr =        0
     vha 
 ierr =        0
     shec
 ierr =        0
     igra
 ierr =        0
 errs =  0.72214E-07
     ****
     ****
 icas =        2
     **es
     shai
 ierr =        0
     sha 
 ierr =        0
     vhsi
 ierr =        0
     grad
 ierr =        0
 errv =  0.63749E-06
 errw =  0.49424E-06
     vhai
 ierr =        0
     vha 
 ierr =        0
     shes
 ierr =        0
     igra
 ierr =        0
 errs =  0.72214E-07
     ****
     ****
 icas =        3
     **gc
     shai
 ierr =        0
     sha 
 ierr =        0
     vhgc
 ierr =        0
     grad
 ierr =        0
 errv =  0.70052E-06
 errw =  0.55800E-06
     vhai
 ierr =        0
     vha 
 ierr =        0
     shgc
 ierr =        0
     igra
 ierr =        0
 errs =  0.82393E-07
     ****
     ****
 icas =        4
     **gs
     shai
 ierr =        0
     sha 
 ierr =        0
     vhgs
 ierr =        0
     grad
 ierr =        0
 errv =  0.69951E-06
 errw =  0.55800E-06
     vhai
 ierr =        0
     vha 
 ierr =        0
     shgs
 ierr =        0
     igra
 ierr =        0
 errs =  0.80516E-07
rm -f tidvt.exe
gfortran  tidvt.f -o tidvt.exe -L../lib -l spherepack
./tidvt.exe
 nlat =       25
 nlon =       16
   nt =        3
     gaqd
  ier =        0
 thtg
  0.9430E-01 0.2165E+00 0.3393E+00 0.4624E+00 0.5855E+00 0.7086E+00 0.8318E+00 0.9549E+00
  0.1078E+01 0.1201E+01 0.1324E+01 0.1448E+01 0.1571E+01 0.1694E+01 0.1817E+01 0.1940E+01
  0.2063E+01 0.2187E+01 0.2310E+01 0.2433E+01 0.2556E+01 0.2679E+01 0.2802E+01 0.2925E+01
  0.3047E+01
     ****
     ****
 icas =        1
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     shsi
 ierr =        0
     div 
 ierr =        0
     vrt 
 ierr =        0
 errd =  0.12991E-05
 errv =  0.95240E-06
     **ec
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =   150000
 lwrk =    50000
     idvi
 ierr =        0
     idvt
 ierr =        0
 prtd
  0.0000E+00 0.4204E-44 0.2242E-43
 prtv
  0.0000E+00 0.0000E+00 0.4204E-44
 errv =  0.14045E-06
 errw =  0.99068E-07
     ****
     ****
 icas =        2
     **es
     vhai
 ierr =        0
     vha 
 ierr =        0
     shsi
 ierr =        0
     div 
 ierr =        0
     vrt 
 ierr =        0
 errd =  0.12991E-05
 errv =  0.95240E-06
     **es
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =   150000
 lwrk =    50000
     idvi
 ierr =        0
     idvt
 ierr =        0
 prtd
  0.0000E+00 0.4204E-44 0.2242E-43
 prtv
  0.0000E+00 0.0000E+00 0.4204E-44
 errv =  0.14045E-06
 errw =  0.99068E-07
     ****
     ****
 icas =        3
     **gc
     vhai
 ierr =        0
     vha 
 ierr =        0
     shsi
 ierr =        0
     div 
 ierr =        0
     vrt 
 ierr =        0
 errd =  0.15374E-05
 errv =  0.12650E-05
     **gc
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =   150000
 lwrk =    50000
     idvi
 ierr =        0
     idvt
 ierr =        0
 prtd
  0.0000E+00 0.4204E-44 0.2242E-43
 prtv
  0.0000E+00 0.0000E+00 0.4204E-44
 errv =  0.11749E-06
 errw =  0.83465E-07
     ****
     ****
 icas =        4
     **gs
     vhai
 ierr =        0
     vha 
 ierr =        0
     shsi
 ierr =        0
     div 
 ierr =        0
     vrt 
 ierr =        0
 errd =  0.11936E-05
 errv =  0.11925E-05
     **gs
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =   150000
 lwrk =    50000
     idvi
 ierr =        0
     idvt
 ierr =        0
 prtd
  0.0000E+00 0.4204E-44 0.2242E-43
 prtv
  0.0000E+00 0.0000E+00 0.4204E-44
 errv =  0.12064E-06
 errw =  0.11935E-06
rm -f tsha.exe
gfortran  tsha.f -o tsha.exe -L../lib -l spherepack
./tsha.exe
 nlat =       15
 nlon =       18
   nt =        3
     gaqd
  ier =        0
 thtg
  0.1551E+00 0.3561E+00 0.5582E+00 0.7606E+00 0.9631E+00 0.1166E+01 0.1368E+01 0.1571E+01
  0.1773E+01 0.1976E+01 0.2178E+01 0.2381E+01 0.2583E+01 0.2786E+01 0.2986E+01
     ****
     ****
 icas =        1
     **ec
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     shs 
 ierr =        0
 err2 =  0.12874E-07
     ****
     ****
 icas =        2
     **es
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     shs 
 ierr =        0
 err2 =  0.12874E-07
     ****
     ****
 icas =        3
     **gc
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     shs 
 ierr =        0
 err2 =  0.16190E-07
     ****
     ****
 icas =        4
     **gs
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     shs 
 ierr =        0
 err2 =  0.16190E-07
rm -f tshpe.exe
gfortran  tshpe.f -o tshpe.exe -L../lib -l spherepack
./tshpe.exe

case nlat =    6 and mtrunc =    5


 error =   1.490116E-07
 tusl  =   4.200000E-05
 toe   =   1.799999E-05

case nlat =    6 and mtrunc =    4


 error =   1.788139E-07
 tusl  =   1.800002E-05
 toe   =   1.299998E-05

case nlat =    6 and mtrunc =    3


 error =   1.788139E-07
 tusl  =   1.700001E-05
 toe   =   1.199997E-05

case nlat =    7 and mtrunc =    6


 error =   2.384186E-07
 tusl  =   3.100000E-05
 toe   =   1.400005E-05

case nlat =    7 and mtrunc =    5


 error =   2.086163E-07
 tusl  =   2.500002E-05
 toe   =   1.699990E-05

case nlat =    7 and mtrunc =    4


 error =   1.490116E-07
 tusl  =   2.400007E-05
 toe   =   1.400011E-05

case nlat =    8 and mtrunc =    7


 error =   2.458692E-07
 tusl  =   4.399999E-05
 toe   =   2.799998E-05

case nlat =    8 and mtrunc =    6


 error =   2.309680E-07
 tusl  =   3.699993E-05
 toe   =   2.899999E-05

case nlat =    8 and mtrunc =    5


 error =   2.384186E-07
 tusl  =   3.700005E-05
 toe   =   2.699997E-05
rm -f tshpg.exe
gfortran  tshpg.f -o tshpg.exe -L../lib -l spherepack
./tshpg.exe

case nlat =    6 and mtrunc =    5


 error =   2.980232E-07
 tusl  =   4.099999E-05
 toe   =   1.500000E-05

case nlat =    6 and mtrunc =    4


 error =   1.490116E-07
 tusl  =   2.000004E-05
 toe   =   1.099997E-05

case nlat =    7 and mtrunc =    6


 error =   2.086163E-07
 tusl  =   3.399997E-05
 toe   =   1.399999E-05

case nlat =    7 and mtrunc =    5


 error =   2.384186E-07
 tusl  =   2.599997E-05
 toe   =   1.400005E-05

case nlat =    8 and mtrunc =    7


 error =   2.384186E-07
 tusl  =   4.500005E-05
 toe   =   2.399995E-05

case nlat =    8 and mtrunc =    6


 error =   1.788139E-07
 tusl  =   4.000007E-05
 toe   =   2.399995E-05
rm -f tslap.exe
gfortran  tslap.f -o tslap.exe -L../lib -l spherepack
./tslap.exe
 nlat =       15
 nlon =       22
   nt =        3
     gaqd
  ier =        0
 thtg
  0.1551E+00 0.3561E+00 0.5582E+00 0.7606E+00 0.9631E+00 0.1166E+01 0.1368E+01 0.1571E+01
  0.1773E+01 0.1976E+01 0.2178E+01 0.2381E+01 0.2583E+01 0.2786E+01 0.2986E+01
     ****
     ****
 icas =        1
     **ec
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     slap
 ierr =        0
 err2 =  0.10921E-04
     shsi
 ierr =        0
     isla
 ierr =        0
 ptrb
  0.0000E+00 0.0000E+00 0.0000E+00
 errs =  0.20262E-06
     ****
     ****
 icas =        2
     **es
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     slap
 ierr =        0
 err2 =  0.10921E-04
     shsi
 ierr =        0
     isla
 ierr =        0
 ptrb
  0.0000E+00 0.0000E+00 0.0000E+00
 errs =  0.20262E-06
     ****
     ****
 icas =        3
     **gc
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     slap
 ierr =        0
 err2 =  0.12798E-04
     shsi
 ierr =        0
     isla
 ierr =        0
 ptrb
  0.0000E+00 0.0000E+00 0.0000E+00
 errs =  0.22429E-06
     ****
     ****
 icas =        4
     **gs
     shai
 ierr =        0
     sha 
 ierr =        0
     shsi
 ierr =        0
     slap
 ierr =        0
 err2 =  0.12798E-04
     **gs
     shsi
 ierr =        0
     isla
 ierr =        0
 ptrb
  0.0000E+00 0.0000E+00 0.0000E+00
 errs =  0.22429E-06
rm -f tvha.exe
gfortran  tvha.f -o tvha.exe -L../lib -l spherepack
./tvha.exe
 nlat =       25
 nlon =       19
   nt =        2
     gaqd
  ier =        0
 thtg
  0.9430E-01 0.2165E+00 0.3393E+00 0.4624E+00 0.5855E+00 0.7086E+00 0.8318E+00 0.9549E+00
  0.1078E+01 0.1201E+01 0.1324E+01 0.1448E+01 0.1571E+01 0.1694E+01 0.1817E+01 0.1940E+01
  0.2063E+01 0.2187E+01 0.2310E+01 0.2433E+01 0.2556E+01 0.2679E+01 0.2802E+01 0.2925E+01
  0.3047E+01
     ****
     ****
 icas =        1
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vhs 
 ierr =        0
 errv =  0.18833E-06
 errw =  0.20575E-06
     ****
     ****
 icas =        2
     **es
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vhs 
 ierr =        0
 errv =  0.18833E-06
 errw =  0.20575E-06
     ****
     ****
 icas =        3
     **gc
     vhgi
 nlat =       25
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vhs 
 ierr =        0
 errv =  0.14764E-06
 errw =  0.21103E-06
     ****
     ****
 icas =        4
     **gs
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vhs 
 ierr =        0
 errv =  0.16451E-06
 errw =  0.21459E-06
rm -f tvlap.exe
gfortran  tvlap.f -o tvlap.exe -L../lib -l spherepack
./tvlap.exe
 nlat =       29
 nlon =       16
   nt =        1
     gaqd
  ier =        0
 thtg
  0.8152E-01 0.1871E+00 0.2933E+00 0.3997E+00 0.5061E+00 0.6125E+00 0.7190E+00 0.8255E+00
  0.9319E+00 0.1038E+01 0.1145E+01 0.1251E+01 0.1358E+01 0.1464E+01 0.1571E+01 0.1677E+01
  0.1784E+01 0.1890E+01 0.1997E+01 0.2103E+01 0.2210E+01 0.2316E+01 0.2423E+01 0.2529E+01
  0.2635E+01 0.2742E+01 0.2848E+01 0.2954E+01 0.3060E+01
     ****
     ****
 icas =        1
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vlap
 ierr =        0
 errv =  0.15730E-04
 errw =  0.17057E-04
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     ivlp
 ierr =        0
 errv =  0.99577E-07
 errw =  0.12249E-06
     ****
     ****
 icas =        2
     **es
     vhai
 ierr =        0
     vhae
 ierr =        0
     vhsi
 ierr =        0
     vlap
 ierr =        0
 errv =  0.15730E-04
 errw =  0.17057E-04
     **es
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     ivlp
 ierr =        0
 errv =  0.99577E-07
 errw =  0.12249E-06
     ****
     ****
 icas =        3
     **gc
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vlap
 ierr =        0
 errv =  0.20355E-04
 errw =  0.22735E-04
     **gc
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     ivlp
 ierr =        0
 errv =  0.82410E-07
 errw =  0.58983E-07
     ****
     ****
 icas =        4
     **gs
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     vlap
 ierr =        0
 errv =  0.16595E-04
 errw =  0.20736E-04
     **gs
     vhai
 ierr =        0
     vha 
 ierr =        0
     vhsi
 ierr =        0
     ivlp
 ierr =        0
 errv =  0.10692E-06
 errw =  0.75480E-07
rm -f tvrt.exe
gfortran  tvrt.f -o tvrt.exe -L../lib -l spherepack
./tvrt.exe
 nlat =       24
 nlon =       14
   nt =        3
     gaqd
  ier =        0
 thtg
  0.9815E-01 0.2253E+00 0.3532E+00 0.4813E+00 0.6094E+00 0.7375E+00 0.8657E+00 0.9939E+00
  0.1122E+01 0.1250E+01 0.1378E+01 0.1507E+01 0.1635E+01 0.1763E+01 0.1891E+01 0.2019E+01
  0.2148E+01 0.2276E+01 0.2404E+01 0.2532E+01 0.2660E+01 0.2788E+01 0.2916E+01 0.3043E+01
     ****
     ****
 icas =        1
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     vrti
 ierr =        0
     vrt 
 ierr =        0
 nlat =       24
 nlon =       14
 err2 =  0.95622E-06
     **ec
     vhsi
 ierr =        0
     vhs 
 ierr =        0
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    40320
 lwrk =    40320
     vhsi
 ierr =        0
     ivrt
 ierr =        0
 prtb = -0.48126E-08
 errv =  0.85499E-07
 errw =  0.56451E-07
     ****
     ****
 icas =        2
     **es
     vrti
 ierr =        0
     vrt 
 ierr =        0
 err2 =  0.95622E-06
     **es
     vhsi
 ierr =        0
     vhs 
 ierr =        0
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    40320
 lwrk =    40320
     ivti
 ierr =        0
     ivrt
 ierr =        0
 prtb = -0.48126E-08
 errv =  0.85499E-07
 errw =  0.56451E-07
     ****
     ****
 icas =        3
     **gc
     vrti
 ierr =        0
     vrt 
 ierr =        0
 err2 =  0.10029E-05
     **gc
     vhsi
 ierr =        0
     vhs 
 ierr =        0
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    40320
 lwrk =    40320
     ivti
 ierr =        0
     ivrt
 ierr =        0
 prtb =  0.34079E-09
 errv =  0.36932E-07
 errw =  0.29839E-07
     ****
     ****
 icas =        4
     **gs
     vrti
 ierr =        0
     vrt 
 ierr =        0
 err2 =  0.10029E-05
     **gs
     vhsi
 ierr =        0
     vhs 
 ierr =        0
     shai
 ierr =        0
     sha 
 ierr =        0
 lsav =    40320
 lwrk =    40320
     ivti
 ierr =        0
     ivrt
 ierr =        0
 prtb =  0.34079E-09
 errv =  0.36932E-07
 errw =  0.29839E-07
rm -f tvts.exe
gfortran  tvts.f -o tvts.exe -L../lib -l spherepack
./tvts.exe
 nlat =       25
 nlon =       19
   nt =        3
     gaqd
  ier =        0
 thtg
  0.9430E-01 0.2165E+00 0.3393E+00 0.4624E+00 0.5855E+00 0.7086E+00 0.8318E+00 0.9549E+00
  0.1078E+01 0.1201E+01 0.1324E+01 0.1448E+01 0.1571E+01 0.1694E+01 0.1817E+01 0.1940E+01
  0.2063E+01 0.2187E+01 0.2310E+01 0.2433E+01 0.2556E+01 0.2679E+01 0.2802E+01 0.2925E+01
  0.3047E+01
     ****
     ****
 icas =        1
     **ec
     vhai
 ierr =        0
     vha 
 ierr =        0
     vtsi
 ierr =        0
     vts 
 ierr =        0
 errv =  0.20209E-05
 errw =  0.23206E-05
     ****
     ****
 icas =        2
     **es
     vhai
 ierr =        0
     vha 
 ierr =        0
     vtsi
 ierr =        0
     vts 
 ierr =        0
 errv =  0.20209E-05
 errw =  0.23206E-05
     ****
     ****
 icas =        3
     **gc
     vhgi
 nlat =       25
     vhai
 ierr =        0
     vha 
 ierr =        0
     vtsi
 ierr =        0
     vts 
 ierr =        0
 errv =  0.21040E-05
 errw =  0.24129E-05
     ****
     ****
 icas =        4
     **gs
     vhai
 ierr =        0
     vha 
 ierr =        0
     vtsi
 ierr =        0
     vts 
 ierr =        0
 errv =  0.19493E-05
 errw =  0.22770E-05
