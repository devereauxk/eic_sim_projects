
*$ CREATE DT_RESNCL.FOR
*COPY DT_RESNCL
*
*===resncl=============================================================*
*
      SUBROUTINE DT_RESNCL(EPN,NLOOP,MODE)

************************************************************************
* Treatment of residual nuclei and nuclear effects.                    *
*         MODE = 1     initializations                                 *
*              = 2     treatment of final state                        *
* This version dated 16.11.95 is written by S. Roesler.                *
*                                                                      *
* Last change 05.01.2007 by S. Roesler.                                *
* 2020-03-06 M.D. Baker - allow mode =-1 for x>1 events.               *
************************************************************************

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE

      INCLUDE 'beagle.inc'
      INCLUDE 'bea_pyqm.inc'

      PARAMETER ( LINP = 5 ,
     &            LOUT = 6 ,
     &            LDAT = 9 )

      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,TINY3=1.0D-3,
     &           TINY2=1.0D-2,TINY1=1.0D-1,TINY4=1.0D-4,TINY10=1.0D-10,
     &           ONETHI=ONE/THREE)
      PARAMETER (AMUAMU = 0.93149432D0,
     &           FM2MM  = 1.0D-12,
     &           RNUCLE = 1.12D0)
      PARAMETER ( EMVGEV = 1.0                D-03 )
      PARAMETER ( AMUGEV = 0.93149432         D+00 )
      PARAMETER ( AMPRTN = 0.93827231         D+00 )
      PARAMETER ( AMNTRN = 0.93956563         D+00 )
      PARAMETER ( AMELCT = 0.51099906         D-03 )
      PARAMETER ( HLFHLF = 0.5D+00 )
      PARAMETER ( FERTHO = 14.33       D-09 )
      PARAMETER ( BEXC12 = FERTHO * 72.40715579499394D+00 )
      PARAMETER ( AMUNMU = HLFHLF * AMELCT - BEXC12 / 12.D+00 )
      PARAMETER ( AMUC12 = AMUGEV - AMUNMU )

* event history

      PARAMETER (NMXHKK=200000)

      COMMON /DTEVT1/ NHKK,NEVHKK,ISTHKK(NMXHKK),IDHKK(NMXHKK),
     &                JMOHKK(2,NMXHKK),JDAHKK(2,NMXHKK),
     &                PHKK(5,NMXHKK),VHKK(4,NMXHKK),WHKK(4,NMXHKK)

* extended event history
      COMMON /DTEVT2/ IDRES(NMXHKK),IDXRES(NMXHKK),NOBAM(NMXHKK),
     &                IDBAM(NMXHKK),IDCH(NMXHKK),NPOINT(10),
     &                IHIST(2,NMXHKK)

* particle properties (BAMJET index convention)
      CHARACTER*8  ANAME
      COMMON /DTPART/ ANAME(210),AAM(210),GA(210),TAU(210),
     &                IICH(210),IIBAR(210),K1(210),K2(210)

* nuclear potential
      LOGICAL LFERMI
      COMMON /DTNPOT/ PFERMP(2),PFERMN(2),FERMOD,
     &                EBINDP(2),EBINDN(2),EPOT(2,210),
     &                ETACOU(2),EUNBIND,ICOUL,LFERMI

* properties of interacting particles
      COMMON /DTPRTA/ IT,ITZ,IP,IPZ,IJPROJ,IBPROJ,IJTARG,IBTARG,ITMODE,
     &     ITMMOD,MODHYP,NHYPER,IDHYP(5)

* properties of photon/lepton projectiles
      COMMON /DTGPRO/ VIRT,PGAMM(4),PLEPT0(4),PLEPT1(4),PNUCL(4),IDIREC

* kinematics at lepton-gamma vertex in TRF/IRF
      COMMON /DTLGVX/ PPL0(4),PPL1(4),PPG(4),PPA(4)

* Lorentz-parameters of the current interaction
      COMMON /DTLTRA/ GACMS(2),BGCMS(2),GALAB,BGLAB,BLAB,
     &                UMO,PPCM,EPROJ,PPROJ

* treatment of residual nuclei: wounded nucleons
      COMMON /DTWOUN/ NPW,NPW0,NPCW,NTW,NTW0,NTCW,IPW(210),ITW(210)

* treatment of residual nuclei: 4-momenta
      LOGICAL LRCLPR,LRCLTA
      COMMON /DTRNU1/ PINIPR(5),PINITA(5),PRCLPR(5),PRCLTA(5),
     &                TRCLPR(5),TRCLTA(5),LRCLPR,LRCLTA

      DIMENSION PFSP(4),PSEC(4),PSEC0(4)
      DIMENSION PMOMB(5000),IDXB(5000),PMOMM(10000),IDXM(10000),
     &          IDXCOR(15000),IDXOTH(NMXHKK)

Cc...added by liang & Mark to include pythia energy loss datas
      double precision PAUX, DPF
      COMMON /PFAUX/ PAUX(4), DPF(4)
C... Local added by Mark to avoid double-counting mean-field FSI for GCF
      LOGICAL USEPOTDRG

      GOTO (1,2) MODE

*------- initializations
    1 CONTINUE

* initialize arrays for residual nuclei
      DO 10 K=1,5
         IF (K.LE.4) THEN
            PFSP(K)     = ZERO
         ENDIF
         PINIPR(K) = ZERO
         PINITA(K) = ZERO
         PRCLPR(K) = ZERO
         PRCLTA(K) = ZERO
         TRCLPR(K) = ZERO
         TRCLTA(K) = ZERO
   10 CONTINUE
      SCPOT = ONE
      NLOOP = 0

* correction of projectile 4-momentum for effective target pot.
* and Coulomb-energy (in case of hadron-nucleus interaction only)
      IF ((IP.EQ.1).AND.(IT.GT.1).AND.LFERMI) THEN
         EPNI = EPN
*   Coulomb-energy:
*     positively charged hadron - check energy for Coloumb pot.
         IF (IICH(IJPROJ).EQ.1) THEN
            THRESH = ETACOU(2)+AAM(IJPROJ)
            IF (EPNI.LE.THRESH) THEN
               WRITE(LOUT,1000)
 1000          FORMAT(/,1X,'KKINC:  WARNING!  projectile energy',
     &                ' below Coulomb threshold - event rejected',/)
               ISTHKK(1) = 1
               RETURN
            ENDIF
*     negatively charged hadron - increase energy by Coulomb energy
         ELSEIF (IICH(IJPROJ).EQ.-1) THEN
            EPNI = EPNI+ETACOU(2)
         ENDIF
         IF ((IJPROJ.EQ.1).OR.(IJPROJ.EQ.8)) THEN
*   Effective target potential
*sr 6.6. binding energy only (to avoid negative exc. energies)
C           EPNI = EPNI+EPOT(2,IJPROJ)
            EBIPOT = EBINDP(2)
            IF ((IJPROJ.NE.1).AND.(ABS(EPOT(2,IJPROJ)).GT.5.0D-3))
     &         EBIPOT = EBINDN(2)
            EPNI = EPNI+ABS(EBIPOT)
* re-initialization of DTLTRA
            DUM1 = ZERO
            DUM2 = ZERO
c...modified by liang to use uniform dtltra parameter in event generation      C...Removed by Mark - MCGENE 1-4 no longer used. MCGENE 6 same as 5.
C            IF(MCGENE.NE.5)CALL DT_LTINI(IJPROJ,IJTARG,EPNI,DUM1,DUM2,0)
         ENDIF
      ENDIF

* projectile in n-n cms
      IF ((IP.LE.1).AND.(IT.GT.1)) THEN
         PMASS1 = AAM(IJPROJ)
C* VDM assumption
C         IF (IJPROJ.EQ.7) PMASS1 = AAM(33)
         IF (IJPROJ.EQ.7) PMASS1 = AAM(IJPROJ)-SQRT(VIRT)
         IF (MCGENE.EQ.7) THEN
C Use IRF-G frame. Naive HCMS does not exist for x>1 e.g.
            PINIPR(3) = SQRT(PPG(1)**2+PPG(2)**2+PPG(3)**2)
            PINIPR(4) = PPG(4)
            PINIPR(5) = PMASS1
         ELSE
            PMASS2 = AAM(idNucBAM)
            PM1 = SIGN(PMASS1**2,PMASS1)
            PM2 = SIGN(PMASS2**2,PMASS2)
            PINIPR(4) = (UMO**2-PM2+PM1)/(TWO*UMO)
            PINIPR(5) = PMASS1
            IF (PMASS1.GT.ZERO) THEN
               PINIPR(3) = SQRT((PINIPR(4)-PINIPR(5))
     &              *(PINIPR(4)+PINIPR(5)))
            ELSE
               PINIPR(3) = SQRT(PINIPR(4)**2-PM1)
            ENDIF
         ENDIF
C         AIT  = DBLE(IT)
C         AITZ = DBLE(ITZ)
C        PINITA(5) = AIT*AMUAMU+1.0D-3*ENERGY(AIT,AITZ)
C         PINITA(5) = AIT*AMUC12+EMVGEV*EXMSAZ(AIT,AITZ,.TRUE.,IZDUM)
         IF (ITMMOD.LT.0) THEN
            PINITA(5) = AZMASS(IT,ITZ,3)
         ELSE
            PINITA(5) = AZMASS(IT,ITZ,ITMMOD)
         ENDIF

         CALL DT_LTNUC(ZERO,PINITA(5),PINITA(3),PINITA(4),3)
      ELSEIF ((IP.GT.1).AND.(IT.LE.1)) THEN
         STOP 'BeAGLE Error in DT_RESNCL: Illegal IP>1'
      ELSEIF ((IP.GT.1).AND.(IT.GT.1)) THEN
         STOP 'BeAGLE Error in DT_RESNCL: Illegal IP>1'
      ENDIF

      RETURN

*------- treatment of final state
    2 CONTINUE
*     @kdebug
      NLOOP = NLOOP+1
      IF (NLOOP.GT.1) SCPOT = 0.10D0
C     WRITE(LOUT,*) 'event ',NEVHKK,NLOOP,SCPOT

      JPW  = NPW
      JPCW = NPCW
      JTW  = NTW
      JTCW = NTCW
      DO 40 K=1,4
         PFSP(K)   = ZERO
   40 CONTINUE

      NOB = 0
      NOM = 0
      DO 900 I=NPOINT(4),NHKK
         IDXOTH(I) = -1
         IF (ISTHKK(I).EQ.1) THEN
            IF (IDBAM(I).EQ.7) GOTO 900
            IPOT = 0
            IOTHER = 0
* particle moving into forward direction
            IF (PHKK(3,I).GE.ZERO) THEN
*   most likely to be effected by projectile potential
               IPOT = 1
*     there is no projectile nucleus, try target
               IF ((IP.LE.1).OR.((IP-NPW).LE.1)) THEN
                  IPOT   = 2
                  IF (IP.GT.1) IOTHER = 1
*       there is no target nucleus --> skip
*       Include diproton or dineutron as "no target nucleus" for A=3
                  IF ((IT.LE.1).OR.((IT-NTW).LE.1) .OR.
     &              (IT.EQ.3.AND.((IT-NTW).EQ.2).AND.((ITZ-NTCW).NE.1)))
     &                 GOTO 900
               ENDIF
* particle moving into backward direction
            ELSE
*   most likely to be effected by target potential
               IPOT = 2
*     there is no target nucleus, try projectile
*     Include diproton or dineutron as "no target nucleus" for A=3
               IF ((IT.LE.1).OR.((IT-NTW).LE.1) .OR.
     &              (IT.EQ.3.AND.((IT-NTW).EQ.2).AND.((ITZ-NTCW).NE.1)))
     &              THEN
                  IPOT   = 1
                  IF (IT.GT.1) IOTHER = 1
*       there is no projectile nucleus --> skip
                  IF ((IP.LE.1).OR.((IP-NPW).LE.1)) GOTO 900
               ENDIF
            ENDIF
            IFLG = -IPOT
* nobam=3: particle is in overlap-region or neither inside proj. nor target
*      =1: particle is not in overlap-region AND is inside target (2)
*      =2: particle is not in overlap-region AND is inside projectile (1)
* flag particles which are inside the nucleus ipot but not in its
* overlap region
            IF ((NOBAM(I).NE.IPOT).AND.(NOBAM(I).LT.3)) IFLG = IPOT
            IF (IDBAM(I).NE.0) THEN
* baryons: keep all nucleons and all others where flag is set
               IF (IIBAR(IDBAM(I)).NE.0) THEN
                  IF ((IDBAM(I).EQ.1).OR.(IDBAM(I).EQ.8).OR.(IFLG.GT.0))
     &                                                              THEN
                     NOB = NOB+1
                     PMOMB(NOB) = PHKK(3,I)
                     IDXB(NOB)  = SIGN(10000000*IABS(IFLG)
     &                           +1000000*IOTHER+I,IFLG)
                  ENDIF
* mesons: keep only those mesons where flag is set
               ELSE
                  IF (IFLG.GT.0) THEN
                     NOM = NOM+1
                     PMOMM(NOM) = PHKK(3,I)
                     IDXM(NOM)  = 10000000*IFLG+1000000*IOTHER+I
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
  900 CONTINUE
*
* sort particles in the arrays according to increasing long. momentum
      CALL DT_SORT1(PMOMB,IDXB,NOB,1,NOB,1)
      CALL DT_SORT1(PMOMM,IDXM,NOM,1,NOM,1)
*
* shuffle indices into one and the same array according to the later
* sequence of correction
      NCOR = 0
      IF (IT.GT.1) THEN
         DO 910 I=1,NOB
            IF (PMOMB(I).GT.ZERO) GOTO 911
            NCOR = NCOR+1
            IDXCOR(NCOR) = IDXB(I)
  910    CONTINUE
  911    CONTINUE
         IF (IP.GT.1) THEN
            DO 912 J=1,NOB
               I = NOB+1-J
               IF (PMOMB(I).LT.ZERO) GOTO 913
               NCOR = NCOR+1
               IDXCOR(NCOR) = IDXB(I)
  912       CONTINUE
  913       CONTINUE
         ELSE
            DO 914 I=1,NOB
               IF (PMOMB(I).GT.ZERO) THEN
                  NCOR = NCOR+1
                  IDXCOR(NCOR) = IDXB(I)
               ENDIF
  914       CONTINUE
         ENDIF
      ELSE
         DO 915 J=1,NOB
            I = NOB+1-J
            NCOR = NCOR+1
            IDXCOR(NCOR) = IDXB(I)
  915    CONTINUE
      ENDIF
      DO 925 I=1,NOM
         IF (PMOMM(I).GT.ZERO) GOTO 926
         NCOR = NCOR+1
         IDXCOR(NCOR) = IDXM(I)
  925 CONTINUE
  926 CONTINUE
      DO 927 J=1,NOM
         I = NOM+1-J
         IF (PMOMM(I).LT.ZERO) GOTO 928
         NCOR = NCOR+1
         IDXCOR(NCOR) = IDXM(I)
  927 CONTINUE
  928 CONTINUE
*
*
      DO 20 ICOR=1,NCOR
         IPOT   = IABS(IDXCOR(ICOR))/10000000
         IOTHER = IABS(IDXCOR(ICOR))/1000000-IPOT*10
         I = IABS(IDXCOR(ICOR))-IPOT*10000000-IOTHER*1000000
         IDXOTH(I) = 1

         IDSEC  = IDBAM(I)

* reduction of particle momentum by corresponding nuclear potential
* (this applies only if Fermi-momenta are requested)

         IF (LFERMI) THEN

*   modification factor for nuclear potential correction,
*   it reduces the correction for particles produced with small
*   momenta in the n-n cms, i.e. far away from the original nuclei
*   and avoids a somewhat unphysical dip in the cos(theta) distribution
*   around zero caused by the cos(theta) shift in the n-n cms after
*   energy reduction in the rest frame of the colliding nuclei
*   Skip it for quasielastic where the HCMS PPCM is ill-defined
            XSCPOT = ONE
            IF (MCGENE.NE.7) THEN
               XSEC   = MAX(ABS(PHKK(3,I)/PPCM),TINY10)
               IF (XSEC.LT.ONE) XSCPOT = MIN(ONE,ONE/LOG(XSEC)**2.0D0)
            ENDIF

*   Lorentz-transformation into the rest system of the selected nucleus
            IMODE = -IPOT-1
            CALL DT_LTRANS(PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                  PSEC(1),PSEC(2),PSEC(3),PSEC(4),IDSEC,IMODE)
            PSECO  = SQRT(PSEC(1)**2+PSEC(2)**2+PSEC(3)**2)
            AMSEC  = SQRT(ABS((PSEC(4)-PSECO)*(PSEC(4)+PSECO)))
            JPMOD  = 0

            CHKLEV = TINY3
            IF ((EPROJ.GE.1.0D4).AND.(IDSEC.EQ.7)) CHKLEV = TINY1
            IF (EPROJ.GE.2.0D6) CHKLEV = 1.0D0
            IF (ABS(AMSEC-AAM(IDSEC)).GT.CHKLEV) THEN
               IF (IOULEV(3).GT.0)
     &            WRITE(LOUT,2000) I,NEVHKK,IDSEC,AMSEC,AAM(IDSEC)
 2000          FORMAT(1X,'RESNCL: inconsistent mass of particle',
     &                ' at entry ',I5,' (evt.',I8,')',/,' IDSEC: ',
     &                I4,'   AMSEC: ',E12.3,'  AAM(IDSEC): ',E12.3,/)
               GOTO 23
            ENDIF

            DO 21 K=1,4
               PSEC0(K) = PSEC(K)
   21       CONTINUE

*   the correction for nuclear potential effects is applied to as many
*   p/n as many nucleons were wounded; the momenta of other final state
*   particles are corrected only if they materialize inside the corresp.
*   nucleus (here: NOBAM = 1 part. outside proj., = 2 part. outside targ
*   = 3 part. outside proj. and targ., >=10 in overlapping region)
            IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) THEN
               IF (IPOT.EQ.1) THEN
                  IF ((JPW.GT.0).AND.(IOTHER.EQ.0)) THEN
*      this is most likely a wounded nucleon
**test
C                    RDIST = SQRT((VHKK(1,IPW(JPW))/FM2MM)**2
C    &                           +(VHKK(2,IPW(JPW))/FM2MM)**2
C    &                           +(VHKK(3,IPW(JPW))/FM2MM)**2)
C                    RAD   = RNUCLE*DBLE(IP)**ONETHI
C                    FDEN  = 1.4D0*DT_DENSIT(IP,RDIST,RAD)
C                    PSEC(4) = PSEC(4)-XSCPOT*SCPOT*FDEN*EPOT(IPOT,IDSEC)
**
                     PSEC(4) = PSEC(4)-XSCPOT*SCPOT*EPOT(IPOT,IDSEC)
                     JPW = JPW-1
                     JPMOD = 1
                  ELSE
*      correct only if part. was materialized inside nucleus
*      and if it is ouside the overlapping region
                     IF ((NOBAM(I).NE.1).AND.(NOBAM(I).LT.3)) THEN
                        PSEC(4) = PSEC(4)-XSCPOT*SCPOT*EPOT(IPOT,IDSEC)
                        JPMOD = 1
                     ENDIF
                  ENDIF
               ELSEIF (IPOT.EQ.2) THEN
C  In the case of GCF, the mean field FSI between the SRC pair and the
C  nuclear remnant is in the spectral function. We shouldn't double
C  count. IUSEPOT=1 forces double counting anyway. IUSEPOT=-1 in GCF
C  will force NO potential drag.
                  USEPOTDRG = (MCGENE.NE.7 .OR. IUSEPOT.EQ.1 .OR.
     &                 (IUSEPOT.EQ.0 .AND. ISTHKK(JMOHKK(1,I)).NE.3))
                  IF ((JTW.GT.0).AND.(IOTHER.EQ.0)) THEN
*      this is most likely a wounded nucleon
**test
C                    RDIST = SQRT((VHKK(1,ITW(JTW))/FM2MM)**2
C    &                           +(VHKK(2,ITW(JTW))/FM2MM)**2
C    &                           +(VHKK(3,ITW(JTW))/FM2MM)**2)
C                    RAD   = RNUCLE*DBLE(IT)**ONETHI
C                    FDEN  = 1.4D0*DT_DENSIT(IT,RDIST,RAD)
C                    PSEC(4) = PSEC(4)-XSCPOT*SCPOT*FDEN*EPOT(IPOT,IDSEC)
**
                     JTW = JTW-1
                     IF (USEPOTDRG) THEN
                        PSEC(4) = PSEC(4)-XSCPOT*SCPOT*EPOT(IPOT,IDSEC)
                        JPMOD = 1
                     ENDIF
                  ELSE
*      correct only if part. was materialized inside nucleus
                     IF ((NOBAM(I).NE.2).AND.(NOBAM(I).LT.3)
     &                    .AND.USEPOTDRG) THEN
                        PSEC(4) = PSEC(4)-XSCPOT*SCPOT*EPOT(IPOT,IDSEC)
                        JPMOD = 1
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               IF ((NOBAM(I).NE.IPOT).AND.(NOBAM(I).LT.3)) THEN
                  PSEC(4) = PSEC(4)-SCPOT*EPOT(IPOT,IDSEC)
                  JPMOD = 1
               ENDIF
            ENDIF

            IF (NLOOP.EQ.1) THEN
* Coulomb energy correction:
* the treatment of Coulomb potential correction is similar to the
* one for nuclear potential
               IF (IDSEC.EQ.1) THEN
                  IF ((IPOT.EQ.1).AND.(JPCW.GT.0)) THEN
                     JPCW = JPCW-1
                  ELSEIF ((IPOT.EQ.2).AND.(JTCW.GT.0)) THEN
                     JTCW = JTCW-1
                  ELSE
                     IF ((NOBAM(I).EQ.IPOT).OR.(NOBAM(I).EQ.3)) GOTO 25
                  ENDIF
               ELSE
                  IF ((NOBAM(I).EQ.IPOT).OR.(NOBAM(I).EQ.3)) GOTO 25
               ENDIF
               IF (IICH(IDSEC).EQ.1) THEN
*    pos. particles: check if they are able to escape Coulomb potential
                  IF (PSEC(4).LT.AMSEC+ETACOU(IPOT)) THEN
                     ISTHKK(I) = 14+IPOT
                     IF (ISTHKK(I).EQ.15) THEN
                        DO 26 K=1,4
                           PHKK(K,I) = PSEC0(K)
                           TRCLPR(K) = TRCLPR(K)+PSEC0(K)
   26                CONTINUE
                        IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NPW = NPW-1
                        IF (IDSEC.EQ.1) NPCW = NPCW-1
                     ELSEIF (ISTHKK(I).EQ.16) THEN
                        DO 27 K=1,4
                           PHKK(K,I) = PSEC0(K)
                           TRCLTA(K) = TRCLTA(K)+PSEC0(K)
   27                   CONTINUE
                        IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NTW = NTW-1
                        IF (IDSEC.EQ.1) NTCW = NTCW-1
                     ENDIF
                     GOTO 20
                  ENDIF
               ELSEIF (IICH(IDSEC).EQ.-1) THEN
*    neg. particles: decrease energy by Coulomb-potential
                  PSEC(4) = PSEC(4)-ETACOU(IPOT)
                  JPMOD = 1
               ENDIF
            ENDIF

   25       CONTINUE

            IF (PSEC(4).LT.AMSEC) THEN
               IF (IOULEV(6).GT.0)
     &            WRITE(LOUT,2001) I,IDSEC,PSEC(4),AMSEC
 2001          FORMAT(1X,'KKINC: particle at DTEVT1-pos. ',I5,
     &                ' is not allowed to escape nucleus',/,
     &                8X,'id : ',I3,'   reduced energy: ',E15.4,
     &                '   mass: ',E12.3)
               ISTHKK(I) = 14+IPOT
               IF (ISTHKK(I).EQ.15) THEN
                  DO 28 K=1,4
                     PHKK(K,I) = PSEC0(K)
                     TRCLPR(K) = TRCLPR(K)+PSEC0(K)
   28             CONTINUE
                  IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NPW = NPW-1
                  IF (IDSEC.EQ.1) NPCW = NPCW-1
               ELSEIF (ISTHKK(I).EQ.16) THEN
                  DO 29 K=1,4
                     PHKK(K,I) = PSEC0(K)
                     TRCLTA(K) = TRCLTA(K)+PSEC0(K)
   29             CONTINUE
                  IF ((IDSEC.EQ.1).OR.(IDSEC.EQ.8)) NTW = NTW-1
                  IF (IDSEC.EQ.1) NTCW = NTCW-1
               ENDIF
               GOTO 20
            ENDIF

            IF (JPMOD.EQ.1) THEN
               PSECN  = SQRT( (PSEC(4)-AMSEC)*(PSEC(4)+AMSEC) )
* 4-momentum after correction for nuclear potential
               DO 22 K=1,3
                  PSEC(K) = PSEC(K)*PSECN/PSECO
   22          CONTINUE

* store recoil momentum from particles escaping the nuclear potentials
               DO 30 K=1,4
                  IF (IPOT.EQ.1) THEN
                     TRCLPR(K) = TRCLPR(K)+PSEC0(K)-PSEC(K)
                  ELSEIF (IPOT.EQ.2) THEN
                     TRCLTA(K) = TRCLTA(K)+PSEC0(K)-PSEC(K)
                  ENDIF
   30          CONTINUE

* transform momentum back into n-n cms
               IMODE = IPOT+1
               CALL DT_LTRANS(PSEC(1),PSEC(2),PSEC(3),PSEC(4),
     &                     PHKK(1,I),PHKK(2,I),PHKK(3,I),PHKK(4,I),
     &                     IDSEC,IMODE)
            ENDIF

         ENDIF

   23    CONTINUE
         DO 31 K=1,4
            PFSP(K) = PFSP(K)+PHKK(K,I)
   31    CONTINUE

   20 CONTINUE

      DO 33 I=NPOINT(4),NHKK
         IF ((ISTHKK(I).EQ.1).AND.(IDXOTH(I).LT.0)) THEN
            PFSP(1) = PFSP(1)+PHKK(1,I)
            PFSP(2) = PFSP(2)+PHKK(2,I)
            PFSP(3) = PFSP(3)+PHKK(3,I)
            PFSP(4) = PFSP(4)+PHKK(4,I)
         ENDIF
   33 CONTINUE

      DO 34 K=1,5
         PRCLPR(K) = TRCLPR(K)
         PRCLTA(K) = TRCLTA(K)
   34 CONTINUE

C-TEMP-TEMP-TEMP
C      CALL DT_PYOUTEP(4)

      IF ((IP.EQ.1).AND.(IT.GT.1).AND.LFERMI) THEN
* hadron-nucleus interactions: get residual momentum from energy-
* momentum conservation
C... 2017-02-23 Added by Mark - the Pythia event is MISSING the wounded
C... nucleon's original momentum and kinetic energy. Add these back in.
C... Note: DPF(mu) = P(mu)_true - P(mu)_naive is a 4-momentum too.
C... We boost this into the HCMS in DT_PYEVENTEP now

C         DPF(1) = PXF
C         DPF(2) = PYF
C         CALL DT_LTNUC(PZF,EKF,DPF(3),DPF(4),3)

         IF (IOULEV(3).GT.1) THEN
            WRITE(*,*)
            WRITE(*,*)'DT_RESNCL E* Calculation'
            WRITE(*,*)'PQRECF,PYQ_IEG,IFERPY,FERFAC:',
     &           PQRECF,PYQ_IEG,IFERPY,FERFAC
            WRITE(*,*)'PAUX:  ',PAUX(1),PAUX(2),PAUX(3),PAUX(4)
            WRITE(*,*)'PFSP i:',PFSP(1),PFSP(2),PFSP(3),PFSP(4)
            WRITE(*,*)'PYQREC:',PYQREC(1),PYQREC(2),PYQREC(3),PYQREC(4)
            WRITE(*,*)'PINIPR:',PINIPR(1),PINIPR(2),PINIPR(3),PINIPR(4)
            WRITE(*,*)'PINITA:',PINITA(1),PINITA(2),PINITA(3),PINITA(4)
         ENDIF
C         WRITE(*,*)'PFSP(pre): ',PFSP(1),PFSP(2),PFSP(3),PFSP(4)
         DO 32 K=1,4
            PRCLPR(K) = ZERO
C...added by liang to consider quenching effect from PyQM routines
C            IF(MCGENE.EQ.5) PFSP(K)=PFSP(K)+PAUX(K)
C...DPF added by Mark to fix Pythia/DPMJet p-nonconservation
C            PFSP(K)=PFSP(K)+PAUX(K)+DPF(K)
C...Use PYQREC rather than PAUX for q-hat recoil.
C...Adding PYQREC to PFSP **SUPPRESSES** the recoil absorption
C...    so PQRECF=0 means recoil not taken
C...If pF is added back into the event, don't add it here
            IF (IFERPY.GT.1) THEN
               FERFAC=0.0D0
            ELSE
               FERFAC=1.0D0
            ENDIF
            IF (PYQ_IEG.EQ.0) THEN
               PFSP(K)=PFSP(K)+(1.0D0-PQRECF)*PYQREC(K)+FERFAC*DPF(K)
            ELSE
               PFSP(K)=PFSP(K)+FERFAC*DPF(K)
            ENDIF
            PRCLTA(K) = PINIPR(K)+PINITA(K)-PFSP(K)
   32    CONTINUE
         IF (IOULEV(3).GT.1) THEN
            WRITE(*,*)'PFSP f:',PFSP(1),PFSP(2),PFSP(3),PFSP(4)
            WRITE(*,*)'PRCLTA:',PRCLTA(1),PRCLTA(2),PRCLTA(3),PRCLTA(4)
            WRITE(*,*)
            WRITE(*,*)'Current event:'
C            CALL DT_PYOUTEP(4)
            WRITE(*,*)
         ENDIF
      ELSE
* nucleus-hadron, nucleus-nucleus: get residual momentum from
* accumulated recoil momenta of particles leaving the spectators
*   transform accumulated recoil momenta of residual nuclei into
*   n-n cms
         PZI = PRCLPR(3)
         PEI = PRCLPR(4)
         CALL DT_LTNUC(PZI,PEI,PRCLPR(3),PRCLPR(4),2)
         PZI = PRCLTA(3)
         PEI = PRCLTA(4)
         CALL DT_LTNUC(PZI,PEI,PRCLTA(3),PRCLTA(4),3)
C        IF (IP.GT.1) THEN
            PRCLPR(3) = PRCLPR(3)+PINIPR(3)
            PRCLPR(4) = PRCLPR(4)+PINIPR(4)
C        ENDIF
         IF (IT.GT.1) THEN
            PRCLTA(3) = PRCLTA(3)+PINITA(3)
            PRCLTA(4) = PRCLTA(4)+PINITA(4)
         ENDIF
      ENDIF

* check momenta of residual nuclei
      IF (LEMCCK) THEN
         CALL DT_EVTEMC(-PINIPR(1),-PINIPR(2),-PINIPR(3),-PINIPR(4),
     &               1,IDUM,IDUM)
         CALL DT_EVTEMC(-PINITA(1),-PINITA(2),-PINITA(3),-PINITA(4),
     &               2,IDUM,IDUM)
         CALL DT_EVTEMC(PRCLPR(1),PRCLPR(2),PRCLPR(3),PRCLPR(4),
     &               2,IDUM,IDUM)
         CALL DT_EVTEMC(PRCLTA(1),PRCLTA(2),PRCLTA(3),PRCLTA(4),
     &               2,IDUM,IDUM)
         CALL DT_EVTEMC(PFSP(1),PFSP(2),PFSP(3),PFSP(4),2,IDUM,IDUM)
**sr 19.12. changed to avoid output when used with phojet
C        CHKLEV = TINY3
         CHKLEV = TINY1
         CALL DT_EVTEMC(DUM,DUM,DUM,CHKLEV,-1,501,IREJ1)
C        IF ((NEVHKK.EQ.409).OR.(NEVHKK.EQ.460).OR.(NEVHKK.EQ.765))
C    &      CALL DT_EVTOUT(4)
         IF (IREJ1.GT.0) RETURN
      ENDIF

      RETURN
      END
