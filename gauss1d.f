C**********************************************************************C
      SUBROUTINE GAUSST(XAT,NATOMS,STR,DELTA,POT,IWORK,WORK,LENW,
     1                XTARG,NTARG,TOL,IOUT1,IOUT2,INFORM)
C**********************************************************************C
C
C     DESCRIPTION :
C
C     This subroutine computes the fast Gauss transform, that is
C     the sum of Gaussians of the form
C
C                M
C     POT(j) =  SUM  str(i) * exp[-((xat(i) - xtarg(j))^2)/delta ]
C               i=1
C
C     The source positions and targets must lie within
C     the unit interval [0,1].
C
C     TOL  is the specified precision.
C          The precision which can be achieved by a double precision
C          version of this program can not be set higher than about
C          1.0D-14. If too high a precision is requested, the
C          subroutine returns with an error code described below.
C
C     NTERMS = number of terms desired. If NTERMS = 0, the program
C          resets it to 8.
C
C     WORK is a workspace array dimensioned as REAL *8.
C          The amount of workspace required is approximately
C          ( 25 + log_2(tol) )*NATOMS.
C          If insufficient workspace is provided, the subroutine
C          returns with an error code described below.
C
C     IWORK is a workspace array dimensioned as INTEGER *4.
C
C     IOUT1,IOUT2 are logical unit numbers for printing output.
C          All reports generated during the execution of the subroutine
C          are printed on two fortran units IOUT1, IOUT2 specified
C          by the user. By setting either (or both) of these two numbers
C          to zero, the user can supress the corresponding output.
C
C     On  input :
C
C     XAT(i) = center of ith Gaussian
C     STR(i) = height of ith Gaussian.
C     NATOMS = number of Gaussians.
C     XTARG =  x-coord of ith target
C     NTARG = number of targets
C     WORK   = real work array
C     IWORK  = integer work array
C     LENW       = length of work arrays
C     TOL        = desired precision of calculation
C     NTERMS = number of terms desired
c     IOUT1, IOUT2 = options for printing outputs
C
C     On  output :
C
C     POT(I) = potential at position of ith Gaussian.
C
C     INFORM   = information returned to user
C
C              INFORM(1) is the number of levels actually used
C                 in the calculations.
C              INFORM(2) is the total number of boxes used in the
C                 calculations.
C              INFORM(3) is the amount of workspace actually used
C                 in the calculations.
C              INFORM(4) = amount of workspace which MUST BE KEPT
C                 FIXED  between successive calls to GAUSST
C              INFORM(5) is the number of terms used in the
C                 Hermite expansions to achieve desired accuracy.
C
C     IER        = error code array
C
C     IER(1) = 0  => no errors encountered
C
C     IER(1) = 4  => TOL set too high. Too many terms in Hermite
C                    expansions needed.
C                    IER(2) returns the number of terms in expansion
C                    required to satisfy TOL, which should not exceed
C                    20.
C
C     IER(1) = 8  => Insufficient workspace allotted.
C                    IER(2) returns amount of workspace needed
C                    IER(3) returns amount of workspace provided.
C
C     (local variables)
C
C     NTERMS     = number of terms used in expansions
C
C**********************************************************************C
      IMPLICIT REAL *8  (A-H,O-Z)
      INTEGER *4 NATOMS,LENW,IOUT1,IOUT2,IER(1)
      INTEGER *4 IFLAG,IFLAG2,INFORM(1),IWORK(1)
      REAL *8  XAT(1),STR(1),POT(1),TOL,WORK(1),DELTA
      REAL *8  XTARG(1)
C
C----------------------------------------------
C
C----- determine number of Hermite expansion terms needed for
C      specified precision
C
      IER(1) = 0
      CALL GET_NTERMS(TOL,DELTA,NTERMS)

      IF (NTERMS .EQ. 0) NTERMS = 8
      INFORM(5) = NTERMS
Cjw      write(6,*)' NUMBER OF TERMS IS =',NTERMS
Cjw      write(13,*)' NUMBER OF TERMS IS =',NTERMS
      IF (NTERMS .GT. 50) THEN
	 IER(1) = 4
	 IER(2) = NTERMS
      ENDIF
C
C----- determine total number of boxes and set boxsize to 64:
C
      BOXDIM = DSQRT(0.5d0*DELTA)
      NMAX = MAX(NATOMS,NTARG)
      NSIDE = INT(1.0d0/BOXDIM)+1
      KDIS = DSQRT(-2.d0*DLOG(TOL))+1
      NFMAX = NTERMS
      NLMAX = NTERMS
ccc      write(6,*) ' nmax,nterms ',nmax,nterms
ccc      write(6,*) ' enter nfmax,nlmax '
ccc      read(5,*) nfmax,nlmax

cjw      write(6,*) ' nfmax,nlmax ', nfmax,nlmax
cjw      write(13,*) ' nfmax,nlmax ', nfmax,nlmax
cjw      write(6,*)'KDIS = ',KDIS
cjw      write(13,*)'KDIS = ',KDIS
      NALLBX = NSIDE
      INFORM(2) = NALLBX
cjw      write(6,*)' TOTAL NUMBER OF BOXES IS ',NALLBX
cjw      write(13,*)' TOTAL NUMBER OF BOXES IS ',NALLBX
C
C----- allocate arrays used by algorithm from workspace.
C
      LFARXP = (NTERMS+1)
      LICNT  = NALLBX
      LICNT2 = NALLBX
      LLOC   = NALLBX
      LTOFST = NALLBX+1
      LOFFST = NALLBX+1
      LCENTR = NALLBX
      LIBOX  = NMAX
      LTRADR = NTARG
      LXATMS = NATOMS
      LCHRG2 = NATOMS
      LBEXP = (NTERMS+1)
      LTEMP = (NTERMS+1)
      LNBORS = NALLBX
C
      ICNT = 1
      ICNT2 = ICNT + LICNT
      IOFFST = ICNT2 + LICNT2
      ITOFST = IOFFST + LOFFST
      IBOX = ITOFST + LTOFST
      ITRADR = IBOX + LIBOX
      ILOC = ITRADR + LTRADR
      INBORS = ILOC + LLOC
      ITOT = ILOC + LNBORS
c
      CALL LENCHK(NSIDE,XTARG,NTARG,IWORK(ILOC),NTBOX)
c
      LLOCXP = (NTERMS+1)*NTBOX
      INFORM(4) = 0
      IFARXP =  1
      ILOCXP = IFARXP + LFARXP
      IXATMS = ILOCXP + LLOCXP
      ICHRG2 = IXATMS + LXATMS
      ICENTR = ICHRG2 + LCHRG2
      IBEXP = ICENTR + LCENTR
      ITEMP = IBEXP + LBEXP
      ITOTR = IBEXP + LTEMP
C
      IF ( ITOTR .GE. LENW ) THEN
	  IER(1) = 8
	  IER(2) = ITOTR
	  IER(3) = LENW
      ENDIF
      IF ( ITOT .GE. LENW ) THEN
	  IER(1) = 8
	  IER(2) = ITOT
	  IER(3) = LENW
      ENDIF
      INFORM(3) = ITOT
cjw      write(6,*)' TOTAL REAL STORAGE NEEDED IS = ',ITOTR
cjw      write(13,*)' TOTAL REAL STORAGE NEEDED IS = ',ITOTR
cjw      write(6,*)' TOTAL INT STORAGE NEEDED IS = ',ITOT
cjw      write(13,*)' TOTAL INT STORAGE NEEDED IS = ',ITOT
C
C-----if error encountered, return
C
      IF (IER(1) .NE. 0) THEN
	 CALL ERROUT(IER)
	 RETURN
      ENDIF

      DO 100 J = 1, NTARG
	 POT(J) = 0
100   CONTINUE
C
C----- call GAFEXP to create all expansions on grid,
C      evaluate all appropriate far field expansions,
C      and evaluate all appropriate direct interactions.
C
cc      call prinf(' calling GAFEXP, ntarg = *',ntarg,1)
      CALL GAFEXP(XAT,STR,NATOMS,XTARG,NTARG,DELTA,NSIDE,
     1   WORK(IXATMS),WORK(ICHRG2),WORK(IFARXP),
     2   WORK(ILOCXP),WORK(IBEXP),WORK(ITEMP),WORK(ICENTR),NTERMS,
     3   TOL,IWORK(IOFFST),IWORK(IBOX),IWORK(ITRADR),
     4   IWORK(ITOFST),IWORK(ICNT),IWORK(ICNT2),
     5   IWORK(INBORS),IWORK(ILOC),NTBOX,NFMAX,NLMAX,POT,KDIS)
ccc   call prinf(' finished GAFEXP *',nterms,1)
C
C-----evaluate all appropriate local expansions
C     (that is, for targets lying in boxes with more than NLMAX
C     target points).
C
      CALL GAEVAL(XTARG,POT,NTARG,WORK(ILOCXP),NTERMS,
     1   IWORK(IBOX),WORK(ICENTR),DELTA,IWORK(ILOC),NLMAX,NSIDE,
     2   IWORK(ITOFST),IWORK(ITRADR))
      RETURN
      END
C *********************************************************************C
      SUBROUTINE ERROUT(IER)
C
C     error handling messages.
C
C
C *********************************************************************C
      INTEGER *4 IER(1)
C
      IF (IER(1) .EQ. 4) THEN
	write(6,*)' ERROR IN GAUSS TRANSFORM'
	write(13,*)' ERROR IN GAUSS TRANSFORM'
      ELSE IF (IER(1) .EQ. 8) THEN
	write(6,*)' ERROR IN GAUSS TRANSFORM'
	write(13,*)' ERROR IN GAUSS TRANSFORM'
      ENDIF
      RETURN
      END
C****************************************************************
      SUBROUTINE GAFEXP(XAT,STR,NATOMS,XTARG,NTARG,
     1       DELTA,NSIDE,XATOMS,CHARG2,FFEXP,LOCEXP,B,
     2       TEMP,CENTER,NTERMS,TOL,IOFFST,IBOX,ITRADR,ITOFST,
     3       ICNT,ICNT2,NBORS,ILOC,NTBOX,NFMAX,NLMAX,POT,KDIS)
C
C     The main subroutine of fast Gauss algorithm.
C

C     STEP 1)  Sources are assigned to boxes.
C
C     Boxes containing sources are then dealt with sequentially.
C
C     STEP 2)  Far field expansions are formed for each box
C              containing more than NFMAX sources.
C
C     The neighboring boxes (a square of side 2*KDIS)
C     are then searched for targets and processed by means of
C     following decision analysis:
C
C     IF  (far field expansion has been created, i.e. there are
C          more than NFMAX sources in source box) THEN
C
C         IF (NTARG > NLMAX) .. convert far field expansion to
C                               local expansion.
C
C         IF (NTARG .LE. NLMAX) .. evaluate far field expansion
C                                  directly at target positions.
C     ELSE (there are fewer than NFMAX sources in source box
C           and no far field expansion was created)
C
C         IF (NTARG > NLMAX) .. convert each source to a
C                               local expansion.
C
C         IF (NTARG .LE. NLMAX) .. evaluate Gaussian field due to
C                                  sources directly.
C     ENDIF
C
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS,NATOMS,IOFFST(1),NTBOX
      INTEGER *4 IBOX(1),ICNT(1),ICNT2(1),NSIDE,NTARG,ILOC(1)
      INTEGER *4 NFMAX,NBORS(1),ITOFST(1),ITRADR(1)
C
      REAL *8  XAT(1),STR(1),DELTA
      REAL *8  XTARG(1),POT(1)
      REAL *8  XATOMS(1),CHARG2(1),CENTER(1),CENT
      REAL *8  FFEXP(0:NTERMS),LOCEXP(0:NTERMS,1)
      REAL *8  B(0:NTERMS),TEMP(0:NTERMS)
      DATA ZERO/0.0d0/
C
C-----Step 1) assign particles to boxes
C
      NBOXES = NSIDE
cc    call prinf(' calling assign *',nterms,0)
cc      call prinf(' ntarg =  *',ntarg,1)
      CALL ASSIGN(NSIDE,XAT,STR,XATOMS,CHARG2,IOFFST,
     1            XTARG,NTARG,IBOX,NATOMS,ICNT,ICNT2,
     2            CENTER,ITRADR,ITOFST)
cc      call prinf(' finished assign *',nterms,0)
cc      call prinf(' ntarg =  *',ntarg,1)
cc      call prinf(' ibox array is *',ibox,natoms)
C
C-----initialize local expansions to zero.
C
      DO 70 I = 1,NTBOX
	 DO 60 J = 0,NTERMS
	     LOCEXP(J,I) = ZERO
 60      CONTINUE
 70   CONTINUE
cc    call prinf(' set all locexps to zero *',nterms,0)
c
c---- process all boxes
c
      DO 800 I = 1,NBOXES
	 IOFF =   IOFFST(I)
	 NINBOX = IOFFST(I+1) - IOFF
cc	 call prinf(' ninbox = *',ninbox,1)
	 IF ( NINBOX .LE. 0 ) THEN
	    GOTO 800
	 ELSE IF ( NINBOX .LE. NFMAX ) THEN
c
c---- do not create far field expansion
c
	    CALL MKNBOR(I,NBORS,NNBORS,NSIDE,KDIS)
	    DO 650 J = 1,NNBORS
	       IADR = NBORS(J)
cc	 print *,' jth nbor is ',iadr
	       INOFF = ITOFST(IADR)
	       NINNBR = ITOFST(IADR+1) - INOFF
cc	 print *,' ninnbr is ',ninnbr
	       IF (NINNBR .LE. NLMAX) THEN
cc	  print *,'direct interaction'
c
c---- direct interaction
c
		  DO 620 K = 1,NINNBR
		     JT = ITRADR(INOFF+K)
		     XP = XTARG(JT)
		     CALL GDIR(XP,XATOMS(IOFF+1),
     1                    NINBOX,CHARG2(IOFF+1),POT(JT),DELTA)
620               CONTINUE
	       ELSE
cc        print *,'convert sources into Taylor series'
c
c---- convert each source into Taylor series and accumulate
c
		  CENT = CENTER(IADR)
cc        write(6,*)'cent = *',cent
cc        write(13,*)'cent = *',cent
		  IADLOC = ILOC(IADR)
cc        print *,' iadloc is ',iadloc
		  DO 640 K = 1,NINBOX
		     CALL GALKXP(CENT,XATOMS(IOFF+K),
     1                    CHARG2(IOFF+K),B,NTERMS,DELTA)
		     CALL ADDEXP(B,LOCEXP(0,IADLOC),NTERMS)
cc           print *,' b = ',b
cc           write(13,*)' b = ',b
640                CONTINUE
	       ENDIF
650         CONTINUE
	 ELSE
c
c---- create far field expansion
c
	    CENT = CENTER(I)
	    CALL GAMKXP(CENT,XATOMS(IOFF+1),
     1              CHARG2(IOFF+1),NINBOX,FFEXP,NTERMS,DELTA)
	    CALL MKNBOR(I,NBORS,NNBORS,NSIDE,KDIS)
	    DO 750 J = 1,NNBORS
	       IADR = NBORS(J)
	       INOFF = ITOFST(IADR)
	       NINNBR = ITOFST(IADR+1) - INOFF
	       IF (NINNBR .LE. NLMAX) THEN
cc        print *,'evaluate Hermite series'
c
c---- evaluate far field expansion
c
		  DO 720 K = 1,NINNBR
		     JT = ITRADR(INOFF+K)
		     XP = XTARG(JT)
		     CALL GAMVAL(FFEXP,NTERMS,CENT,XP,POT(JT),DELTA)
720               CONTINUE
	       ELSE
cc        print *,'convert Hermite into Taylor series'
c
C     convert far field expansion to local expansion in neighbor box
c
		  CALL GSHIFT(FFEXP,CENT,CENTER(IADR),
     1                  B,NTERMS,TEMP,DELTA)
		  IADLOC = ILOC(IADR)
		  CALL ADDEXP(B,LOCEXP(0,IADLOC),NTERMS)
	       ENDIF
750         CONTINUE
	 ENDIF
800   CONTINUE
cc      call prin2(' pot array is *',pot,natoms)
      RETURN
      END
C
C**********************************************************************
      SUBROUTINE GDIR(XP,XAT,NATOMS,STR,DPOT,DELTA)
C**********************************************************************
C
C---- increment DPOT by direct evaluation of effect of sources
C     (XAT(i),YAT(i),STR(i))    i = 1,...,NATOMS.
C
      IMPLICIT REAL *8  (A-H,O-Z)
      REAL *8 XAT(1),STR(1),DPOT,XP,DELTA
C
      DO 30 J = 1,NATOMS
	 X = XP - XAT(J)
	 DPOT = DPOT + DEXP(- (X**2)/DELTA)*STR(J)
30    CONTINUE
      RETURN
      END



C**********************************************************************
      SUBROUTINE ADDEXP(B,A,NTERMS)
C --------------------------------------------------------------------
C     INPUT :     Taylor Expansions  A and B of length NTERMS
C     OUTPUT:     A is over written by (A+B).
C***********************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4  NTERMS
      REAL *8 B(0:NTERMS),A(0:NTERMS)
      INTEGER *4 I
C
      DO 200 I = 0,NTERMS
	 A(I) = A(I) + B(I)
 200  CONTINUE
      RETURN
      END
C
C*********************************************************************C
      SUBROUTINE ASSIGN(NSIDE,XAT,STR,XATOMS,CHARG2,IOFFST,
     1              XTARG,NTARG,IBOX,NATOMS,ICNT,ICNT2,
     2              CENTER,ITRADR,ITOFST)
C*********************************************************************C
C     This subroutine assigns sources to boxes.
C
C   INPUT    XAT(i) is position of ith source.
C            STR(i) is srtength of ith source.
C
C            ICNT,ICNT2 are used as workspace arrays.
C
C   OUTPUT   IOFFST(i) indicates where in the arrays XATOMS
C            CHARG2  the listing
C            of particles contained in box i begins.
C
C            ITRADR is an ordered list of addresses of atoms.
C            The first group of addresses correspond to atoms which
C            lie in the first box, etc.
C
C            ITOFST(i) indicates where in the array ITRADR the listing
C            of particles contained in box i begins.
C
C            IBOX(i) is the address of the box containing the ith target
C
C ------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4  NATOMS,NSIDE,IOFFST(*),IBOX(NATOMS)
      INTEGER *4  ICNT(*),ICNT2(*),ITRADR(*),ITOFST(*)
      REAL *8     STR(NATOMS),CHARG2(NATOMS),CENTER(1)
      REAL *8     XAT(NATOMS)
      REAL *8     XTARG(NTARG)
      REAL *8     XATOMS(NATOMS)
C--------------------------------------------------------
      REAL *8     H
C
C----- initialize counting arrays
C
      NBOXES = NSIDE
      DO 100 J = 1,NBOXES
	 ICNT2(J) = 1
	 ICNT(J) = 0
 100  CONTINUE
      H = 1.0D0/NSIDE
C
C-----find box in which jth source lies and increment counter array
C
      DO 200 J = 1, NATOMS
	 IXH = XAT(J)/H
	 IF (IXH .GE. NSIDE) IXH = NSIDE-1
	 IF (IXH .LT. 0) IXH = 0
	 IADR = 1 + IXH
	 ICNT(IADR) = ICNT(IADR) + 1
	 IBOX(J) = IADR
 200  CONTINUE
C
C-----compute the array IOFFST.
C
      IOFFST(1) = 0
      DO 300 J = 2,NBOXES+1
	 IOFFST(J) = IOFFST(J-1) + ICNT(J-1)
 300  CONTINUE
C
C-----reorder sources in arrays XATOMS,CHARG2
C
      DO 400 J = 1,NATOMS
	 IADR = IBOX(J)
	 INDX = IOFFST(IADR) + ICNT2(IADR)
ccc      IATADR(INDX) = J
ccc      IATADR(J) = INDX
	 XATOMS(INDX) = XAT(J)
	 CHARG2(INDX) = STR(J)
	 ICNT2(IADR) = ICNT2(IADR)+1
 400  CONTINUE
C
C     ICNT will be used to count number of targets in boxes
C
C
      DO 450 J = 1,NBOXES
	 ICNT2(J) = 1
	 ICNT(J) = 0
 450  CONTINUE
      DO 500 J = 1, NTARG
	 IXH = XTARG(J)/H
	 IF (IXH .GE. NSIDE) IXH = NSIDE-1
	 IF (IXH .LT. 0) IXH = 0
	 IADR = 1 + IXH
c        print *,'iadr = ',iadr
	 ICNT(IADR) = ICNT(IADR) + 1
	 IBOX(J) = IADR
 500  CONTINUE
C
C-----compute the array IOFFST.
C
      ITOFST(1) = 0
c      print *,' itofst(1) = ',itofst(1)
      DO 550 J = 2,NBOXES+1
	 ITOFST(J) = ITOFST(J-1) + ICNT(J-1)
 550  CONTINUE
c      print *,' itofst(1) = ',itofst(1)
c      print *,'itofst = ',(itofst(i),i=1,nboxes+1)
C
C-----reorder addresses of sources in array ITRADR
C
      DO 600 J = 1,NTARG
	 IADR = IBOX(J)
	 INDX = ITOFST(IADR) + ICNT2(IADR)
	 ITRADR(INDX) = J
	 ICNT2(IADR) = ICNT2(IADR)+1
 600  CONTINUE
C
C---- create array of centers for all boxes
C
      DO 700 J = 1, NBOXES
	 ICOL = 1 + MOD(J-1,NSIDE)
	 CENTER(J) =  ICOL*H - H/2
 700  CONTINUE
ccc   call prin2(' center array is *',center,nboxes)
cc     write(6,*)' center array is ',(center(i),i=1,nboxes)
cc      write(13,*)' center array is ',(center(i),i=1,nboxes)
      RETURN
      END
C
C
C*********************************************************************C
      SUBROUTINE GAEVAL(XTARG,POT,NTARG,LOCEXP,NTERMS,IBOX,
     1                  CENTER,DELTA,ILOC,NLMAX,NSIDE,ITOFST,ITRADR)
C
C
C
C     evaluate local expansion for each particle.
C
C*****************************************************************
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4  NTERMS,NTARG,ITOFST(1),ITRADR(1)
      INTEGER *4  ILOC(1)
C
      REAL *8  XTARG(1),POT(1)
      REAL *8  CENTER(1),CENT
      REAL *8  LOCEXP(0:NTERMS,1)
      DATA ZERO/0.0d0/
C
      NBOXES = NSIDE
ccc      print *,'itofst = ',(itofst(i),i=1,nboxes+1)
      DO 200 I = 1,NBOXES
	 INOFF = ITOFST(I)
	 NINBOX = ITOFST(I+1) - INOFF
ccc      print *,'i = ',i
ccc      print *,'ninbox =',ninbox
ccc      print *,'nlmax =',nlmax
	 IF (NINBOX .LE. NLMAX) THEN
	    GOTO 200
	 ELSE
	    JADR = ILOC(I)
ccc         print *,' jadr = ',jadr
	    DO 100 K = 1,NINBOX
	       JT = ITRADR(INOFF+K)
ccc            XP = XTARG(JT)
ccc            print *,' k,jt,xp,yp =',k,jt,xp,yp
ccc            nt2 = (nterms+1)
ccc               print *,' locexp = ',(locexp(ii,jadr),ii=0,nt2)
ccc            write(13,*)' locexp = ',(locexp(ii,jadr),ii=0,nt2)
	       CALL GALVAL(LOCEXP(0,JADR),NTERMS,CENTER(I),
     1               XTARG(JT),POT(JT),DELTA)
100         CONTINUE
	 ENDIF
200   CONTINUE
      RETURN
      END
C*********************************************************************C
      SUBROUTINE GAMKXP(CENT,XATOMS,STR,NINBOX,FFEXP,
     1            NTERMS,DELTA)
C
C     This subroutine computes the far field expansion about
C     the center CENT due to the sources at locations
C     XATOMS(i) of strength STR.
C
C     INPUT
C
C     CENT  = center of the expansion
C     XATOMS(i) - x-coordinate of jth source
C     STR(j)) - strength of ith source
C     NINBOX - number of sources
C     NTERMS = number of terms in expansion
C
C     OUTPUT:
C
C     FFEXP(i,j) = (i,j)th coefficient of far field expansion
C                  due to sources
C
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NINBOX,NTERMS
      REAL *8 CENT,XATOMS(NINBOX),STR(NINBOX)
      REAL *8 FFEXP(0:NTERMS)
      REAL *8 X,XP
C
C     initialize coefficients to zero.
C
      DO 200 I=0,NTERMS
	  FFEXP(I) = 0.0D0
 200  CONTINUE
      DSQ = 1.0d0/DSQRT(DELTA)
C
C     accumulate expansion due to each source.
C
      DO 1000 I=1,NINBOX
	 X = (XATOMS(I) - CENT)*DSQ
250      CONTINUE
	 XP = STR(I)
	 DO 400 J=0,NTERMS
	    FFEXP(J) = FFEXP(J) + XP
	    XP = XP*X/(J+1)
 400     CONTINUE
1000  CONTINUE
cc    call prin2(' mpole is *',mpole,(nterms+1)**2)
      RETURN
      END
C
      SUBROUTINE GAMVAL(FFEXP,NTERMS,CENT,XPOINT,POT,DELTA)
C
C     This subroutine evaluates the far field expansion FFEXP about
C     CENT at location XP.

C     INPUT
C
C     CENT = center of the expansion
C     FFEXP(i) = ith coefficient of far field expansion
C     NTERMS = number of terms in expansion
C
C     OUTPUT:
C
C     POT = evaluated potential
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      REAL *8 CENT
      REAL *8 FFEXP(0:NTERMS)
      REAL *8 HEXPX(0:100),X,XPOINT
C
C     tabulate Hermite polynomials in HEXPX.
C
      DSQ = 1.0D0/DSQRT(DELTA)
      X = (XPOINT - CENT)*DSQ
cc    call prin2(' cent = *',cent,2)
cc    call prin2(' x = *',x,1)
cc    call prin2(' y = *',y,1)
      FACX = DEXP(-X*X)
      HEXPX(0) = 1.0D0*FACX
      HEXPX(1) = 2.0D0*X*FACX
      DO 100 I = 1,NTERMS-1
	 HEXPX(I+1) = 2.0D0*( X*HEXPX(I) - I*HEXPX(I-1))
100   CONTINUE
c     DO 150 I =0,NTERMS
c        HEXPX(I) = HEXPX(I)*FACX
c        HEXPY(I) = HEXPY(I)*FACY
150   CONTINUE
ccc   call prin2(' hexpx = *',hexpx(0),nterms+1)
ccc   call prin2(' hexpy = *',hexpy(0),nterms+1)
C
c---- evaluate expansion
C
      POTINC = 0.0D0
      DO 300 J = 0,NTERMS
	 POTINC = POTINC + HEXPX(J)*FFEXP(J)
300   CONTINUE
      POT = POT + POTINC
      RETURN
      END
C
C*******************************************************************
      SUBROUTINE GALKXP(CENT,XP,CHARGE,B,NTERMS,DELTA)
C
C     This subroutine converts the source at (XP) of strength CHARGE
C     into the Taylor expansion B about CENT.
C
C     INPUT:
C
C     XP  center of Gaussian
C     CHARGE strength of Gaussian
C     CENT = center of the box for Taylor series expansion
C     NTERMS = number of terms in expansion
C     DELTA = variance
C
C     OUTPUT:
C
C     B = Taylor expansion coefficients
C
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      REAL *8 CENT,XP,CHARGE,DELTA
      REAL *8 B(0:NTERMS),X
      REAL *8 HEXPX(0:100)
      REAL *8 FACX,PROD,FAC(0:50)
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
      DSQ = 1.0d0/DSQRT(delta)
      X = (CENT - XP)*DSQ
      FACX = DEXP(-X*X)*CHARGE
      HEXPX(0) = 1.0D0*FACX
      HEXPX(1) = 2.0D0*X*FACX
      DO 100 I = 1,2*NTERMS-1
	 HEXPX(I+1) = 2.0D0*( X*HEXPX(I) - I*HEXPX(I-1))
100   CONTINUE
cc    call prin2(' hexpx = *',hexpx(0),2*nterms+1)
cc    call prin2(' hexpy = *',hexpy(0),2*nterms+1)
      FAC(0) = 1.0d0
      DO 150 I = 1,NTERMS
	 FAC(I) = FAC(I-1)/I
150   CONTINUE
c
c---- compute local expansion
c
      DO 1600 L = 0,NTERMS
	 B(L) = HEXPX(L)*FAC(L)*(-1)**(L)
1600  CONTINUE
      RETURN
      END
C
C**********************************************************************
      SUBROUTINE GSHIFT(FFEXP,CENT1,CENT2,LOCAL,NTERMS,TEMP,DELTA)
C
C     This subroutine converts the far field expansion FFEXP about
C     the center CENT1 to a local expansion LOCAL about CENT2.
C
C     INPUT
C
C     CENT1 = center of the original expansion MPOEL
C     FFEXP = far field expansion coefficients
C     NTERMS = number of terms in expansion
C     CENT2 = center about which local expansion is to be formed
C
C     OUTPUT:
C
C     LOCAL(i) = ith coefficient of shifted local expansion
C                  due to sources
C
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      REAL *8 CENT1,CENT2
      REAL *8 FFEXP(0:NTERMS),LOCAL(0:NTERMS)
      REAL *8 HEXPX(0:100),X,XPOINT
      REAL *8 FACX,SUM,FAC(0:50),TEMP(0:NTERMS)
C
C     tabulate Hermite polynomials in HEXPX and HEXPY.
C
ccc      if (2.ne.3) return
      DSQ = 1.0d0/DSQRT(delta)
      X = (CENT2 - CENT1)*DSQ
cc    call prin2(' inside GSHIFT, shift component x  = *',x,1)
cc    call prin2(' inside GSHIFT, shift component y  = *',y,1)
      FACX = DEXP(-X*X)
      HEXPX(0) = 1.0D0*FACX
      HEXPX(1) = 2.0D0*X*FACX
      DO 100 I = 1,2*NTERMS-1
	 HEXPX(I+1) = 2.0D0*( X*HEXPX(I) - I*HEXPX(I-1))
100   CONTINUE
cc    call prin2(' hexpx = *',hexpx(0),2*nterms+1)
cc    call prin2(' hexpy = *',hexpy(0),2*nterms+1)
      FAC(0) = 1.0d0
      DO 150 I = 1,NTERMS
	 FAC(I) = FAC(I-1)/I
150   CONTINUE
c
c---- compute local expansion
c
      DO 1600 L = 0,NTERMS
	    SUM = 0.0D0
	    DO 1400 J= NTERMS,0,-1
		  SUM = SUM + FFEXP(J)*HEXPX(L+J)
1400        CONTINUE
	    LOCAL(L) = SUM*FAC(L)*(-1)**(L)
1600  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE GALVAL(LOCAL,NTERMS,CENT,XPOINT,POT,DELTA)
C
C     This subroutine evaluates the local expansion LOCAL about
C     CENT at location XPOINT.
C
C     INPUT
C
C     CENT = center of the expansion
C     LOCAL(i) = ith coefficient of local expansion
C     NTERMS = number of terms in expansion
C
C     OUTPUT:
C
C     POT = evaluated potential
C
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4 NTERMS
      REAL *8 CENT,POT,PARSUM,POTINC
      REAL *8 LOCAL(0:NTERMS),X
C
      DSQ = 1.0d0/DSQRT(DELTA)
      X = (XPOINT - CENT)*DSQ
ccc   call prin2(' inside GALVAL cent = *',cent,2)
ccc   call prin2(' inside GALVAL x = *',x,1)
ccc   call prin2(' inside GALVAL y = *',y,1)
ccc   call prin2(' inside GALVAL, local(j,0) = *',local(0,0),nterms+1)
C
c---- evaluate expansion by Horner's rule
c
      POTINC = LOCAL(NTERMS)
      DO 300 J = NTERMS,1,-1
	 POTINC = X*POTINC + LOCAL(J-1)
300   CONTINUE
      POT = POT + POTINC
      RETURN
      END
C
      SUBROUTINE LENCHK(NSIDE,XTARG,NTARG,ILOC,NTBOX)
C*********************************************************************C
C     This subroutine determines the amount of workspace to be reserved
C     for local expansions and sets up indirect addresses for
C     relevant boxes...
C
C     INPUT    XTARG(i),YTARG(i) is position of ith target.
C              NSIDE = number of boxes on a side.
C
C     OUTPUT   NTBOX is number of boxes containing targets.
C              ILOC(J)  is pointer to location in memory of local
C              expansion for box J.
C ------------------------------------------------------
      IMPLICIT REAL *8 (A-H,O-Z)
      INTEGER *4  NTBOX,NSIDE,NTARG
      INTEGER *4  ILOC(*)
      REAL *8     XTARG(NTARG)
C--------------------------------------------------------
      REAL *8     H
C
C----- initialize counting arrays
C
      NBOXES = NSIDE
      H = 1.0D0/NSIDE
      DO 100 I = 1,NBOXES
	 ILOC(I) = 0
100   CONTINUE
C
C-----find box in which jth target lies, increment counter
C     and compute ILOC(J)
C
      NTBOX = 0
      DO 500 J = 1, NTARG
	 IXH = XTARG(J)/H
	 IF (IXH .GE. NSIDE) IXH = NSIDE-1
	 IF (IXH .LT. 0) IXH = 0
	 IADR = 1 + IXH
	 IF (ILOC(IADR) .EQ. 0) THEN
	    NTBOX = NTBOX+1
	    ILOC(IADR) = NTBOX
	 ENDIF
 500  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE MKNBOR(IBOX,NBORS,NNBORS,NSIDE,KDIS)
C
C     creates list of neighbors of box IBOX.
C
C     NSIDE = number of boxes on a side (sqrt of total number)
C     KDIS = max number of neighbors to consider in each direction.
C     IBOX = box number
C     NBORS = array of neighbor addresses
C     NNBORS = number of neighbors to consider (avoid nonexistent
C              neighbors at boundary)
C
C*********************************************************
      INTEGER *4 IBOX,NBORS(1),NNBORS
C
      NNBORS = 0
C
C---- compute actual COL and ROW number of IBOX
C
      ICOLB = IBOX
ccc   print *,' icolb =',icolb
ccc   print *,' irowb =',irowb
C
C---- compute actual range (a subset of -KDIS,KDIS X -KDIS,KDIS)
C
      IMIN = MAX(ICOLB-KDIS,1)
      IMAX = MIN(ICOLB+KDIS,NSIDE)
ccc   print *,' imin =',imin
ccc   print *,' imax =',imax
ccc   print *,' jmin =',jmin
ccc   print *,' jmax =',jmax
C
C----- determine NBORS list (nearest and second nearest shells).
C
      DO 600 I= IMIN,IMAX
	    NNBORS = NNBORS+1
	    NBORS(NNBORS) = I
600   CONTINUE
      RETURN
      END
C


c----------------------------------------------------------

        subroutine get_nterms(tol,delta,nterms)
          implicit real*8 (a-h,o-z)
          real*8 tol,delta
          integer nterms

          TOLER=TOL
          BOXDIM = DSQRT(0.5d0*DELTA)
          NSIDE = INT(1.0d0/BOXDIM)
          boxside=1.d0/(float(nside)*sqrt(2.0*delta))
          IP=1
          Z= 4.d0*1.09*boxside**2/(1.-boxside)
645       CONTINUE
          IP=IP+1
          Z=boxside*Z/sqrt(FLOAT(IP))
          IF(Z.GT.TOLER)GOTO 645
cc          WRITE(6,*) 'Z= ',Z,' TOLER = ',TOLER,' Q= ',Q,' IP = ',IP
          NTERMS=IP
ccc          print *,' we recommend NTERMS =  ',NTERMS
ccc          write(13,*) ' we recommend NTERMS =  ',NTERMS
ccc          print *,' enter NTERMS'
ccc          read *, NTERMS
c          print *, ' NTERMS = ',NTERMS
        end subroutine



