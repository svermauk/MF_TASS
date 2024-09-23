! Mean force based TASS reweighting code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!                    Contribution                              !
! WTMTD: Shalini Aswasthi and others                           !
! Integration using Trapezoidal Method: Asit Pal, Subhendu Pal !
! Debugging and Interpolation: Shivani Verma                   !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


PROGRAM mean_force_2D
USE bspline_module
USE bspline_kinds_module, ONLY: wp, ip

IMPLICIT NONE
REAL*8 :: gridmin1, gridmax1, griddif1, dummy2,v, &
             gridmin2, gridmax2, griddif2,num,num2, num3, &
             gridmin3, gridmax3, griddif3,&
             den,alpha,&
             kt,kt0,ktb,bias_fact, deltaT,&
             diff_s2,diff_s,ds2,ss,hh,dum,dummy11
REAL*8 :: dmin , dmax , drange, min_dum

REAL*8, ALLOCATABLE ::cv1(:,:),cv2(:,:),ht(:,:),vbias(:,:),ct(:,:),hill(:,:),width(:,:),&
           prob(:,:),fes1(:),fes2(:),dfds(:,:),av_dfds1(:),pcons(:),&
kcons(:), norm(:),s1(:),s2(:),fes(:,:)

REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: kj_to_kcal = 0.239006 
REAL*8, PARAMETER :: deg_to_rad = 4.d0*atan(1.d0)/180.d0
REAL*8, PARAMETER :: pi=4.d0*atan(1.d0)

INTEGER :: dummy1,i,j,index1,index2,k,&
           t_min,t_max, &
           nbin1, nbin2,&
           i_mtd, i_md, i_s2, i_s1, ir, &
           mtd_max, nr, &
           narg,ios
INTEGER :: w_cv,w_hill,mtd_steps, md_steps 

CHARACTER(LEN=50) :: filename_loc
CHARACTER(LEN=50), ALLOCATABLE :: filename(:), filename_mtd(:)
CHARACTER*120 :: arg 

LOGICAL :: pmf,inpgrid,read_ct, read_vbias

! Variable declaration for spline interploation
INTEGER(ip),PARAMETER :: kx     = 4    !! x bspline order
INTEGER(ip),PARAMETER :: ky     = 4    !! y bspline order
INTEGER(ip),PARAMETER :: idx    = 0    !! [[db2val]] input
INTEGER(ip),PARAMETER :: idy    = 0    !! [[db2val]] input
INTEGER(ip),PARAMETER :: iknot  = 0    !! automatically select the knots
INTEGER(ip) :: nx_new
REAL(wp),ALLOCATABLE :: fcn_new(:,:)  !! new grid function evaluations
REAL(wp),ALLOCATABLE :: x_new(:)    !! new grid x points
REAL(wp),ALLOCATABLE :: tx(:)       !! x knots
REAL(wp),ALLOCATABLE :: ty(:)       !! y knots
REAL(wp) :: val,tru,err,errmax
INTEGER(ip) :: i_umb
INTEGER(ip) :: iflag  !! STATUS flag
INTEGER(ip) :: inbvx,inbvy,iloy
REAL(wp),DIMENSION(ky)           :: w1 !! work array
REAL(wp),DIMENSION(3*max(kx,ky)) :: w2 !! work array
!
dmin=-pi
dmax=pi
drange=2.d0*pi

OPEN(11,FILE='COLVAR_1',STATUS='unknown')
OPEN(12,FILE='HILLS_1',STATUS='unknown')
CALL get_steps(11,md_steps)
CALL get_steps(12,mtd_steps)
PRINT *, "MD steps =", md_steps
PRINT *, "MTD steps =", mtd_steps
CLOSE(11)
CLOSE(12)

kt0=300.d0
kt=3000.D0
deltaT=900.D0
t_min=1
t_max=md_steps
pmf=.FALSE.
inpgrid=.false.
narg = IARGC()
nr=1
!md_steps=0
!mtd_steps=0
w_hill=0
w_cv=0
!read_vbias=.false.
!read_ct=.false.

DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-T0').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSE IF(INDEX(arg,'-dT').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)deltaT
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
  ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_max
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)gridmin1
      CALL GETARG(i+2,arg)
      READ(arg,*)gridmax1
      CALL GETARG(i+3,arg)
      READ(arg,*)griddif1
      CALL GETARG(i+4,arg)
      READ(arg,*)gridmin2
      CALL GETARG(i+5,arg)
      READ(arg,*)gridmax2
      CALL GETARG(i+6,arg)
      READ(arg,*)griddif2
      inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_hill
  ELSE IF(INDEX(arg,'-nr').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)nr
  ELSE IF(INDEX(arg,'-mtd_steps').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)mtd_steps
  ELSE IF(INDEX(arg,'-md_steps').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)md_steps
  ELSE IF(INDEX(arg,'-read_ct').NE.0)THEN
      read_ct=.true.
  ELSE IF(INDEX(arg,'-read_vbias').NE.0)THEN
      read_vbias=.true.
  END IF
END DO

bias_fact=(kt0+deltaT)/kt0  !gamma factor in PLUMED manual

IF(nr.EQ.1)PRINT *, 'Warning! nr=1 is taken!'
IF(.NOT.inpgrid) PRINT *,"error grids have to be mentioned in the input!"

IF(md_steps.EQ.0)STOP 'md_steps cannot be zero'
IF(mtd_steps.EQ.0)STOP 'mtd_steps cannot be zero'
IF(w_hill.EQ.0)STOP 'w_hill cannot be zero'
IF(w_cv.EQ.0)STOP 'w_cv cannot be zero'

WRITE(*,'(A,F9.2)')'System Temp (K)        =',kt0
WRITE(*,'(A,F9.2)')'CV Temp (K)            =',kt
WRITE(*,'(A,F9.2)')'DeltaT Factor (K)      =',deltaT
WRITE(*,'(A,F9.2)')'Bias Factor (K)        =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
WRITE(*,'(A,I10)')'No of replicas          =',nr

nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1

WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'BLM  COORD:', &
               gridmin1,gridmax1,griddif1,nbin1
WRITE(*,'(A10,3F8.4,I10)')'MTD COORD:', &
               gridmin2,gridmax2,griddif2,nbin2

IF(nbin1.EQ.0.or.nbin2.EQ.0) STOP 'number of bins cannot be zero'
IF(nr.NE.nbin1)STOP 'nr and nbin1 have to be the same'

ALLOCATE(pcons(nr))
ALLOCATE(kcons(nr))
ALLOCATE(filename(nr))
ALLOCATE(filename_mtd(nr))
OPEN(44,FILE='replica.inp',STATUS='old')
DO ir=1,nr
  READ(44,*)pcons(ir),kcons(ir)
  kcons(ir)=kcons(ir)*kj_to_kcal
  READ(44,'(a)')filename(ir)
  READ(44,'(a)')filename_mtd(ir)
END DO

!alpha=(kt0+deltaT)/deltaT
alpha=(kt+deltaT)/deltaT   ! Use this if TEMP keyword is used with METAD keyword
                           ! in plumed.dat
kt=kb*kt                   ! kB T is calculated

!ktb=(alpha-1.d0)/kt       ! (alpha-1)/kT factor for c(t) computation

!alpha=bias_fact/(bias_fact-1.D0) ! This is gamma in the paper   
                                  ! Using kt0 as PLUMED USEs kt0 instead of kt
!bias_fact=(kt+bias_fact)/kt      ! Using correct bias_fact for c(t)
                                  ! calculation

PRINT *, "NOTE: alpha parameter         = ", alpha
PRINT *, "NOTE: kB T0 (kcal/mol)        = ", kt
PRINT *, "NOTE: (alpha-1)/kT            = ", ktb

ALLOCATE(cv1(nr,md_steps),cv2(nr,md_steps),dfds(nr,md_steps))
DO ir=1,nr
  OPEN(11,FILE=filename(ir),STATUS='old')
  DO i_md=1,md_steps
    READ(11,*)dummy11,cv1(ir,i_md),dummy11,cv2(ir,i_md),dummy11,dummy11,dummy11,dummy11,dummy11
    IF( cv1(ir,i_md) .GT.  dmax)  cv1(ir,i_md) = cv1(ir,i_md) - drange
    IF( cv1(ir,i_md) .LT.  dmin ) cv1(ir,i_md) = cv1(ir,i_md) + drange
    IF( cv2(ir,i_md) .GT.  dmax)  cv2(ir,i_md) = cv2(ir,i_md) - drange
    IF( cv2(ir,i_md) .LT. dmin ) cv2(ir,i_md) = cv2(ir,i_md) + drange
    diff_s=cv1(ir,i_md)-pcons(ir)  ! TODO pcons should be in radians
    IF (diff_s .GT. dmax ) diff_s =diff_s - drange 
    IF (diff_s .LT. dmin ) diff_s =diff_s + drange   
    dfds(ir,i_md)=-diff_s*kcons(ir) ! TODO kcons should be in kcal/mol
    ! TODO in the above, the sign is not clear; need to change if
    ! necessary
  END DO
  CLOSE(11)
END DO
      
ALLOCATE(ht(nr,mtd_steps),ct(nr,mtd_steps),hill(nr,mtd_steps),width(nr,mtd_steps))
DO ir=1,nr
  OPEN(12,FILE=filename_mtd(ir),STATUS='old')
  DO i_mtd=1,mtd_steps
    !READ(12,*) dummy11,ht(ir,i_mtd), width(ir,i_mtd),hill(ir,i_mtd)
    READ(12,*) dummy11,hill(ir,i_mtd), width(ir,i_mtd),ht(ir,i_mtd)
    ht(ir,i_mtd)=ht(ir,i_mtd)*kj_to_kcal
    IF(hill(ir,i_mtd).GT. dmax)hill(ir,i_mtd) = hill(ir,i_mtd) - drange
    IF(hill(ir,i_mtd).LT. dmin)hill(ir,i_mtd) = hill(ir,i_mtd) + drange
  END DO
  CLOSE(12)
END DO

IF(read_ct)THEN
  OPEN(77,FILE='ct.rst',FORM='unformatted')
  PRINT *, '...reading ct.rst'
  DO i_mtd=1,mtd_steps
   READ(77,IOSTAT=ios)ct(1:nr,i_mtd) 
   IF(ios.NE.0) THEN
     PRINT *, 'error reading ct.rst ; i_mtd=', i_mtd, '  nr=',nr
     STOP
   END IF
  END DO
  CLOSE(77)
ELSE
  ! Calculate c(t)
  ALLOCATE(fes1(nbin2))
  DO ir=1,nr
    CALL get_filename('ct.dat_',filename_loc,ir)
    OPEN(21,FILE=filename_loc,STATUS='unknown')
    fes1=0.d0
    DO i_mtd=1,mtd_steps
      ds2=width(ir,i_mtd)*width(ir,i_mtd)
      ss=hill(ir,i_mtd)
      hh=ht(ir,i_mtd)/alpha
      num=0.D0
      den=0.D0
      DO i_s2=1,nbin2
        diff_s=gridmin2+DFLOAT(i_s2-1)*griddif2-ss
        IF (diff_s .GT. dmax ) diff_s =diff_s - drange 
        IF (diff_s .LT. dmin ) diff_s =diff_s + drange   
        diff_s2=0.5d0*diff_s**2
        fes1(i_s2)=fes1(i_s2)+hh*DEXP(-diff_s2/ds2)
        num=num+DEXP(alpha*fes1(i_s2)/kt)
        den=den+DEXP((alpha-1.d0)*fes1(i_s2)/kt)
      END DO
      ct(ir,i_mtd)=kt*DLOG(num/den)
      WRITE(21,'(I10,F16.8)')i_mtd,ct(ir,i_mtd)
    END DO
    CLOSE(21)
    WRITE(*,'(A)')'CV values written in ',filename_loc
  END DO
  DEALLOCATE(fes1)

  OPEN(77,FILE='ct.rst',FORM='unformatted')
  DO i_mtd=1,mtd_steps
    WRITE(77)ct(1:nr,i_mtd) 
  END DO
  CLOSE(77)
END IF
PRINT *, "allocating..vbias"

ALLOCATE(vbias(nr,md_steps))
IF(read_vbias)THEN
  OPEN(77,FILE='vbias.rst',FORM='unformatted')
  DO i_md=1,md_steps
    READ(77)vbias(1:nr,i_md) 
  END DO
  CLOSE(77)
ELSE
  PRINT *, "computing vbias"
  !Calculate v(s,t)
  DO ir=1,nr
    DO i_md=1,md_steps
    !mtd_max=(i_md*w_cv/w_hill)+1
      mtd_max=((i_md-1)*w_cv/w_hill)+1 
      IF( MOD((i_md-1)*w_cv,w_hill).EQ.0)  mtd_max=mtd_max-1 
      dum=0.d0
      !DO i_mtd=1,mtd_max-1    ! Till the previous hill added till the current
                               ! time
      DO i_mtd=1,mtd_max
        ds2=width(ir,i_mtd)*width(ir,i_mtd)
        hh=ht(ir,i_mtd)/alpha
        diff_s=cv2(ir,i_md)-hill(ir,i_mtd) 
        IF (diff_s .GT. dmax ) diff_s =diff_s - drange 
        IF (diff_s .LT. dmin ) diff_s =diff_s + drange   
        diff_s2=diff_s*diff_s*0.5D0
        dum=dum+hh*DEXP(-diff_s2/ds2)
      END DO
      vbias(ir,i_md)=dum
    END DO
    PRINT *, "DOne...vbias ir=", ir
  END DO

  OPEN(77,FILE='vbias.rst',FORM='unformatted')
  DO i_md=1,md_steps
    WRITE(77)vbias(1:nr,i_md) 
  END DO
  CLOSE(77)
END IF

DEALLOCATE(hill)
PRINT *, 'hill deallocated'


! Calculate reweighted mean force < dF/ds exp(+beta [V^b (s2,t)-c(t)] ) >
! bug <totaling should be DOne within each replica>
!den=0.d0
!num=0.d0
ALLOCATE(av_dfds1(nr))
DO ir=1,nr
! bug <totaling should be DOne within each replica>
  den=0.d0
  num=0.d0
!DO i_md=1,md_steps
  DO i_md=t_min,t_max
    i_mtd=((i_md-1)*w_cv/w_hill) + 1
    ! Since the bias is only felt when md_step is greater than first mtd
    ! step, the following is added !bug
    !if(i_mtd*w_cv.gt.w_hill)THEN  !bug
     dum=vbias(ir,i_md) - ct(ir,i_mtd)
     num=num+dfds(ir,i_md)*DEXP(dum/kt)
     den=den+DEXP(dum/kt)
     !WRITE(44,*) ir, i_md, dfds(ir,i_md), dum
     !END IF !bug
  END DO
  av_dfds1(ir)=num/den
  !WRITE(43,*) ir,av_dfds1(ir),num,den
END DO
PRINT *, 'av_dfds1 is computed'
OPEN(42,FILE="av_dfds.dat")
DO ir=1,nr
  WRITE(42,*)pcons(ir), av_dfds1(ir)
END DO
CLOSE(42)

DEALLOCATE(dfds)

! Calculate projected free energy along s1
OPEN(42,FILE="free_energy_1.dat")
ALLOCATE(fes1(nr))
!ALLOCATE(fes2(nr))
fes1(1)=0.d0
num=0.d0
DO ir=1,nr-1
  dum=pcons(ir+1)-pcons(ir)
  !if (dum .GT. dmax ) dum =dum - drange  
  !if (dum .LT. dmin ) dum =dum + drange   
   num=num+dum* &
          (av_dfds1(ir+1)+av_dfds1(ir))
  fes1(ir+1)=num*0.5d0
  WRITE(42,*) pcons(ir+1),fes1(ir+1)
END DO
PRINT *, 'fes1 computed'

DEALLOCATE(pcons)
DEALLOCATE(av_dfds1)

PRINT *, 'deallocated pcons avdfds1'

ALLOCATE(prob(nbin1,nbin2))
ALLOCATE(norm(nr))
! Calculate prob (unbiased from MTD potential)       
prob=0.d0
DO ir=1,nr
  norm(ir)=0.d0 ! Normalizing each slices of probability
  den=0.d0
  DO i_md=1,md_steps
    IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
      index1 = nint((cv1(ir,i_md)-gridmin1)/griddif1) +1
      index2 = nint((cv2(ir,i_md)-gridmin2)/griddif2) +1
      IF(index1.GT.0.AND.index2.GT.0.AND. &
                index1.LE.nbin1.AND.index2.LE.nbin2)THEN
        !i_mtd=(i_md*w_cv/w_hill) + 1  ! bug
        !i_mtd=(i_md*w_cv/w_hill) 
        i_mtd=((i_md-1)*w_cv/w_hill) + 1
        !IF(i_md*w_cv.gt.w_hill)THEN   ! bug <ONLY after the first mtd
                                       ! step, the bias is experienced>
        dum=vbias(ir,i_md) - ct(ir,i_mtd)
        prob(index1,index2)=prob(index1,index2)+DEXP(dum/kt)
        den=den+DEXP(dum/kt)
        !END IF                        ! bug
      END IF
    END IF
  END DO
  norm(ir)=1.d0/(den*griddif1*griddif2)
END DO
PRINT *, 'prob computed'

DEALLOCATE(cv1)
DEALLOCATE(cv2)
DEALLOCATE(vbias)
DEALLOCATE(ht)
DEALLOCATE(ct)

ALLOCATE(s1(nbin1))
ALLOCATE(s2(nbin2))
ALLOCATE(fes(nbin1,nbin2))
DO i_s1=1,nbin1
  s1(i_s1)=DFLOAT(i_s1-1)*griddif1+gridmin1
  DO i_s2=1,nbin2
    s2(i_s2)=DFLOAT(i_s2-1)*griddif2+gridmin2
    num=prob(i_s1,i_s2)*norm(i_s1)     ! TODO: here it is assumed that i_s1
                                       ! and i_r are correctly tallied
    fes(i_s1,i_s2)=-kt*dlog(max(num,1E-32))+fes1(i_s1)
  END DO
END DO
!
! Spline Interpolation
!
nx_new=nbin2-3                         ! Change -1 to -2 or -3 if there is some
                                        ! problem in interpolation
!nx_new=nbin2

ALLOCATE(fcn_new(nx_new,nbin2))
ALLOCATE(x_new(nx_new))
ALLOCATE(tx(nbin1+kx))
ALLOCATE(ty(nbin2+ky))

! regrid:

inbvx = 1
inbvy = 1
iloy  = 1

CALL db2ink(s1,nbin1,s2,nbin2,fes,kx,ky,iknot,tx,ty,fes,iflag)
IF (iflag/=0) error STOP 'error calling db2ink'
errmax = 0.0_wp
OPEN(UNIT=10, FILE='free_energy.dat',STATUS='unknown')
DO i_umb=1,nx_new
  x_new(i_umb) = dfloat(i_umb-1)*griddif2+gridmin1   
  DO i_s2=1,nbin2
    s2(i_s2) = s2(i_s2)
    CALL db2val(x_new(i_umb),s2(i_s2),idx,idy,tx,ty,nbin1,nbin2,kx,ky,fes,val,iflag,&
               inbvx,inbvy,iloy,w1,w2)
    IF (iflag/=0) error STOP 'error calling db2val'
    fcn_new(i_umb,i_s2) = val
    !WRITE(10,'(3E16.8)') x_new(i_umb),  s2(i_s2), fcn_new(i_umb,i_s2)
  END DO
  !WRITE(10,*)
END DO

! Finding minimun.

min_dum=fcn_new(1,1)

DO i_umb=1,nx_new
  DO i_s2=1,nbin2
     IF(fcn_new(i_umb,i_s2) .LT. min_dum) min_dum=fcn_new(i_umb,i_s2)
  END DO
END DO

! Writing free energy file.

DO i_umb=1,nx_new
  DO i_s2=1,nbin2
    fcn_new(i_umb,i_s2)=fcn_new(i_umb,i_s2)-min_dum
    WRITE(10,'(3E16.8)') x_new(i_umb),s2(i_s2),fcn_new(i_umb,i_s2)
  END DO
  WRITE(10,*)
END DO

CLOSE(10)

DEALLOCATE(fcn_new)
DEALLOCATE(x_new)
DEALLOCATE(tx)
DEALLOCATE(ty)
!
! End Spline Interpolation
!
DEALLOCATE(prob)
DEALLOCATE(fes1)
DEALLOCATE(s1)
DEALLOCATE(s2)
DEALLOCATE(fes)
!DEALLOCATE(fes2)

WRITE(*,'(A)')'Free-energy and Unbiased distribution are in free_energy.dat'

END PROGRAM mean_force_2D 
!---------------------------------------------------------------------!

!------------------------** SUBROUTINES **----------------------------!

SUBROUTINE get_steps(iunit,nsteps)
IMPLICIT NONE
INTEGER iunit, nsteps
INTEGER ios1,ios2
nsteps=0
REWIND(iunit)
Read_Loop: DO
  READ(iunit,*,IOSTAT=ios1)
  IF(ios1.NE.0)EXIT Read_Loop
  ! READ(iunit,*,IOSTAT=ios2)
  ! IF(ios2.ne.0)STOP 'error! expecting even no. of steps in
  ! trajectory of cv'
  nsteps=nsteps+1
END DO Read_Loop 
REWIND(iunit)
END 
!---------------------------------------------------------------------!
SUBROUTINE check_files(iunit,dt)
IMPLICIT NONE
INTEGER iunit, dt
INTEGER ii, jj,i,ios
dt=0
i=2
REWIND(iunit)
READ(iunit,*)ii
READ(iunit,*)jj
dt=jj-ii
ii=jj
RLoop: DO 
  i=i+1
  READ(iunit,*,IOSTAT=ios)jj
  IF(ios.NE.0)EXIT RLoop
  IF(jj.NE.ii+dt)THEN
    PRINT *, '!!ERROR: Steps are not at constant stride!!'
    PRINT *, '!!       Unit No:',iunit,'!!'
    PRINT *, '!!       Line No:',i,'!!'
    PRINT *, '!! Expected stride =', dt,'!!'
    PRINT *, '!! Actual stride =', jj-ii,'!!'
    STOP
  END IF
  ii=jj
END DO RLoop
REWIND(iunit)
END 
!---------------------------------------------------------------------!
SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1,&
                          gridmin2,gridmax2,griddif2)
IMPLICIT NONE 
INTEGER :: iunit 
REAL*8  :: gridmin1, gridmax1, griddif1, & 
          gridmin2, gridmax2, griddif2
INTEGER :: ii, ios
REAL*8  :: cv1, cv2
INTEGER, PARAMETER :: Def_Grid_Size=101
REWIND(iunit)
READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
IF(ios.NE.0)STOP 'ERROR reading CV.dat'
gridmin1=cv1
gridmax1=cv1
gridmin2=cv2
gridmax2=cv2
RLoop: DO 
  READ(iunit,*,IOSTAT=ios)ii,cv1,cv2
  if(ios.NE.0)EXIT RLoop
  gridmin1=MIN(gridmin1,cv1)
  gridmin2=MIN(gridmin2,cv2)
  gridmax1=MAX(gridmax1,cv1)
  gridmax2=MAX(gridmax2,cv2)
END DO RLoop
griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
END
!---------------------------------------------------------------------!
SUBROUTINE get_filename(head,filename,ir)
IMPLICIT NONE
CHARACTER (LEN=*) :: head
CHARACTER (LEN=*) :: filename
INTEGER :: ir

IF(ir.LT.10)THEN
  WRITE(filename,'(A,I1)')trim(head),ir
ELSE IF(ir.LT.100)THEN
  WRITE(filename,'(A,I2)')trim(head),ir
ELSE
  PRINT *,'ERROR! get_filaname error: ir=',ir
  STOP 
END IF
END SUBROUTINE get_filename
