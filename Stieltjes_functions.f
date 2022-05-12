!__________________________________________________________________________________
!
! Author: Yury Grabovsky https://www.math.temple.edu/~yury
! Last modified: Feb 14, 2021
! This code is distributed under the GNU LGPL license.
! It is the FORTRAN implementation of the algorithm described in
! https://www.math.temple.edu/~yury/Stieltjes.pdf
!
!__________________________________________________________________________________

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c The subroutine StCaprini is the driver routine. The actual main program only 
c acquires data, calls this driver routine, and writes the results to text files.
c
c StCaprini(z,w,zeta,nz,n_zeta,NMC,npd,mar,mnd,p,W_extr,WMC,
c           status,g,s0,rtk,sigma,ns,tC,C,ntC)
c z: nz points in the complex UHP
c w: on input, the approximate values of f(z) for a Stieltjes function f(.)
c on output, the "alternative data" representing a different realization of
c noise in the data for which our results are certifiably optimal
c zeta: n_zeta points in the UHP where the vlues of f(zeta) are wanted
c NMC: number of Monte-Carlo tries to quantify the uncertainty in f(zeta)
c npd: number of points per decade for discretization of [0,+infinity)
c mar: margin for the above discretization
c mnd: maximal number of decades so that we can reserve enough computer memory
c p: the closest feasible point to the "alternative data" w
c W_extr: computed values of f(zeta)
c WMC: rows of this matrix are possible W_extr under different realization of noise
c status: status=false means don't show the Caprini function because 
c either the optimality could not be achieved or the data is exact (no noise)
c In addition to computed values f(zeta) we also compute the spectral representation
c of f(z)=g-s0/z+sum(sigma(k)/(rtk(k)-z),k=1,ns)
c tC,C are lists of length ntC representing the Caprini function C(j)=C(tC(j))

      subroutine StCaprini(z,w,zeta,nz,n_zeta,NMC,npd,mar,mnd,
     $                     p,W_extr,WMC,status,
     $                     g,s0,rtk,sigma,ns,tC,C,ntC)
! Inputs:z,w,zeta,nz,n_zeta,NMC,npd,mar,mnd
! Outputs: w,NMC,p,W_extr,WMC,status,g,s0,rtk,sigma,ns,tC,C,ntC
! Both inputs and outputs: w,NMC
      implicit none
      complex*16 p(nz),p_direct(nz),p_indirect(nz),W_extr(n_zeta)
      complex*16 WMC(NMC,n_zeta),z(nz),w(nz),zeta(n_zeta),w_orig(nz)
      real*8 sigma(nz),rtk(nz),g,s0,mar,sgm,Caprinitmax
      real*8 tmin,tmax,distCn,znorm,dist_direct,dist_indirect
      real*8 tC(mnd*npd),C(mnd*npd)!mnd=max number of decades to reserve enough memory
      integer nz,n_zeta,NMC,npd,ns,ntC,numlm,lmc,mnd,i,nt4C,size_t4C
      integer pad,start_bar0,start_bar,start_ar,ind_bar,ind
      logical status,Nval,Pval
      real*8, allocatable::t4C(:)
      real*8, allocatable::tk(:)
      real*8, allocatable::sk(:)
      real*8, allocatable::ALLC(:)
      real*8, allocatable::lmns(:)
      real*8, allocatable::Cmns(:)

      call Valcheck(z,w,nz,Nval,Pval)
      if (Nval .and. Pval) then 
         NMC=0 !If we have no noise we cannot run Monte-Carlo nor use the Caprini method
         call Stieltjes(zeta,z,w,n_zeta,nz,NMC,W_extr,WMC,p,npd,mar)
         status=.false. ! do not compute the Caprini function
         write(6,*) 'Data is exact'
      else
         call Stieltjes(zeta,z,w,n_zeta,nz,NMC,W_extr,WMC,p,npd,mar) 
         tmin=minval(imagpart(z))
         tmax=Caprinitmax(z,w,p,nz)
         nt4C=size_t4C(tmin,tmax,mar,npd)
         allocate(t4C(nt4C))
         allocate(ALLC(nt4C))
         allocate(tk(nt4C))
         allocate(sk(nt4C))
         
         w_orig=w

         call t4Caprini(tmin,tmax,mar,npd,nt4C,t4C)
         call FullCaprinifix(t4C,z,w,p,nz,nt4C,tk,lmc,sk,ALLC,status)
      end if

! sk is of length 1 more than size(tk)
! If status=false because we were not able to achieve optimality, 
! the program will report it.
! tk,sk give the spectral measure indirectly from the Caprini function
! We can also get the spectral measure directly (see the paper)
      call Spectr(z,p,nz,npd,mar,g,s0,sigma,rtk,ns) 

! We'll see which spectral measure does better and use the best one

      sgm=distCn(w_orig,p,nz)/znorm(w_orig,nz)
      if (distCn(w_orig,w,nz)/znorm(w_orig,nz)>sgm) then
         status=.false.
          write(6,*) 'Optimality was not reached by small perturbations
     $of original data'
      end if
      write(6,*) 'Caprini plot status is',status

      if (status) then
! Compute values of C(t) at the local minima
         if (tk(1).eq.0.0) then
            numlm=lmc-1
            allocate(lmns(numlm))
            lmns=tk(2:lmc)
         else
            numlm=lmc
            allocate(lmns(numlm))
            lmns=tk(1:lmc)
         end if
         allocate(Cmns(numlm))

         call CapriniC(lmns,z,w,p,nz,numlm,Cmns)

! Merge local minima of C(t) with other valued of C(t)
         pad=200 ! We don't want a loc min to be at the very edge of the graph of C(t)
         start_bar0=min(npd,count(t4C<lmns(1))-pad) ! There can be very small locmin
         start_bar0=max(start_bar0,2) ! t4C(1)=0, which we must exclude
         ntC=nt4C-start_bar0+1+numlm ! the size of the meaningful part of C and tC
         start_bar=start_bar0
         start_ar=1
      do i=1,numlm
         ind_bar=count(t4C<lmns(i)) ! there are ind_bar-old_ind_bar elements to insert
         ind=ind_bar-start_bar0+i
         tC(start_ar:ind)=t4C(start_bar:ind_bar)
         C(start_ar:ind)=ALLC(start_bar:ind_bar) 
         tC(ind+1)=lmns(i)
         C(ind+1)=Cmns(i)
         start_ar=ind+2
         start_bar=ind_bar+1
      end do
      tC(start_ar:ntC)=t4C(start_bar:nt4C)
      C(start_ar:ntC)=ALLC(start_bar:nt4C)

! We have two ways of computing spectral measure. Let us decide which is better
! Spectral measure from the Caprini function
         p_indirect=sk(1)
         do i=1,lmc
            p_indirect=p_indirect+sk(i+1)/(tk(i)-z)
         end do
         dist_indirect=distCn(p_indirect,p,nz)
! Spectral measure computed directly
         p_direct=g-s0/z
         do i=1,ns
            p_direct=p_direct+sigma(i)/(rtk(i)-z)
         end do
         dist_direct=distCn(p_direct,p,nz)
!         write(6,*) '|p-p_dr|=',dist_direct,'|p-p_indr|=',dist_indirect
! It looks like in most cases directly computed spectral measure is better,

         if (dist_direct>dist_indirect) then ! Caprini SM is better
            write(6,*) 'Spectral measure comes from Caprini'
            ns=numlm
            g=sk(1)
            if (tk(1)==0) then
               s0=sk(2)
               if (numlm>0) then
                  rtk(1:numlm)=lmns
                  sigma(1:numlm)=sk(3:lmc+1)
               end if
            else
               s0=0.0
               if (numlm>0) then
                  rtk(1:numlm)=lmns
                  sigma(1:numlm)=sk(2:lmc+1)
               end if
            end if
         end if ! If direct SM is better we don't change ns,g,s0,rtk and sigma
      end if ! If status=false, only the direct SM is computed
      end subroutine StCaprini

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Given the noisy measurements w(i) of a Stieltjes function f(z) at z(i)
! we compute the Stieltjes function such that |w-f(z)| is as small as possible
! What is returned is the list of values f(zeta(k)) for a given list zeta of
! points in the upper half-plane. The outputs are: p(j)=f(z(j)), F(k)=f(zeta(k)).
! The presense of noise make the reconstruction uncertain.
! We also compute the degree of uncertaintly of the extrapolation.
! We use the discrepancy between p and w to estimate the noise level.
! We then use that noise level to polute p. We do it Nr times.
! Each of the Nr rows of W_extr contain F corresponding to each Monte-Carlo try.
! Additional inputs are npd=number of points per decade for discretization
! of [0,+infinity) and mar govern this discretization. Here these parameters
! are carried passively to be passed to other subroutines down the stream.

      subroutine Stieltjes(zeta,z,w,nzeta,nz,Nr,F,W_extr,p,npd,mar)
! Inputs: zeta,z,w,nzeta,nz,Nr,npd,mar
! Outputs: F,W_extr,p
      implicit none
      complex*16 zeta(nzeta),z(nz),w(nz),p(nz),F(nzeta),wk(nz)
      complex*16 Wexp(Nr,nz),Wall(Nr,nz+nzeta),W_extr(Nr,nzeta)
      real*8 distCn,sigma,nrm,new_nrm,mar
      integer nzeta,nz,Nr,npd,k,kbest
      logical Nval,Pval

      call Valcheck(z,w,nz,Nval,Pval)
      if (Nval .and. Pval) then 
         Nr=0 !If we have no noise we cannot run Monte-Carlo
         p=w
      else
! Compute our first candidate in V(z) for closest point to w
         call myt_proj(z,w,nz,p,npd,mar) 
      end if
      
! Uncertainty quantification
      if (Nr>0) then   
         sigma=distCn(w,p,nz)/sqrt(2.0*nz-1.0) ! estimate noise level
! Run Monte-Carlo
         call NPmonte_carlo(zeta,z,p,sigma,nz,nzeta,Nr,Wall,npd,mar)
         Wexp=Wall(:,nzeta+1:nzeta+nz) ! Wexp(s*,last nz elements)=p for best s*
         W_extr=Wall(:,1:nzeta) ! W(s*,first nzeta elements)=F
      end if

! Choose best data fit among all Monte-Carlo runs
      nrm=distCn(p,w,nz)
      kbest=0
      do k=1,Nr
         wk=Wexp(k,:)
         new_nrm=distCn(wk,w,nz)
        if (new_nrm<nrm) then
           kbest=k
           nrm=new_nrm
        end if
      end do
      if (kbest>0) then
         p=Wexp(kbest,:)
      else
      end if
    
! compute Stieljes lsq interpolant using best fit
      call Sinterp(z,p,zeta,nz,nzeta,F,npd,mar)
      end subroutine Stieltjes
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compute the Caprini function C(t). Here t is a list of values of t of length nt, 
c z,w,p are complex vectors of length nz. C(t)=Re[sum((p(i)-w(i))/(t-conjg(z(i)))]

      subroutine CapriniC(t,z,w,p,nz,nt,C)
! Inputs: t,z,w,p,nz,nt
! Outputs: C
      implicit none
      real*8 C(nt),t(nt)
      complex*16 z(nz),w(nz),p(nz)
      integer nz,nt,i
! we assume that array t(nt) is huge, z,w,p (nz) are relatively short
! we assume that array operations are faster in FORTRAN compared to loops
      C=realpart((p(1)-w(1))/(t-conjg(z(1))))
      do i=2,nz
         C=C+realpart((p(i)-w(i))/(t-conjg(z(i))))
      end do
      return
      end subroutine CapriniC

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compute all points of local minima of the Caprini function C(t)
c They are returned as first lmc components of the vector tk
c The subroutine also computes a vector of values C(t(j)) stored in C
c Inputs are vector t of length nt, discretizing [0,+infty)
c z,w,p of length nz, where w are measured values of f(z) for some
c Stieltjes function f(zeta). p is in the feasible set V(z)
c and is supposed to be as close to w as possible
c We use BrentRoots module taken from J-P Moreau, www.jpmoreau.fr

      subroutine CapriniCT(t,z,w,p,nz,nt,C,tk,lmc)
! Inputs: t,z,w,p,nz,nt
! Outputs: C,tk,lmc (only tk(1:lmc) contains the output)
      implicit none
      real*8 t(nt),C(nt),AC(nt),tol,VAT,TolX,locmin,eps
      real*8 tk(nt),t3(3),Cpr(3),BrentRoots,Minimum,x1,x2
      real*8, allocatable::DC(:)
      complex*16 z(nz),w(nz),p(nz)
      integer nz,nt,i,inds(nt),er,ni,MIT,k,j,num_locmins,nt3,lmc
      integer RootBracketed
      integer, allocatable::locmin_ind(:)
      integer, allocatable::allinds(:)
      logical SC(nt)
      logical, allocatable::locmin_logic(:)
      logical, allocatable::incr(:)
      logical, allocatable::decr(:)
      
      parameter(eps=1.0d-17,MIT=200,TolX=1.0d-15) ! for BrentRoots

      call CapriniC(t,z,w,p,nz,nt,C) 

      do i=1,nt
         inds(i)=i
      end do

      allocate(DC(nt-1))
      allocate(incr(nt-1))
      allocate(decr(nt-1))
      allocate(allinds(nt))
      allocate(locmin_logic(nt-2))
      call diff(C,nt,DC)
      incr=DC>0 ! incr(j)=T iff C(j+1)>C(j)
      decr=DC<0 ! decr(j)=T iff C(j+1)<C(j)
      locmin_logic=incr(2:nt-1).and.decr(1:nt-2) 
c     locmin_ind(j)=T iff C(j+2)>C(j+1) and C(j+1)<C(j)
      allinds=inds(1:nt-2)
      num_locmins=count(locmin_logic)
      allocate(locmin_ind(num_locmins))
      locmin_ind=pack(allinds,locmin_logic) 
      lmc=0
      do k=1,num_locmins
         j=locmin_ind(k)
         t3(1)=t(j)
         t3(2)=t(j+1)
         t3(3)=t(j+2)
         nt3=3
         call vecCprime(t3,z,w,p,nz,nt3,Cpr)
         if ((Cpr(2)>0).NEQV.(Cpr(1)>0)) then
            x1=t(j)
            x2=t(j+1)
            locmin=BrentRoots(x1,x2,TolX,MIT,VAT,ni,er,z,w,p,nz)
            lmc=lmc+1       
            tk(lmc)=locmin
         else
            if ((Cpr(2)>0).NEQV.(Cpr(3)>0)) then
               x1=t(j+1)
               x2=t(j+2)
               locmin=BrentRoots(x1,x2,TolX,MIT,VAT,ni,er,z,w,p,nz)
               lmc=lmc+1
               tk(lmc)=locmin
            end if
         end if
      end do
      return
      end subroutine CapriniCT

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c subroutine Caprinifix(t,z,w,p,nz,nt,tk0,C,lmc,repeat)
c Inputs: t(1:nt)=the discretization of [0,+infty);
c z(1:nz) and w(1:nz) experimental data, p(1:nz) projection of w onto V(z)
c p(j)=f(z(j)), w(j)=f(z(j)) + noise, where f is a Stieltjes function
c Outputs: lmc=number of points in support of the spectral measure
c tk(1:lmc)=the support of the spectral measure, w(1:nz)=alterntive data
c C(1:nt)=Caprini function C(t)=Re[sum(p(j)-w(j))/(t-conj(z(j)))].
c
c In general, C(tk(k)) is not zero. The idea is to fix the noisy data w.
c We look for a point w~ near w s.t. the nodes tk are exactly optimal.
c We output C(s)=Re[sum(p(j)-w~(j))/(t(s)-conj(z(j)))].
c If p above is optimal, then
c C(t) must be a nonnegative function with zeros at tk
c We look for w~=w+dw.
c To preserve C'(tk)=0 we require Re[sum(dw(j)/(tk-conj(z(j)))^2)]=0
c To ensure that C(tk)=0 we require Re[sum(dw(j)/(tk-conj(z(j))))]=C(tk)
c If gamma>0 then Re[sum(dw(j))]=Re[sum(p(j)-w(j))]
c If the original C(0)<0, then we prepend 0 to the vector tk and require
c that -Re[sum(dw(j)/conj(z(j))=C(0)

      subroutine Caprinifix(t,z,w,p,nz,nt,tk0,C,lmc,repeat)
! Inputs: t,z,w,p,nz,nt
! Outputs: w,tk0,C,lmc,repeat
! Both inputs and outputs: w
      implicit none
      real*8 t(nt),tk0(nt),C(nt),m0,m1,C0,binf,tol,mag_C,RCOND
      real*8 re_dlt(nz),im_dlt(nz),abs_dlt(nz),QWORK(1)
      complex*16, allocatable::Mt(:,:)
      complex*16, allocatable::MtH(:,:)
      complex*16, allocatable::M2(:,:)
      real*8, allocatable::tk(:)
      real*8, allocatable::WORK(:)
      real*8, allocatable::Mr(:,:)
      real*8, allocatable::M2r(:,:)
      real*8, allocatable::Ct(:)
      real*8, allocatable::A(:,:)
      real*8, allocatable::b(:)
      real*8, allocatable::x(:)
      real*8, allocatable::BB(:,:)
      complex*16, allocatable::Ctc(:)
      complex*16 z(nz),w(nz),p(nz),dlt(nz),dw(nz)
      integer nt,nz,lmc,i,j,RANK,LW,INFO,nrhs,ldb,nr,ncc,lmc_plus
      integer zero,JPVT(2*nz),tmp(1)
      logical repeat,gamma,fix_zero

      parameter(tol=1.0d-11,RCOND=1.0d-14) 

      dlt=conjg(p-w)
      re_dlt=realpart(dlt)
      im_dlt=imagpart(dlt)
      abs_dlt=sqrt(re_dlt*re_dlt+im_dlt*im_dlt)
      mag_C=sum(abs_dlt/imagpart(z))/nz ! Average size of |C(t)|, t>0.

c Decide if we want to enforce optimality at infinity (gamma>0 or gamma=0)
      gamma=.false.
      m0=sum(re_dlt)
      m1=sum(realpart(dlt*z))
c If m0 is near 0 (indicating that gamma>0), but m1 is substantially negative, 
c then forcing m0 to be 0 will still violate optimality and so there is no
c point doing this anyway.
      if ((m0<tol*mag_C).and.(m1>-mag_C)) then
         gamma=.true.
      end if

c Compute the original Caprini function and its local minima
      call CapriniCT(t,z,w,p,nz,nt,C,tk0,lmc)
c There are lmc local minima stored in tk0

c Decide if t=0 should be in the support of the spectral measure
c if it is we add it to the list of local minima tk
      C0=-sum(realpart(dlt/z)) ! C(0)
      fix_zero=.false.
      zero=0
      if (C0<tol*mag_C) then
         fix_zero=.true.
         zero=1
      end if

c Compute equations for dw
      lmc_plus=lmc+zero
      allocate(Mt(nz,lmc_plus))
      allocate(MtH(lmc_plus,nz))
      allocate(M2(nz,lmc))
      allocate(Mr(2*nz,lmc_plus))
      allocate(M2r(2*nz,lmc))
      allocate(Ctc(lmc_plus))
      allocate(Ct(lmc_plus))
      allocate(tk(lmc_plus))
      if (fix_zero) then
         Mt(:,1)=-1.0/z
      end if
      do i=1,nz
         do j=1,lmc
            Mt(i,j+zero)=1.0/(tk0(j)-z(i)) 
            M2(i,j)=Mt(i,j+zero)**2
         end do
      end do      
      MtH=transpose(Mt)
      call zmatvect_prod(MtH,dlt,lmc_plus,nz,Ctc)
      Ct=realpart(Ctc) ! Ct=C(tk0)
      call mcomplex2real(Mt,Mr,nz,lmc_plus) ! Mr is 2nz x lmc or 2nz x (lmc+1)
      call mcomplex2real(M2,M2r,nz,lmc);   ! M2r= is 2nz x lmc      
      
      if (gamma) then
         nr=2*lmc+1+zero
         allocate(A(nr,2*nz))
         binf=sum(realpart(dlt))
         allocate(b(nr))
         b(nr)=binf
         do i=1,nz
            A(nr,2*i-1)=1.0
            A(nr,2*i)=0.0
         end do
      else
         nr=2*lmc+zero
         allocate(A(nr,2*nz))
         allocate(b(nr))
      end if
      A(1:lmc_plus,:)=transpose(Mr) ! check if an assignment can be done like this
      A(lmc_plus+1:2*lmc+zero,:)=transpose(M2r)         
      b(1:lmc_plus)=Ct
      do i=1,lmc
         b(lmc_plus+i)=0.0
      end do

c Solve equations for dw using min norm lsq: LAPACK's DGELSY
      ncc=2*nz
      nrhs=1
      LDB=max(nr,ncc)
      allocate(BB(LDB,nrhs))
      BB(1:nr,1)=b
      do i=1,2*nz
         JPVT(i)=0
      end do


      LW=-1  ! workspace query
      call dgelsy(nr,ncc,NRHS,A,nr,BB,LDB,JPVT,RCOND,RANK,QWORK,LW,INFO)
      LW=QWORK(1)

      allocate(WORK(LW))
      call dgelsy(nr,ncc,NRHS,A,nr,BB,LDB,JPVT,RCOND,RANK,WORK,LW,INFO)
      allocate(x(ncc))
      x=BB(1:ncc,1)
      call real2complex(x,dw,nz)
      w=w+dw
      
c Recompute the Caprini function for the alternative data

      call CapriniC(t,z,w,p,nz,nt,C)

      if (fix_zero) then
         lmc=lmc+1
         tk(1)=0.0
         tk(2:lmc)=tk0(1:lmc-1)
      else
         tk=tk0(1:lmc) 
      end if

c Figure out if Caprini fix should be repeated using new values of C
      repeat=.false.
      if (minval(C)<-tol*mag_C) then
         write(6,*) 'defect = ', minval(C)/mag_C
         repeat=.true.
      end if
      tk0(1:lmc)=tk ! return the support of spectral measure
      return
      end subroutine Caprinifix

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Returns the value t_max, such that the Caprini function
c C(t)=Re[sum((p(i)-w(i))/(t-conjg(z(i))) has no local 
c minima on [t_max,+infinity). The analysis is in the paper
c by Yury Grabovsky "Reconstructing Stieltjes functions from 
c their approximate values: a search for a needle in a hay stack." 

      real*8 function Caprinitmax(z,w,p,n)
      implicit none
      integer n
      complex*16 z(n), w(n), p(n), dlta(n)
      real*8 Nz,N0,M0,D0,absz(n),absdlta(n)

      absz=sqrt((realpart(z)**2+imagpart(z)**2))
      Nz=maxval(absz)
      dlta=p-w
      absdlta=sqrt((realpart(dlta)**2+imagpart(dlta)**2))
      N0=sum(absz*absdlta)
      D0=sum(realpart(dlta*conjg(z)))
      M0=2.0d0*N0/abs(D0);
      Caprinitmax=(M0+1.0d0)*Nz
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Returns the value of C'(t) for a given scalar t, where
c C'(t)=-Re[sum((p(i)-w(i))/(t-conjg(z(i))^2. It is used
c by BrentRoots for finding its zeros.
c
      real*8 function Cprime(t,z,w,p,n)
      implicit none
      complex*16 z(n),w(n),p(n)
      integer n,i
      real*8 t
      Cprime=realpart((w(1)-p(1))/(t-conjg(z(1)))**2)
      do i=2,n
         Cprime=Cprime+realpart((w(i)-p(i))/(t-conjg(z(i)))**2)
      end do
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes the "alternative data" stored in w on the output for which
c p is the closest point in V(z), and which is also very close to the
c original w (input). Once this is done it treats the local minima of
c of the Caprini function (possibly with t=0) as the support of the
c optimal spectral measure (stored in the first lmc elements of the
c array tk0 on the output) and computes weights (stored in first lmc+1
c elements of wts), so that 
c p(j)=wts(1)+sum(wts(i+1)/(tk0(i)-z(j)),i=1,lmc)
c The input t of length nt is our discretization of [0,+infinity) and 
c C=C(t), where C(t)=Re[sum((p(i)-w(i))/(t-conjg(z(i)))] is the Caprini
c function. Status=true, if we succeed in finding the appropriate 
c "alternative data" w, and false, otherwise.

      subroutine FullCaprinifix(t,z,w,p,nz,nt,tk0,lmc,wts,C,status)
! Inputs: t,z,w,p,nz,nt
! Outputs: w,tk0,lmc,wts,C,status
! Both inputs and outputs: w
      USE nonneg_leastsq
      implicit none
      complex*16 z(nz),w(nz),p(nz)
      real*8 t(nt),wts(nt),C(nt),tk0(nt)
      real*8, allocatable::sigmas(:)
      real*8, allocatable::tk(:)
      integer i,j,nz,nt,lmc,lmc_plus
      integer num_repeats,max_repeats,zero
      logical repeat,status

      status=.true.
      repeat=.true.
      num_repeats=0
      max_repeats=7
      do while (repeat.and.(num_repeats<=max_repeats))
         call Caprinifix(t,z,w,p,nz,nt,tk0,C,lmc,repeat)
         num_repeats=num_repeats+1
      end do
      if (repeat) then
c Cannot get the spectral measure after num2str(max_repeats) Caprini fixes
c Compute and display the defect             
         status=.false.
         write(6,*) 'Maximum of',max_repeats,'repeats is reached'
         write(6,*) 'Increasing npd may help.' 
         write(6,*) 'Even 1% increase may be effective.'
      end if
      
c Compute weights for the spectral measure. sigmas(1) is gamma. 
      zero=0
      if (tk0(1).eq.0.0) then
         zero=1
      end if
      call CapriniCT(t,z,w,p,nz,nt,C,tk0,lmc)
      lmc_plus=lmc+zero
      allocate(tk(lmc_plus))
      if (zero.eq.1) then
         tk(1)=0.0
      end if
      tk(zero+1:lmc_plus)=tk0(1:lmc)

      allocate(sigmas(lmc_plus+1))
      call getweights(z,p,tk,nz,lmc_plus,sigmas)
      wts(1:lmc_plus+1)=sigmas
      tk0(1:lmc_plus)=tk
      lmc=lmc_plus
      return
      end subroutine FullCaprinifix

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Given p in V(z) and the support of the spectral measure tk
c we compute the weights sigmas in the spectral representation
c p(j)=sigmas(1)+sum(sigmas(k+1)/(tk(k)-z(j)),k=1,lmc)

      subroutine getweights(z,p,tk,nz,lmc,sigmas)
! Inputs: z,p,tk,nz,lmc
! Outputs: sigmas
      USE nonneg_leastsq
      implicit none
      complex*16, allocatable::AA(:,:)
      real*8, allocatable::CC(:,:)
      complex*16 z(nz),p(nz)
      real*8 tk(lmc),sigmas(lmc+1),d(2*nz),rnorm,v(lmc+1)
      integer i,ncc,lmc,nz,nz2,mode,indx(lmc+1)

      ncc=lmc+1 ! sigmas(1) is a gamma
      allocate(AA(nz,ncc))
      allocate(CC(2*nz,ncc))
      do i=1,nz
         AA(i,1)=1.0
         AA(i,2:ncc)=1.0/(tk-z(i))
      end do

      call mcomplex2real(AA,CC,nz,ncc)
      call vcomplex2real(p,d,nz);
      nz2=2*nz
      call nnls(CC,nz2,ncc,d,sigmas,rnorm,v,indx,mode)
      return
      end subroutine getweights

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compute Nevanlinna-Pick matrices N(z,w) and P(z,w), given by
c N(i,j)=(w(i)-conjg(w(j)))/(z(i)-conjg(z(j))), and
c P(i,j)=(z(i)*w(i)-conjg(z(j)*w(j)))/(z(i)-conjg(z(j)))
c
      subroutine npmat(z,w,n,nmat,pmat)
! Inputs: z,w,n
! Outputs: nmat,pmat
      implicit none
      complex*16 z(n),w(n),nmat(n,n),pmat(n,n)
      integer n,i,j
      
      do i=1,n
         do j=1,n
            nmat(i,j)=(w(i)-conjg(w(j)))/(z(i)-conjg(z(j)))
            pmat(i,j)=(z(i)*w(i)-conjg(z(j)*w(j)))/(z(i)-conjg(z(j)))
         end do
      end do
      return
      end subroutine npmat

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes Nr randomly generated interpolants for the nz feasible data (z0,w0)
c evaluated  at nz points z0 and at additional ne points z. The data w0
c are perturbed by N(0,sgm^2) noise. The output WMC has Nr rows of length
c ne+nz, storing each random interpolant.
c Additional inputs are npd=number of points per decade for discretization
c of [0,+infinity), and ,mar govern this discretization. Here these parameters
c are carried passively to be passed to other subroutines down the stream.
c
      subroutine NPmonte_carlo(z,z0,w0,sgm,nz,ne,Nr,WMC,npd,mar)
! Inputs: z,z0,w0,sgm,nz,ne,Nr,npd,mar
! outputs: WMC
      implicit none
      complex*16 z(ne),z0(nz),w0(nz),zext(ne+nz),eye,w(ne+nz)
      complex*16 WMC(Nr,ne+nz),wn(nz)
      real*8 sgm,rand_re(Nr,nz),rand_im(Nr,nz),mar
      integer nz,ne,Nr,next,j,re_seed,im_seed,npd

      eye=cmplx(0,1)
      next=nz+ne
      zext(1:ne)=z
      zext(ne+1:next)=z0
! I don't know how to vary the seed dynamically
      re_seed=1234 
      im_seed=5678
      call r8mat_normal_01(Nr,nz,re_seed,rand_re)
      call r8mat_normal_01(Nr,nz,im_seed,rand_im)
      do j=1,Nr
         wn=w0+sgm*rand_re(j,:)+eye*sgm*rand_im(j,:)
         call Sinterp(z0,wn,zext,nz,next,w,npd,mar)
         WMC(j,:)=w
      end do
      end subroutine NPmonte_carlo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes the ad-hoc basis of V(z) stored in the rows of the matrix a
c The input is is a list tk of ts>0. The ad-hoc basis corresponds to 
c spectral measures delta(tk(i)) and chi_{[tk(i),m(i)]}(t), and 
c chi_{[m(i),tk(i+1)]}(t), where m(i)=0.5*(tk(i)+tk(i+1)).
c
      subroutine Sbasis(z,tk,a,nz,nt)
! Inputs: z,tk,nz,nt
! Outputs: a
      implicit none
      complex*16 z(nz),eye
      real*8 tk(nt),tmid
      complex*16 a(3*nt-1,nz) ! a is 3*nt-1 x nz complex matrix
      integer nt,nz,i,j
      
      eye=cmplx(0,1)
      do j=1,nz
         a(1,j)=1.0+0.0*eye         
         do i=2,nt
            tmid=0.5*(tk(i-1)+tk(i))
            a(i,j)=1.0/(tk(i-1)-z(j))
            a(nt+i,j)=log((tk(i)-z(j))/(tmid-z(j)))
            a(2*nt+i-1,j)=log((tmid-z(j))/(tk(i-1)-z(j)))
         end do
         a(nt+1,j)=1.0/(tk(nt)-z(j))
      end do
      return
      end subroutine Sbasis

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes H=f(z), where f is a Stieltjes
c function "satisfing" f(c_ex(i))=h_ex(i). The data may be noisy and we project
c it onto the set of valid data V(c_ex). We then proceed with interpolation,
c projecting onto V(c) at intermediate steps of recursion, if necessary.
c The recursion algorithm is described in a paper by Yury Grabovsky 
c "Reconstructing Stieltjes functions from their approximate values: a search 
c for a needle in a hay stack." Additional inputs npd, and mar govern 
c the discretization of [0,+infinity). Here these parameters are
c carried passively to be passed to myt_proj.

      recursive subroutine Sinterp(c_ex,h_ex,z,nc,ne,H,npd,mar)
! Inputs: c_ex,h_ex,z,nc,ne,npd,mar
! Outputs: H
      implicit none
      complex*16 c_ex(nc),h_ex(nc),z(ne),H(ne)
      complex*16 L11z(ne),L21z(ne),L22z(ne)
      complex*16, allocatable::z0(:)
      complex*16, allocatable::w0(:) 
      complex*16, allocatable::p(:)
      complex*16, allocatable::ck(:)
      complex*16, allocatable::hk(:)
      complex*16, allocatable::h_new(:)
      complex*16, allocatable::L11c(:)
      complex*16, allocatable::L21c(:)
      complex*16, allocatable::L22c(:)
      complex*16, allocatable::exit_check(:)
      complex*16, allocatable::c_rem(:)
      complex*16, allocatable::h_rem(:)
      complex*16 c1,h1
      real*8, allocatable::rep(:)
      real*8, allocatable::imp(:)
      real*8, allocatable::rext(:)
      real*8, allocatable::imxt(:)
      real*8, allocatable::absec(:)
      real*8 MC,rec1,imc1,reh1,imh1,imc1h1,mass,MH,L12,tol,f,mar
      integer nc,ne,nz,npd,i
      logical bad(nc),good(nc),Nval,Pval

      parameter(tol=2.0d-15) ! To determine when a denominator is too close to 0

c First we eliminate zero denominators in the NP matrices
      MC=maxval(imagpart(c_ex))
      bad=abs(imagpart(c_ex)/MC)<tol
      good=.not.bad
      nz=count(good)
      if (nz==0) then
         STOP 'Data is too badly conditioned'
      end if
      allocate(z0(nz))
      allocate(w0(nz))
      allocate(p(nz))
      allocate(rep(nz))
      allocate(imp(nz))
      allocate(ck(nz))
      allocate(hk(nz))

      z0=pack(c_ex,good)
      w0=pack(h_ex,good)
      
c exit clause
      if (nz==1) then
         c1=z0(1)
         h1=w0(1)
         call one_pt(c1,h1,z,ne,H)
c         write(6,*) 'Last exit at N=1'
         return
      end if

c Check validity and project onto V(z) if necessary
      call Valcheck(z0,w0,nz,Nval,Pval)
      p=w0
      if (.NOT.(Nval.AND.Pval)) then
c         write(6,*) 'Data lost validity. Projecting',nz
         call myt_proj(z0,w0,nz,p,npd,mar)
      end if

      rep=realpart(p)
      imp=imagpart(p)
      MH=sqrt(maxval(rep*rep+imp*imp))
      call reorder(z0,p,nz)
      c1=z0(1) 
      h1=p(1)

      allocate(c_rem(nz-1))
      allocate(h_rem(nz-1))
      c_rem=z0(2:nz) 
      h_rem=p(2:nz)
      nz=nz-1
      allocate(L11c(nz))
      allocate(L21c(nz))
      allocate(L22c(nz))
      allocate(exit_check(nz))
      allocate(rext(nz))
      allocate(imxt(nz))
      allocate(absec(nz))

c Apply the recursive procedure, but avoid divisions by zero
      imc1=imagpart(c1)
      rec1=realpart(c1)
      imh1=imagpart(h1)
      reh1=realpart(h1)
      imc1h1=imagpart(c1*h1)

      if (imc1h1/(MC*MH)<tol) then ! imag(c1*h1)=0 => H=-c1*h1/z
c         write(6,*) 'emergency exit 1',nz
         mass=realpart(c1*h1)
         H=mass/z
         return
      end if
      
      L22c=c_rem-(imc1*imc1+rec1*rec1)*imh1/imc1h1
      L21c=-imc1*c_rem/imc1h1

c Prepare remaining parameters
      if (imh1/MH<tol) then ! imag(h1)=0 => H=h1
c         write(6,*) 'emergency exit 2',nz
         H=reh1
         return
      end if
      L12=(reh1*reh1+imh1*imh1)*imc1/imh1
      L11c=c_rem-imc1h1/imh1

      exit_check=L22c+h_rem*L21c ! check if we need to exit early
      rext=realpart(exit_check)
      imxt=imagpart(exit_check)
      absec=sqrt(rext*rext+imxt*imxt)
      if (count(absec/MC<tol)>0) then
c         write(6,*) 'emergency exit 3',nz
         mass=(imc1*imc1+rec1*rec1)*imh1/imc1 ! sigma^{*}
         H=imc1h1/imc1-mass/z
         return
      end if

c Compute the values H(ck)
      allocate(h_new(nz))
      h_new=(L11c*h_rem+L12)/exit_check

c Invoke recursion
      call Sinterp(c_rem,h_new,z,nz,ne,H,npd,mar)
      L22z=z-(imc1*imc1+rec1*rec1)*imh1/imc1h1
      L21z=-imc1*z/imc1h1
      L11z=z-imc1h1/imh1
      H=(L22z*H-L12)/(L11z-H*L21z)
      end subroutine Sinterp

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Choose j such that w(j) that lies as close to the 
c bisector of the angle formed by a postivive real axis
c and the ray {w:Im(w*z(j))=0} that lies in the UHP.
c We regard such a w(j) as the one that "best satisfies"
c the constrains: Im(w(j))>0 and Im(w(j)*z(j))>0.

      subroutine reorder(z,w,nz)
! Inputs: z,w,nz
! Outputs: z,w
! Both inputs and outputs: z,w
      implicit none
      complex*16 z(nz),w(nz),wsqz(nz),savezw
      real*8 cosalpha(nz)
      integer nz,ind(1)
      wsqz=w*sqrt(-z)
      call absangle(wsqz,nz,cosalpha)
      ind=maxloc(cosalpha)
      savezw=z(1)
      z(1)=z(ind(1))
      z(ind(1))=savezw
      savezw=w(1)
      w(1)=w(ind(1))
      w(ind(1))=savezw
      return
      end subroutine reorder

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c The function computes Nevanlinna-Pick matrices N(cexp,hexp) and P(cexp,hexp)
c (see subroutine NPmat) and checks if they are nonnegative definite. 
c Nval=true, if N(cexp,hexp)>=0 and Pval=true if P(cexp,hexp)>=0.
c
      subroutine Valcheck(cexp,hexp,n,Nval,Pval)
! Inputs: cexp,hexp,n
! Outputs: Nval,Pval
      implicit none
      complex*16 cexp(n),hexp(n),Nw(n,n),Pw(n,n),qwork(1)
      complex*16, allocatable :: work(:) 
      integer n,ierr,lwork
      logical Nval,Pval
      real*8 frobnorm,Pnorm,Nnorm,nrmN,nrmP,eps,epsilon,rwork(3*n-2)
      real*8 N_tol,P_tol,eigsNw(n),eigsPw(n),mult
      character*1 jobz, UL

      parameter(mult=100.0d0,eps=1.d-16)

      call npmat(cexp,hexp,n,Nw,Pw)
      Nval=.true. 
      Pval=.true.
      Pnorm=frobnorm(Pw,n,n)
      Nnorm=frobnorm(Nw,n,n)
      nrmN=max(Nnorm,1/Nnorm)
      nrmP=max(Pnorm,1/Pnorm)
      N_tol=mult*nrmN*eps
      P_tol=mult*nrmP*eps

c     compute eigsNw and eigsPw

      jobz='N'
      UL='U'
      lwork=-1 ! workspace query
      call zheev(jobz,UL,n,Nw,n,eigsNw,qwork,lwork,rwork,ierr)
      lwork=qwork(1)
      allocate(work(lwork))
      call zheev(jobz,UL,n,Nw,n,eigsNw,work,lwork,rwork,ierr)
      call zheev(jobz,UL,n,Pw,n,eigsPw,work,lwork,rwork,ierr)
      if (eigsNw(1)<-N_tol) Nval=.false.
      if (eigsPw(1)<-P_tol) Pval=.false.
      return
      end subroutine Valcheck

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for a list of complex number z(k)=|z(k)|e^(i*phi(k)), |phi(k)|<=pi,
c the program computes the list cos(phi(k)). It is used to reorder data
c during the recursion algorithm (see subrounine Sinterp()).
c
      subroutine absangle(z,nz,cosphi)
! Inputs: z,nz
! Outputs: cosphi
      implicit none
      complex*16 z(nz)
      real*8 cosphi(nz),rez(nz),imz(nz),modz(nz)
      integer nz,i,numzeros
      logical zrs(nz)
      rez=realpart(z)
      imz=imagpart(z)
      modz=sqrt(rez*rez+imz*imz)
      do i=1,nz
         if (modz(i)==0) then
            cosphi(i)=-1.0 ! will not be chosen when we search for the maximum
         else
            cosphi(i)=rez(i)/modz(i)
         end if
      end do
      return
      end subroutine absangle
      
cccccccccccccccccccccccccccccccccccc
c
c Compute y(i)=x(i+1)-x(i), i=1,n-1
c
      subroutine diff(x,n,y)
! Inputs: x,n
! Outputs: y
      implicit none
      integer n,i
      real*8 y(n-1)
      real*8 x(n)
      do i=2,n
         y(i-1)=x(i)-x(i-1)
      end do
      return
      end subroutine diff

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Compute |x-y|, where x,y are two complex vectors
c
      real*8 function distCn(x,y,n)
      implicit none
      complex*16 x(n),y(n)
      real*8 dsq, re(n), im(n)
      integer n
      re=realpart(x-y)
      im=imagpart(x-y)
      dsq=sum(re*re+im*im)
      distCn=sqrt(dsq)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c The Frobenius norm of a complex m x n matric a: |a|^2=Tr(aa*)
c
      real*8 function frobnorm(a,m,n)
      implicit none
      complex*16 a(m,n)
      integer m,n,i,j
      frobnorm=0.0d0
      do i=1,m
         do j=1,n
            frobnorm=frobnorm+realpart(a(i,j))**2+imagpart(a(i,j))**2
         end do
      end do
      frobnorm=sqrt(frobnorm)
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c Equispaced discretization of an inerval [x1,x2] using n points.
c The list p of length n is the output
c
      subroutine linspace(x1,x2,n,p)
! Inputs: x1,x2,n
! Outputs: p
      implicit none
      integer n,i,m
      real*8 x1,x2,p(n),step
      step=(x2-x1)/(n-1)
      m=(n+1)/2 ! ceil(n/2)
      do i=1,m
         p(i)=x1+(i-1)*step
         p(n+1-i)=x2-(i-1)*step
      end do
      return
      end subroutine linspace

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c Equispaced discretization of an inerval [10^x1,10^x2] using n points
c on the logarithmic scale. The list p of length n is the output
c
      subroutine logspace(x1,x2,n,p)
! Inputs: x1,x2,n
! Outputs: p
      implicit none
      integer n,i,m
      real*8 x1,x2,p(n),step
      step=(x2-x1)/(n-1)
      m=(n+1)/2 ! ceil(n/2)
      do i=1,m
         p(i)=10.0**(x1+(i-1)*step)
         p(n+1-i)=10.0**(x2-(i-1)*step)
      end do
      return
      end subroutine logspace

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Multiply real nr x nc matrix C by a real vector x of length nr, returning
c a real vector p=Cx of length nc.
c
      subroutine matvect_prod(C,x,nr,nc,p)
! Inputs: C,x,nr,nc
! Outputs: p
      implicit none
      real*8 C(nr,nc),p(nr),x(nc)
      integer i,j,nc,nr
      do i=1,nr
         p(i)=0
         do j=1,nc
            p(i)=p(i)+C(i,j)*x(j)
         end do
      end do
      return
      end subroutine matvect_prod

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c converts columns of the complex nr x nc matrix uc into 
c real columns of 2*nr x nc matrix ur

      subroutine mcomplex2real(uc,ur,nr,nc)
! Inputs: uc,nr,nc
! Outputs: ur
      implicit none
      integer nr,nc,i
      real*8 ur(2*nr,nc)
      complex*16 uc(nr,nc)
      do i=1,nr
         ur(2*i-1,:)=realpart(uc(i,:))
         ur(2*i,:)=imagpart(uc(i,:))
      end do
      return
      end subroutine mcomplex2real

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes better "projection" p of w onto V(z)
c Additional inputs are npd=number of points per decade for discretization
c of [0,+infinity) and mar that govern this discretization. See comments in 
c t4Caprini
c
      subroutine myt_proj(z,w,nz,p,npd,mar)
! Inputs: z,w,nz,npd,mar
! Outputs: p
      implicit none
      complex*16 z(nz),w(nz),p(nz)
      integer nz,m,nt,ntk,nt4C,numdecades,npd,nat,nt4Clog,size_t4C
      real*8 t0(3*nz+2),mar,tmin,tmax,Caprinitmax
      real*8 t4Clin(npd),zero
      real*8, allocatable::t4C(:)
      real*8, allocatable::t4Clog(:)
      real*8, allocatable::ALLC(:)
      real*8, allocatable::sor_t1(:)
      real*8, allocatable::all_t1(:)
      real*8, allocatable::sor_t2(:)
      real*8, allocatable::all_t2(:)
      real*8, allocatable::t(:)
      real*8, allocatable::tk(:)

      call naivet(z,t0,nz,nt)
      allocate(t(nt))
      t=t0(1:nt)
      call simple_proj(z,w,t,nz,nt,p) ! compute mediocre projection p
      tmin=minval(imagpart(z)) ! Need to write Caprinitmin function
      tmax=Caprinitmax(z,w,p,nz) ! use mediocre p to estimate t interval for C(t)
      nt4C=size_t4C(tmin,tmax,mar,npd)

      allocate(t4C(nt4C))
      allocate(ALLC(nt4C))
      allocate(tk(nt4C))

      call t4Caprini(tmin,tmax,mar,npd,nt4C,t4C)
      call CapriniCT(t4C,z,w,p,nz,nt4C,ALLC,tk,ntk) ! compute all local minima of C(t)
      nat=ntk+nt
      allocate(all_t1(nat))
      all_t1(1:nt)=t
      all_t1(nt+1:nat)=tk(1:ntk)
      call unique(all_t1,nat,m)
      allocate(sor_t1(m))
      sor_t1=all_t1(1:m)
      call quicksort(sor_t1,m)
      call simple_proj(z,w,sor_t1,nz,m,p) ! compute better projection p
c Now repeat the process
      call CapriniCT(t4C,z,w,p,nz,nt4C,ALLC,tk,ntk) ! compute all local minima of C(t)
      nat=m+ntk
      allocate(all_t2(nat))
      all_t2(1:m)=sor_t1
      all_t2(m+1:nat)=tk(1:ntk)
      call unique(all_t2,nat,m)
      allocate(sor_t2(m))
      sor_t2=all_t2(1:m)
      call quicksort(sor_t2,m)
      call simple_proj(z,w,sor_t2,nz,m,p) ! compute even better p
      return
      end subroutine myt_proj

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Create a list of t's from the list of points z in the UHP. 
c The t's will be used to construct an ad-hoc basis for V(z).
c The number of t's is nts and they are stored in t(1:nts)
c nzs=# of elements in the array z. The rule is arcane, but it works.
c
      subroutine naivet(z,t,nzs,nts)
! Inputs: z,nzs
! outputs: t,nts
      implicit none
      integer nzr_plus,nzr_minus,nsrtd,nzs,nts
      complex*16 z(nzs)
      real*8 zr(nzs),absz(nzs),zrn(nzs),t(3*nzs+2),tmax,tol
      real*8, allocatable::zr_plus(:)
      real*8, allocatable::zr_minus(:)
      logical zrpl_ind(nzs),zrm_ind(nzs)

      parameter(tol=1.d-10)

c real parts, if they are positive
      zr=realpart(z)
      zrpl_ind=zr>0
      nzr_plus=count(zrpl_ind)
      zr_plus=pack(zr,zrpl_ind)

c     1/(negative real parts)
      absz=sqrt((realpart(z)**2+imagpart(z)**2))
      zrn=zr/absz
      zrm_ind=(zrn<-tol)
      zr_minus=-1.0/pack(zr,zrm_ind)
      nzr_minus=count(zrm_ind)

c put all values into a single array
      nsrtd=nzs+nzr_plus+nzr_minus+2
      t(1:nzr_plus)=zr_plus
      t(nzr_plus+1:nzr_plus+nzr_minus)=zr_minus

c imaginary parts
      t(nzr_plus+nzr_minus+1:nzr_plus+nzr_minus+nzs)=imagpart(z)
      t(nzr_plus+nzr_minus+nzs+1)=0.0
      tmax=maxval(t(1:nzr_plus+nzr_minus+nzs))
      t(nsrtd)=2.0*tmax

c sort t's and remove repeats
      call quicksort(t,nsrtd)
      call unique(t,nsrtd,nts) 
      return
      end subroutine naivet

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes a fractional-linear Stieltjes function f(zeta) satisfying f(z)=u
c for give complex numbers z and u in the complex UHP.
c The subroutine returns w=f(x), where x is the array of complex values in 
c the UHP.
c
      subroutine one_pt(z,u,x,nx,w)
! Inputs: z,u,w,x,nx
! Outputs: w
      implicit none
      complex*16 x(nx),eye,z,u,w(nx)
      real*8 t,s1,s2,g
      integer i,nx
c     nx=number of elements in x and w; z and u are complex scalars

      eye=cmplx(0,1)
      if (realpart(z*u)>=0 .and. realpart(u)<=0) then
         w=0.0
      else if (realpart(u)>=0 .and. imagpart(u)<=0) then
         w=realpart(u)
      else if (imagpart(z*u)<=0) then
         w=realpart(z*u)/x
      else if (imagpart(u*sqrt(-z))>0) then
         t=imagpart(z*u)/imagpart(u)
         s2=(realpart(u)**2+imagpart(u)**2)*imagpart(z)/imagpart(u)
         w=-s2/(x-t)
      else
         g=imagpart(z*u)/imagpart(z) ! gamma_{*}
         s1=(realpart(z)**2+imagpart(z)**2)*imagpart(u)/imagpart(z)
         w=g-s1/x
      end if
      return
      end subroutine one_pt

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Convert a real array ur of 2*nover2 elements into a 
c complex array uc of nover2 elements.
c
      subroutine real2complex(ur,uc,nover2)
! Inputs: ur,nover2
! Outputs: uc
      implicit none
      real*8 ur(2*nover2)
      complex*16 uc(nover2),eye
      integer nover2,i
      eye=cmplx(0,1)
      do i=1,nover2
         uc(i)=ur(2*i-1)+ur(2*i)*eye
      end do
      return
      end subroutine real2complex

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c computes simple "projection" pc of w onto V(z), using the ad-hoc basis
c phi_k(z) generated from z and t (see subrounine Sbasis). We look for x(k)>=0,
c such that pc=sum(x(k)*phi_k(z)) is as close to w as possible.
c
      subroutine simple_proj(z,w,t,nz,nt,pc)
! Inputs: z,w,t,nz,nt
! Outputs: pc
      USE nonneg_leastsq
      implicit none
      complex*16 z(nz),w(nz),pc(nz),A(3*nt-1,nz),B(nz,3*nt-1)
      real*8 t(nt),C(2*nz,3*nt-1),d(2*nz),x(3*nt-1),pr(2*nz)
      real*8 rnorm,v(3*nt-1),Csave(2*nz,3*nt-1)
      integer nz,nt,mode,indx(3*nt-1),ncc,nz2,i
      call Sbasis(z,t,A,nz,nt)
      ncc=3*nt-1
      nz2=2*nz
      B=transpose(A)
      call mcomplex2real(B,C,nz,ncc)
      call vcomplex2real(w,d,nz)
      Csave=C
      call nnls(C,nz2,ncc,d,x,rnorm,v,indx,mode)
      call matvect_prod(Csave,x,nz2,ncc,pr)
      call real2complex(pr,pc,nz)
      return
      end subroutine simple_proj

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Given data tmin,tmax,mar and npd we compute the number of
c points in the discretization of [0,+infty). We go left 
c mar decades from tmin, go right 1 decade from tmax and 
c use npd points per decade as well as in a linear discretization of
c [0,tmin/10^mar]. The function returns the number of points in such
c a discretization. This number must be computed before the 
c discretization, because Fortran requires to declare array sizes
c of inputs and outputs of subroutines
c
      integer function size_t4C(tmin,tmax,mar,npd)
      implicit none
      integer numdecades,npd
      real*8 tmin,tmax,mar,ltmin,ltmax

      ltmin=log10(tmin)
      ltmax=log10(tmax)
      numdecades=ceiling(ltmax-ltmin+mar+1.0)
      size_t4C=(numdecades+1)*npd
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Given data tmin,tmax,mar and npd we compute the discretization of 
c [0,+infty). We go left mar decades from tmin, go right 1 decade from 
c tmax and use npd points per decade as well as in a linear discretization
c  of [0,tmin/10^mar]. The discretization is returned in the array t4C. 
c
      subroutine t4Caprini(tmin,tmax,mar,npd,nt4C,t4C)
! Inputs: tmin,tmax,mar,npd,nt4C
! Outputs: t4C
      implicit none
      real*8 tmin,tmax,ltmin,ltmax,zero,mar,lintmin
      real*8 t4Clin(npd),t4Clog(nt4C+1-npd),t4C(nt4C)
      integer npd,nt4C,nt4Clog

      ltmin=log10(tmin)
      ltmax=log10(tmax)
      nt4Clog=nt4C+1-npd
      ltmin=ltmin-mar ! ltmin is an empirical bound
      ltmax=ltmax+1.0 ! ltmax is a theoretical bound
      zero=0.0
      lintmin=10.0**ltmin
      call linspace(zero,lintmin,npd,t4Clin)
      call logspace(ltmin,ltmax,nt4Clog,t4Clog) ! discretize t interval on a log scale
      t4C(1:npd)=t4Clin ! [0,tmin] equispaced grid
      t4C(npd:nt4C)=t4Clog ! [tmin,tmax] logarithmically equispaced grid
      return
      end subroutine t4Caprini
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c The input must be a sorted list a of length n. The program
c eliminates repeated values and the output is stored in the
c first m components of a.
c
      subroutine unique(a,n,m)
! Inputs: a,n
! Outputs: a,m
! Both inputs and outputs: a
      implicit none
      logical same
      logical index_a(n)
      integer i,j,n,m
      real*8 a(n)

      i=1
      index_a(1)=.true.
      do while( i<n )
         same=.true.
         j=i+1
         do while ( same .and. j<=n)
            if ( a(i)==a(j) ) then
               index_a(j)=.false.
               j=j+1
            else
               same=.false.
               index_a(j)=.true.
               i=j
            end if
         end do
         if ( j>n ) then
            i=n
         end if
      end do
      m=count(index_a)
      a(1:m)=pack(a,index_a)
      return
      end subroutine unique

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c converts the complex vector uc into the real vector ur of twice the length
c
      subroutine vcomplex2real(uc,ur,n)
! Inputs: uc,n
! Outputs: ur
      implicit none
      integer n,i
      real*8 ur(2*n)
      complex*16 uc(n)
      do i=1,n
         ur(2*i-1)=realpart(uc(i))
         ur(2*i)=imagpart(uc(i))
      end do
      return
      end subroutine vcomplex2real

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Returns the list of values values of C'(t) for a given list t,
c  where C'(t)=-Re[sum((p(i)-w(i))/(t-conjg(z(i))^2.
c
      subroutine vecCprime(t,z,w,p,n,nt,Cprime)
! Inputs: t,z,w,p,n,nt
! Outputs: Cprime
      implicit none
      complex*16 z(n),w(n),p(n)
      integer n,i,nt
      real*8 t(nt),Cprime(nt)
      Cprime=realpart((w(1)-p(1))/(t-conjg(z(1)))**2)
      do i=2,n
         Cprime=Cprime+realpart((w(i)-p(i))/(t-conjg(z(i)))**2)
      end do
      return
      end subroutine vecCprime

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Multiply complex nr x nc matrix C by a complex vector x of length nr, returning
c a complex vector p=Cx of length nc.
c
      subroutine zmatvect_prod(C,x,nr,nc,p)
! Inputs: C,x,nr,nc
! Outputs: p
      implicit none
      complex*16 C(nr,nc),p(nr),x(nc)
      integer i,j,nc,nr
      do i=1,nr
         p(i)=0
         do j=1,nc
            p(i)=p(i)+C(i,j)*x(j)
         end do
      end do
      return
      end subroutine zmatvect_prod

cccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c Computes the length |z| of the complex vector z.
c
      real*8 function znorm(x,n)
      implicit none
      complex*16 x(n)
      real*8 dsq, re(n), im(n)
      integer n
      re=realpart(x)
      im=imagpart(x)
      dsq=sum(re*re+im*im)
      znorm=sqrt(dsq)
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes poles tk and weights sigma of a function 
c f(z)=g-s0/z+sum(sigma(i)/(tk(i)-z),i=1,ns), provided it also satisfies
c f(c_ex(i))=h_ex(i). The data may be noisy and we project
c it onto the set of valid data V(c_ex). We then proceed with interpolation,
c projecting onto V(c) at intermediate steps of recursion, if necessary.
c The recursion algorithm is described in a paper by Yury Grabovsky 
c "Reconstructing Stieltjes functions from their approximate values: a search 
c for a needle in a hay stack." Additional inputs npd, and mar govern 
c the discretization of [0,+infinity). Here these parameters are
c carried passively to be passed to other subroutines down the stream. 

      recursive subroutine Spectr(c,h_ex,nc,npd,mar,g,s0,sigma,tk,ns)
! Inputs: c,h_ex,nc,npd,mar
! Outputs: g,s0,sigma,tk,ns
      implicit none
      complex*16 c(nc),h_ex(nc),c1,h1
      complex*16, allocatable::z0(:)
      complex*16, allocatable::w0(:) 
      complex*16, allocatable::p(:)
      complex*16, allocatable::ck(:)
      complex*16, allocatable::hk(:)
      complex*16, allocatable::h_new(:)
      complex*16, allocatable::L11c(:)
      complex*16, allocatable::L21c(:)
      complex*16, allocatable::L22c(:)
      complex*16, allocatable::exit_check(:)
      complex*16, allocatable::c_rem(:)
      complex*16, allocatable::h_rem(:)
      real*8, allocatable::rep(:)
      real*8, allocatable::imp(:)
      real*8, allocatable::rext(:)
      real*8, allocatable::imxt(:)
      real*8, allocatable::absec(:)
      real*8, allocatable::sgm(:)
      real*8, allocatable::ptk(:)
      real*8, allocatable::tf(:)
      real*8, allocatable::sf(:)
      real*8 MC,rec1,imc1,reh1,imh1,imc1h1,mass,MH,L12,Tol,stol
      real*8 f,mar,sigma(nc),tk(nc)
      real*8 tk1,t,g,s0,gst,tst,s1,s2,sg
      integer nc,ne,nz,npd,i,ns
      logical bad(nc),good(nc),Nval,Pval
      logical, allocatable::mask(:)
      
      parameter(stol=1.d-16,tol=2.0d-15)

c First we eliminate zero denominators in the NP matrices
      MC=maxval(imagpart(c))
      bad=abs(imagpart(c)/MC)<tol
      good=.not.bad
      nz=count(good)
      if (nz==0) then
         STOP 'Data is too badly conditioned'
      end if
      allocate(z0(nz))
      allocate(w0(nz))
      allocate(p(nz))
      allocate(rep(nz))
      allocate(imp(nz))
      allocate(ck(nz))
      allocate(hk(nz))

      z0=pack(c,good) ! z0,w0 are the effective data
      w0=pack(h_ex,good)
      
c exit clause
      if (nz==1) then
         c1=z0(1)
         h1=w0(1)
         call spec_one_pt(c1,h1,g,s0,sg,tk1,ns) ! ns=0 or 1
         if (ns==1) then ! if ns=0 sg and tk1 are not defined
            sigma(1)=sg
            tk(1)=tk1
         end if
c         write(6,*) 'Last exit at N=1'
         return
      end if

c Check validity and project onto V(z) if necessary
      call Valcheck(z0,w0,nz,Nval,Pval)
      p=w0
      if (.NOT.(Nval.AND.Pval)) then
         call myt_proj(z0,w0,nz,p,npd,mar) ! causes problems
      end if

      rep=realpart(p)
      imp=imagpart(p)
      MH=sqrt(maxval(rep*rep+imp*imp))
      call reorder(z0,p,nz) ! Now it is z0 and p that are regarded as data
      c1=z0(1) 
      h1=p(1)

      allocate(c_rem(nz-1))
      allocate(h_rem(nz-1))
      c_rem=z0(2:nz) 
      h_rem=p(2:nz)
      nz=nz-1
      allocate(L11c(nz))
      allocate(L21c(nz))
      allocate(L22c(nz))
      allocate(exit_check(nz))
      allocate(rext(nz))
      allocate(imxt(nz))
      allocate(absec(nz))
ccccccccccccccccccccccccccccccccccccccccccc

c Apply the recursive procedure, but avoid divisions by zero
      imc1=imagpart(c1)
      rec1=realpart(c1)
      imh1=imagpart(h1)
      reh1=realpart(h1)
      imc1h1=imagpart(c1*h1)      

      if (imc1h1/(MC*MH)<tol) then ! imag(c1*h1)=0 => H=-c1*h1/z
c         write(6,*) 'emergency exit 1',nz
         mass=realpart(c1*h1)
c        H=mass/z
         ns=0
         g=0.0
         s0=-mass
         return
      end if
      gst=imc1h1/imc1
      s2=(imc1*imc1+rec1*rec1)*imh1/imc1
      tst=imc1h1/imh1
! f=(L22*g-L12)/(L11-g*L21): L21=-z, L11=z-tst, L12=s1, L22=gst*z-s2


      L22c=gst*c_rem-s2
      L21c=-c_rem

c Prepare remaining parameters
      if (imh1/MH<tol) then ! imag(h1)=0 => H=h1
c         write(6,*) 'emergency exit 2',nz
c         H=reh1
         ns=0
         g=reh1
         s0=0.0
         return
      end if
      L12=(reh1*reh1+imh1*imh1)*imc1/imh1 ! =s1
      L11c=c_rem-tst

      exit_check=L22c+h_rem*L21c ! check if we need to exit early
      rext=realpart(exit_check)
      imxt=imagpart(exit_check)
      absec=sqrt(rext*rext+imxt*imxt)
      if (count(absec/MC<tol)>0) then
c         write(6,*) 'emergency exit 3',nz
c         H=gst-s2/z
         ns=0
         g=gst
         s0=s2
         return
      end if

c Compute the values h(ck)
      allocate(h_new(nz))
      h_new=(L11c*h_rem+L12)/exit_check

c Invoke recursion
      allocate(sgm(nz))
      allocate(ptk(nz))
      call Spectr(c_rem,h_new,nz,npd,mar,g,s0,sgm,ptk,ns)

! The first ns components of sgm and ptk contain positive weights and poles of g(z)
      allocate(tf(ns+1))
      allocate(sf(ns+1))
      do i=1,ns ! if ns=0 then we do not set tf and sf
         tf(i)=ptk(i)
         sf(i)=sgm(i)
      end do
      call getpoles(g,s0,tf,sf,ns,s1,s2,tst,gst) ! ns is now ns+1
      if (ns==1) then ! both weights and poles are determined
         if (sf(1)>stol) then
            tk(1)=tf(1)
            sigma(1)=sf(1)
         else ! the pole is discarded
            ns=0
         end if
         return
      end if

! now ns>1 and we need to determine the weights
      p=p-g+s0/z0
      nz=nz+1
      call getweights0(z0,p,tf,nz,ns,sf)

! Finally we need to discard negligible weights
      mask=sf>stol
      ns=count(mask)
      sigma(1:ns)=pack(sf,mask)
      tk(1:ns)=pack(tf,mask)
      return
      end subroutine Spectr

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Computes a fractional-linear Stieltjes function f(zeta) satisfying f(z)=u,
c where z and u are complex scalars in the UHP.
c The subroutine returns the spectral decomposition of f(zeta):
c f(zeta)=g-s0/zeta+sigma/(tk-zeta). If ns=0, the last term is not present
c and sigma and tk are not computed.

      subroutine spec_one_pt(z,u,g,s0,sigma,tk,ns) 
! Inputs: z,u
! Outputs: g,s0,sigma,tk,ns
      implicit none
      complex*16 eye,z,u
      real*8 t,s1,s2,g,sigma,tk,s0
      integer ns

      eye=cmplx(0,1)
      if (realpart(z*u)>=0 .and. realpart(u)<=0) then
c        w=0
         ns=0
         g=0.0
         s0=0.0
      else if (realpart(u)>=0 .and. imagpart(u)<=0) then
c        w=re(u)
         ns=0
         g=realpart(u)
         s0=0.0
      else if (imagpart(z*u)<=0) then
c        w=re(z*u)/x
         ns=0
         s0=-realpart(z*u)
         g=0.0
      else if (imagpart(u*sqrt(-z))>0) then
c        w=s2/(tk-x)
         ns=1
         g=0.0
         s0=0.0
         tk=imagpart(z*u)/imagpart(u)
         sigma=(realpart(u)**2+imagpart(u)**2)*imagpart(z)/imagpart(u)
      else
c        w=g-s0/x
         ns=0
         g=imagpart(z*u)/imagpart(z) ! gamma_{*}
         s0=(realpart(z)**2+imagpart(z)**2)*imagpart(u)/imagpart(z)
      end if
      return
      end subroutine spec_one_pt

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c See the paper for description of the algorithm.
c We compute the spectral representation of a Stieltjes 
c function f(z), satisfying f(z_j)=w_j, j=1,N. This is done
c recursvively using the formula 
c f(z)=[(gst*z-s2)*g(z)-s1]/[z*g(z)+z-tst]. Assuming that
c g(z) is given by g-s0/z+sum(s(k)/(t(k)-z),k=1,n) we compute
c the corresponding spectral representation of f(z) stored in 
c the same input variables g,s0,t,s. Variable n is n+1 on the output
c 
      subroutine getpoles(g,s0,t,s,n,s1,s2,tst,gst)
! Inputs: g,s0,t,s,n,s1,s2,tst,gst
! Outputs: g,s0,t,s,n
! Both inputs and outputs: g,s0,t,s,n
      implicit none
      real*8 g,s0,t(n+1),s(n+1),s1,s2,tst,gst,nu1,xs(n+2)
      real*8 term1,term2,Tmax,x1,x2,getapole,spectmax,tn
      real*8, allocatable::gpoles(:)
      real*8, allocatable::gwts(:)
      integer n,i
      if (n==0) then
         n=1
         t(1)=(s0+tst)/(g+1)
         term1=(s2*g+s1+gst*s0)/(g+1)
         term2=(g*gst*t(1)+s0*s2/t(1))/(g+1)
         s(1)=term1-term2
         s0=s0*s2/(s0+tst)
         g=gst*g/(g+1)
         return
      end if
! if n>0 then s is not returned, only t is computed
      xs(1)=0
      xs(2:n+1)=t(1:n)
!      Tmax=max(2*t(n),(tst+2*sum(s(1:n))+2*s0)/(g+1)) ! old Tmax
      tn=t(n)
      Tmax=spectmax(s,n,tn,s0,tst,g) ! New much smaller Tmax
      Tmax=Tmax+1.0 ! Make a little room if Tmax is too close to the actual zero of phi(x)
      xs(n+2)=Tmax
      allocate(gpoles(n))
      allocate(gwts(n))
      gpoles=t(1:n)
      gwts=s(1:n)
      do i=1,n+1
         x1=xs(i)
         x2=xs(i+1)
! Compute the zero of phi(z)=z*g(z)+z-tst in (x1,x2)
         t(i)=getapole(x1,x2,g,s0,gpoles,gwts,n,tst)
      end do      
      n=n+1
      s0=s0*s2/(s0+tst)
      g=gst*g/(g+1)
      return
      end subroutine getpoles

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c The poles of f(z)=[(gst*z-s2)*g(z)-s1]/[z*g(z)+z-tst] lie between the poles
c of g(z). For two specific consecutive poles x1 and x2 of g(z) the function 
c computes the pole of f(z) in (x1,x2). The smallest x1 is always 0 (regardless
c of whether g(z) has a pole at 0), and the largest x2 is given by Tmax in the
c subroutine getpoles. It is larger than the largest pole of g(z).

      real*8 function getapole(x1,x2,g,s0,gpoles,gwts,n,tst)
      implicit none
      real*8 x0,x1,x2,g,s0,gpoles(n),gwts(n),tst,phifun,a,b,f,f0,prod,rt
      real*8 BRphi,Minimum,tol,vat
      integer n,er,ni,maxiter,ii,RootBracketed,mit
! phi(x1)=-infinity and phi(x2)=+infinity

      a=0.5*(x1+x2)
      b=a
      f0=phifun(a,g,s0,gpoles,gwts,n,tst)
      prod=1.0
      maxiter=100 ! if maxiter is reached an error will be thrown by BRphi
      ii=0
      if (f0>0) then
         do while ((prod>=0).and.(ii<=maxiter))
            ii=ii+1
            a=0.5*(a+x1)
            f=phifun(a,g,s0,gpoles,gwts,n,tst)
            prod=f0*f
         end do
      else if (f0<0) then
         do while ((prod>=0).and.(ii<=maxiter))
            ii=ii+1
            b=0.5*(b+x2)
            f=phifun(b,g,s0,gpoles,gwts,n,tst)
            prod=f0*f
         end do
      else
         getapole=x0
         return
      end if
      tol=2.0e-17
      mit=200
      rt=BRphi(a,b,tol,mit,vat,ni,er,g,s0,gpoles,gwts,n,tst)
      if (.not.(er==0)) then
c         write(6,*) er,a,b,rt,vat
         STOP 'Failure to find a pole of f(z) between poles of g(z)'
      end if
      getapole=rt
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c The poles of f(z)=[(gst*z-s2)*g(z)-s1]/[z*g(z)+z-tst] are the
c zeros of phifun(z)=z*g(z)+z-tst. All other input arguments are
c parameters in the spectral representation of g(z):
c g(z)=g-s0/z+sum(gwts(k)/(gpoles(k)-z),k=1,n)

      real*8 function phifun(z,g,s0,gpoles,gwts,n,tst)
      implicit none
      real*8 z,g,s0,gpoles(n),gwts(n),tst,s
      integer n
! phi=z*g(z)+z-tst, g(z)=g-s0/z+sum(gwts/gpoles-z)
      s=sum(gwts/(gpoles-z))
      phifun=g*z+z-tst-s0+z*s
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Given p in V(z) and the support of the spectral measure tk
c we compute the weights sigmas in the spectral representation
c p(j)=sum(sigmas(k)/(tk(k)-z(j)),k=1,lmc), where all tk(k)>0
c This subroutine is invoked only when p is known to have the 
c indicated representation.

      subroutine getweights0(z,p,tk,nz,lmc,sigmas)
! Inputs: z,p,tk,nz,lmc
! Outputs: sigmas
      USE nonneg_leastsq
      implicit none
      complex*16 z(nz),p(nz),AA(nz,lmc)
      real*8 tk(lmc),sigmas(lmc),d(2*nz),rnorm,v(lmc),CC(2*nz,lmc)
      integer i,lmc,nz,nz2,mode,indx(lmc)

      do i=1,nz
         AA(i,:)=1.0/(tk-z(i))
      end do
      call mcomplex2real(AA,CC,nz,lmc)
      call vcomplex2real(p,d,nz);
      nz2=2*nz
      call nnls(CC,nz2,lmc,d,sigmas,rnorm,v,indx,mode)
      return
      end subroutine getweights0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c We compute a better estimate for Tmax for the spectral representation
c of f(z) computed in getpoles routine. Better Tmax solves
c (g+1)*T^2-(sum(s)+s0+tst+tn*(g+1))*T+tn*(s0+tst)

      real*8 function spectmax(s,n,tn,s0,tst,g)
      implicit none
      integer n
      real*8 s(n),tn,s0,tst,g,b,c
      b=sum(s)+s0+tst+tn*(g+1)
      c=tn*(s0+tst)
      spectmax=(b+sqrt(b*b-4*(g+1)*c))/(2*(g+1))
      return
      end
