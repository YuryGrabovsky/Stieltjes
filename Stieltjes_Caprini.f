!__________________________________________________________________________________
!
! Author: Yury Grabovsky https://www.math.temple.edu/~yury
! Last modified: Feb 21, 2021
! This code is distributed under the GNU LGPL license.
! It is an example of using the driver routine StCaprini
!
!__________________________________________________________________________________

!***********************************************************************
! Given nz distinct points z in the UHP and approximate values w of
! a Stieltjes function f(z) we compute p in V(z) which is certifiably
! closest to wfix(=w on the output). wfix is a slight perturbation of the actual data w so 
! p is optimal for wfix. We display the graph of Caprini function on macro 
! and micro-scales to certify optimality.
! The local minima tk, of the Caprini function is the support of the 
! spectral measure of the corresponding Stieltjes function f(z) satisfying
! f(z(j))=p(j). The numbers g,s0,sigma describe the spectral representation of
! f(zeta) = g-s0/zeta + sum(sigma(i)/(rtk(i)-zeta)). 
! 
! Vector zeta of points in the UHP is where we want to compute f(z) at. The
! vector W_extr is the vector of desired values. There is an uncertainty
! associated to this reconstruction problem. It is identified by running NMC
! Monte-Carlo simulations that have the same noise level as w-p. The results 
! are stored in the NMC rows of the output matrix WMC. 
!
!*******************************************************************

      program Stieltjes_Caprini
      implicit none
      complex*16 eye
      real*8 pi,mar,tmin,tmax,rew,imw
      real*8 rezeta,imzeta,rez,imz,sgm,g,s0
      integer nz,n_zeta,i,j,NMC,npd,mnd,ntC,numlm,ns,eof,data_max
      logical status
      complex*16, allocatable::z(:)
      complex*16, allocatable::w(:)
      complex*16, allocatable::p(:)
      complex*16, allocatable::zeta(:)
      complex*16, allocatable::W_extr(:)
      complex*16, allocatable::WMC(:,:)
      real*8, allocatable::VD(:,:)
      real*8, allocatable::t4C(:)
      real*8, allocatable::tC(:)
      real*8, allocatable::rtk(:)
      real*8, allocatable::C(:)
      real*8, allocatable::sigma(:)

      parameter(mar=3.0,npd=1000,mnd=30,data_max=100000)
! mar=number of decades to go left past tmim
! npd=number of discretization points per decade
! nz=number of experimental measurements
! n_zeta=number of extrapolation points
! mnd: maximal number of decades to reserve enough computer memory
      NMC=500
      eye=cmplx(0,1)
      pi=4*atan(1.0d0)
      
! Input "experimental" data z,w and the list of extrapolation points zeta
! size(z)=size(w)=nz, size(zeta)=n_zeta
      allocate(VD(data_max,4))
      open(unit=11,file='exp_data.txt') ! 
      do i=1,data_max
         read(11,*,iostat=eof) rez,imz,rew,imw
         if (eof<0) then ! end of file is reached
            nz=i-1
            exit
         end if
         VD(i,1)=rez
         VD(i,2)=imz
         VD(i,3)=rew
         VD(i,4)=imw
      end do
      close(11)
      allocate(z(nz))
      allocate(w(nz))      
      z=VD(1:nz,1)+eye*VD(1:nz,2)
      w=VD(1:nz,3)+eye*VD(1:nz,4)

      open(unit=21,file='extr_zs.txt') !
      do i=1,data_max
         read(21,*,iostat=eof) rezeta,imzeta
         if (eof<0) then ! end of file is reached
            n_zeta=i-1
            exit
         end if
         VD(i,1)=rezeta
         VD(i,2)=imzeta
      end do
      close(21)
      allocate(zeta(n_zeta))
      zeta=VD(1:n_zeta,1)+eye*VD(1:n_zeta,2)
      deallocate(VD)

! Data has been read in. Allocate output arrays

      allocate(p(nz)) ! projection of data onto the interpolation body
      allocate(sigma(nz)) ! weights if the spectral measure
      allocate(rtk(nz)) ! nodes of the spectral measure
      allocate(W_extr(n_zeta)) ! extrapolation via the recursive algorithm
      allocate(WMC(NMC,n_zeta)) ! Monte Carlo uncertainty quantification
      allocate(tC(mnd*npd)) ! discretizartion of the t-space for C(t)
      allocate(C(mnd*npd)) ! C(t)

! Compute results

      call StCaprini(z,w,zeta,nz,n_zeta,NMC,npd,mar,mnd,
     $               p,W_extr,WMC,status,
     $               g,s0,rtk,sigma,ns,tC,C,ntC)

! W_extr(j)=f(zeta(j)), given by an explicit representation
! WMC(i,j)=f_{i}(zeta(j)) for ith random realization f_i(z) of f
! t4C(s) is the discretization of (0,+infty) on log scale
! C(s)=C(t4C(s)) is the Caprini function
! w is wfix, p is the projection of wfix onto V(z), which is optimal as evidenced by the 
! graphs of C(t) on macro and microscales

! Computations are done. Save results to text files
      if (status) then
         open(unit=15,file='tC.txt')
         do i=1,ntC
            write(15,*) tC(i), C(i)
         end do
         close(15)

         open(unit=16,file='wfix.txt')
         do i=1,nz
            write(16,*) realpart(w(i)),imagpart(w(i))
         end do
         close(16)
      end if

      open(unit=12,file='spectral_measure.txt')
      if (status) then ! inform whoever reads the output of status
         write(12,*) -1.0, g
      else
         write(12,*) 0.0, g
      end if
      write(12,*) 0.0, s0
      do i=1,ns
         write(12,*) rtk(i),sigma(i)
      end do
      close(12)
      open(unit=13,file='W_extr.txt')
      open(unit=14,file='WMC.txt')
      do i=1,n_zeta
         write(13,*) realpart(zeta(i)), imagpart(zeta(i)),
     $   realpart(W_extr(i)),imagpart(W_extr(i))
         do j=1,NMC
            write(14,*) realpart(WMC(j,i)),imagpart(WMC(j,i))
         end do
      end do
      close(13)
      close(14)
      open(unit=17,file='pout.txt')
      do i=1,nz
         write(17,*) realpart(p(i)),imagpart(p(i))
      end do
      close(17)
      open(unit=18,file='data_sizes.txt')
      write(18,*) nz
      write(18,*) n_zeta
      write(18,*) NMC ! this could change to 0 if data is exact
      close(18)

      stop
      end 
! gfortran Stieltjes_Caprini.f  Stieltjes_functions.f *.f90 <path to liblapack.a>liblapack.a. <path to librefblas.a>librefblas.a
