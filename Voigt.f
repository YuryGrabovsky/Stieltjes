!__________________________________________________________________________________
!
! Author: Yury Grabovsky https://www.math.temple.edu/~yury
! Last modified: Feb 14, 2021
! This code is distributed under the GNU LGPL license.
!
! It is an example of using the driver routine StCaprini in the 
! context of electrochemichal impedance spectroscopy or for 
! Voigt circuits - circuits made with resistors and capacitors
!
!__________________________________________________________________________________

!***********************************************************************
! Given nz distinct frequencies f_exp and approximate values w of
! Z(f_exp) for a a complex impedance function of an electrochemical cell
! we compute physically admissible function Z(f), such that |w-Z(f_exp)|
! is as small as possible. In addition to the answer we display the
! certificate of optimality in the form of the graph of the Caprini
! function. We display the graph of the Caprini function on the macro
! and the micro-scales.  The local minima tk, of the Caprini function
! must also be its zeros and the global minima, in which case they
! comprise the support of the spectral measure of the corresponding
! integral representation of Z(f): Z(f)=gamma + int_0^infinity d
! sigma(t)/(t+2i*pi*f). The optimal spectral measure must necessarily be
! finite, so that Z(f) = g + s0/(2i*pi*f) + sum(sk(k+1)/(tk(k)+2i*pi*f)). 
! For a given list of intermediate
! frequencies f we also compute the values Z(f) stored in the vector
! W_extr. There is an uncertainty associated with this reconstruction
! problem. It is identified by running NMC=500 Monte-Carlo simulations
! that have the same noise level as w (estimated from the error
! |w-Z(f_exp)|. The results are stored in the NMC rows of the output
! matrix WMC. We also compute and output the Caprini function to display
! the certificate of optimality.
!
!*******************************************************************

      program Voigt
      implicit none
      complex*16 eye
      real*8 pi,mar,tmin,tmax,rew,imw,VD(300,3)
      real*8 rezeta,imzeta,f_exp,fmin,fmax,sgm,g,s0
      integer nz,n_zeta,i,j,NMC,npd,mnd,ntC,numlm,ns,eof
      logical status
      complex*16, allocatable::z(:)
      complex*16, allocatable::w(:)
      complex*16, allocatable::p(:)
      complex*16, allocatable::zeta(:)
      complex*16, allocatable::W_extr(:)
      complex*16, allocatable::WMC(:,:)
      real*8, allocatable::f(:)
      real*8, allocatable::tC(:)
      real*8, allocatable::rtk(:)
      real*8, allocatable::sigma(:)
      real*8, allocatable::C(:)
      real*8, allocatable::Cmins(:)
      real*8, allocatable::locmins(:)

      parameter(mar=4.0,npd=1000,mnd=30)
! mar=number of decades to go left past tmim
! npd=number of discretization points per decade
! nz=number of experimental measurements
! n_zeta=number of extrapolation points

      eye=cmplx(0,1)
      pi=4*atan(1.0d0)
      
! Input "experimental" data f_exp,w and the list of extrapolation frequencies f
! size(f_exp)=size(w)=nz, size(f)=n_zeta
      n_zeta=300
      NMC=500
      open(unit=11,file='Voigt_data.txt') !
      fmin=1.0d24
      fmax=0.0
      do i=1,302 ! We will not handle more that 300 data points
         read(11,*,iostat=eof) f_exp,rew,imw
         if (eof<0) then ! end of file is reached
            nz=i-1
            exit
         end if
         if (f_exp<fmin) then
            fmin=f_exp
         end if
         if (f_exp>fmax) then
            fmax=f_exp
         end if
         VD(i,1)=f_exp 
         VD(i,2)=rew
         VD(i,3)=imw
      end do
      close(11)
      if (i>301) then
         STOP 'Too many measurements'
      end if
c      write(6,*) 'The number of measurements is',nz
      allocate(z(nz))
      allocate(w(nz))
      allocate(zeta(n_zeta))
      allocate(f(n_zeta))
      z=2*pi*eye*VD(1:nz,1)
      w=VD(1:nz,2)-eye*VD(1:nz,3)
      fmin=log10(fmin) ! -1.0 ! for Voigt_minus
      fmax=log10(fmax) ! +4.0 ! for Voigt_plus
      call logspace(fmin,fmax,n_zeta,f)
      zeta=2*pi*eye*f

! Data has been read in. Allocate output arrays

      allocate(p(nz))
      allocate(sigma(nz))
      allocate(rtk(nz))
      allocate(W_extr(n_zeta))
      allocate(WMC(NMC,n_zeta))
      allocate(tC(mnd*npd))
      allocate(C(mnd*npd))
      allocate(locmins(nz))
      allocate(Cmins(nz))

! Compute results

      call StCaprini(z,w,zeta,nz,n_zeta,NMC,npd,mar,mnd,
     $               p,W_extr,WMC,status,
     $               g,s0,rtk,sigma,ns,tC,C,ntC)

! Computations are done. Save results to text files

      if (status) then
         open(unit=15,file='tC.txt')
         do i=1,ntC
            write(15,*) tC(i), C(i)
         end do
         close(15)

         open(unit=16,file='wfix.txt')
         do i=1,nz
            write(16,*) realpart(w(i)),-imagpart(w(i))
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
         write(13,*) f(i),realpart(W_extr(i)),-imagpart(W_extr(i))
         do j=1,NMC
            write(14,*) realpart(WMC(j,i)),-imagpart(WMC(j,i))
         end do
      end do
      close(13)
      close(14)
      open(unit=18,file='data_sizes.txt')
      write(18,*) nz
      write(18,*) n_zeta
      write(18,*) NMC ! this could change to 0 if data is exact
      close(18)

      stop
      end 
! gfortran Voigt.f Stieltjes_functions.f *.f90 <path to liblapack.a>liblapack.a. <path to librefblas.a>librefblas.a
