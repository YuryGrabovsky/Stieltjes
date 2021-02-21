!* ------------------------------------------------- *
!* Reference:  BORLAND MATHEMATICAL LIBRARY          *
!*                                                   *
!*                F90 version by J-P Moreau, Paris.  *
!*                       (www.jpmoreau.fr)           *
!* ------------------------------------------------- *
!*              Brent Method Function                *
!* ------------------------------------------------- *
!* The purpose is to find a real root of a real      *
!* function f(x) using Brent method.                 *
!*                                                   *
!* INPUTS:  x1,x2     : interval of root             *
!*          Tolerance : desired accuracy for root    *
!*          maxIter   : maximum number of iterations *
!*                                                   *
!* OUTPUTS: The function returns the root value      *
!*          ValueAtRoot : value of f(root)           *
!*          niter    : number of done iterations     *
!*          error    : =0, all OK                    *
!*                   : =1, no root found in interval *
!*                   : =2, no more iterations !      *
!*                                                   *
!* The calling routine must have for some reason     *
!* real*8 BrentRoots,Minimum                         *
!* integer RootBracketed                             *
!*****************************************************
real*8 Function BrentRoots( x1, x2, Tolerance,  &
                            maxIterations,      &
		            valueAtRoot,        &
                            niter, error,       &
                            zz,ww,p,nz)         ! extra inputs into f(x)
                                                ! must be included in all 4 function calls

  parameter(FPP = 1.d-11, nearzero = 1.d-20)

  real*8  x1,x2,Tolerance,valueAtRoot,Cprime,Minimum
  integer niter, error,RootBracketed,nz

  real*8 resultat, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm
  integer i, done

  complex*16 zz(nz),ww(nz),p(nz)

  i = 0; done = 0;   error = 0
  AA = x1;  BB = x2;  FA = Cprime(AA,zz,ww,p,nz); FB = Cprime(BB,zz,ww,p,nz)
  if (RootBracketed(FA,FB).eq.0) then 
    error = 1
  else 
    FC = FB;
    do while (done.eq.0.and.i < maxIterations)
      if (RootBracketed(FC,FB).eq.0) then
        CC = AA; FC = FA; DD = BB - AA; EE = DD
      endif
      if (dabs(FC) < dabs(FB)) then
        AA = BB; BB = CC; CC = AA
        FA = FB; FB = FC; FC = FA
      endif
      Tol1 = 2.0 * FPP * dabs(BB) + 0.5 * Tolerance
      xm = 0.5 * (CC-BB)
      if ((dabs(xm) <= Tol1).or.(dabs(FA) < nearzero)) then
        ! A root has been found
        resultat = BB;
        done = 1
        valueAtRoot = Cprime(resultat,zz,ww,p,nz)
      else 
        if ((dabs(EE) >= Tol1).and.(dabs(FA) > dabs(FB))) then
          SS = FB/ FA;
          if (dabs(AA - CC) < nearzero) then
            PP = 2.0 * xm * SS;
            QQ = 1.0 - SS;
          else 
            QQ = FA/FC;
            RR = FB /FC;
            PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
            QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
          endif
          if (PP > nearzero) QQ = -QQ;
          PP = dabs(PP);
          if ((2.0 * PP) < Minimum(3.0*xm *QQ-dabs(Tol1 * QQ), dabs(EE * QQ))) then
            EE = DD;  DD = PP/QQ;
          else 
            DD = xm;   EE = DD;
          endif
        else 
          DD = xm;
          EE = DD;
        endif
        AA = BB;
        FA = FB;
        if (dabs(DD) > Tol1) then 
          BB = BB + DD;
        else 
          if (xm > 0) then 
	    BB = BB + dabs(Tol1)
          else 
	    BB = BB - dabs(Tol1)
          endif
        endif
        FB = Cprime(BB,zz,ww,p,nz)
        i=i+1
      endif
	end do
    if (i >= maxIterations) error = 2
  endif
  niter = i
  BrentRoots = resultat
end ! BrentRoots()

! TRUE if x1*x2 negative
integer Function RootBracketed(x1,x2)
  real*8 x1,x2 
  integer resultat
  if ((x1 > 0.and.x2 > 0).or.(x1 < 0.and.x2 < 0)) then 
    resultat = 0
  else
    resultat = 1
  endif
  RootBracketed = resultat
end

! returns the minimum of two real numbers
real*8 Function Minimum(x1,x2) 
  real*8 x1,x2,resultat
  if (x1 < x2) then
    resultat = x1
  else 
    resultat = x2
  endif
  Minimum = resultat
end

! End of file zbrent.f90
