subroutine pseudopotd(u,nn)
  USE basic_data_types, ONLY: dp
  USE atom, ONLY: namax,lamax,alpha,nalpha,lmax
  USE pspot
  IMPLICIT NONE
  REAL(dp) :: nn(namax,namax,0:lamax)
  REAL(dp) :: u(namax,namax,0:lamax)
  REAL(dp) :: b,bb,a,u1,u2,u3
  REAL(dp) :: GP2(KBPROJmax,namax,0:lamax)
  REAL(dp) :: gint
  INTEGER :: i,j,k,l,kappa,lambda

  do l=0,lmax
    do i=1,nalpha(l)
      do j=i,nalpha(l)
        a=alpha(i,l)+alpha(j,l)
      !.core potential
        u1=0.D0
        do k=1,ERFnr
          b=PPNerf(k)
          bb=PPNerf(k)**2
          if (l.eq.0) then
            u1=u1+PPLerf(k)*  b*(2*bb+3*a) / (4*a**2*(a+bb)**1.5D0)
          elseif (l.eq.1) then
            u1=u1+PPLerf(k)*  b*(8*bb**2+20*a*bb+15*a**2) /&
     						(8*a**3*(a+bb)**2.5D0)
          elseif (l.eq.2) then
            u1=u1+PPLerf(k)* 3*b*(16*bb**3+56*bb**2*a+70*bb*a**2+35*a**3) /&
     						(16*a**4*(a+bb)**3.5D0)
          elseif (l.eq.3) then
            u1=u1+PPLerf(k)*  3*b*(128*bb**4+576*bb**3*a&
                                   +1008*bb**2*a**2+840*bb*a**3+315*a**4) /&
     						(32*a**5*(a+bb)**4.5D0)
          endif
        enddo
      !.nonlocal part
        u2=0.D0
        do k=1,EXPnr(l)
          u2=u2+PPLexp(k,l)*gint(2*l+4+PPRexp(k,l),a+PPNexp(k,l))
        enddo
        u(i,j,l)=(-u1*Zeff+u2)*nn(i,j,l)
        u(j,i,l)=u(i,j,l)
      enddo
    enddo
  enddo

  if (pptype.eq.5) then

  !.calculate the 'overlap' integrals between basis and r**2*KB-projectors
    do l=0,lmax
      do kappa=1,nalpha(l)
        do i=1,KBPROJnr(l)
          u3=0.d0
          do k=1,KBEXPnr(l)
            u3=u3+KBLexp(k,i,l)&
                 *gint(KBRexp(k,i,l)+l+4,KBNexp(k,i,l)+alpha(kappa,l))
          enddo
          GP2(i,kappa,l)=u3*sqrt(nn(kappa,kappa,l))
        enddo
      enddo
    enddo

  !.put all together for V_kb
    do l=0,lmax
      do kappa=1,nalpha(l)
        do lambda=kappa,nalpha(l)
          do i=1,KBPROJnr(l)
            do j=1,KBPROJnr(l)
              u(kappa,lambda,l)&
                =u(kappa,lambda,l)&
                +0.5d0*GP(i,kappa,l)*GP2(j,lambda,l)*KBV(i,j,l)&
                +0.5d0*GP(j,kappa,l)*GP2(i,lambda,l)*KBV(i,j,l)
              if (lambda.ne.kappa) then
                u(lambda,kappa,l)&
                  =u(lambda,kappa,l)&
                  +0.5d0*GP2(i,kappa,l)*GP(j,lambda,l)*KBV(i,j,l)&
                  +0.5d0*GP2(j,kappa,l)*GP(i,lambda,l)*KBV(i,j,l)
              endif
            enddo
          enddo
        enddo
      enddo
    enddo

  endif

end subroutine pseudopotd
