        subroutine dnear(x,x0,x1,deltat,mu0,mu1,nt1,n,dn)
C   the near part of Double Layer Potential (1D) int_{t-2*deltat}^{t-deltat}
C   input: x, deltat: self-explaining
C       X0, x1: boundary points
C       mu0, mu1: mu(x0,t_i) and mu(x1, t_i)
C       nt1: index in mu0 and mu1 where the values of t-deltat lie...
C       n: N-th order Gauss-legendre rule
C   output: dn: the near part of the double layer potential
C-----------------------------------------------------------
        implicit real*8 (a-h,o-z)
        real*8 x,x0,x1,deltat,dn
        real*8 mu0(1),mu1(1)
        integer nt1,n

        real*8 s(10),w(10),mua(10),mub(10),va(10),vb(10)
C---------------------------------
        done=1.0d0
        pi=4*atan(done)
        open(unit=12,file='legquad.dat')
        do i=1,10
          read(12,*) s(i),w(i)
        enddo
        dn=0.0d0
        do i=1,10
          s(i)=(s(i)+1.0d0)/2*deltat+deltat
          w(i)=deltat*w(i)/2;
C           Gauss-legendre nodes and weights
          mub(i)=mu1(nt1-1)+(mu1(nt1)-mu1(nt1-1))/deltat*(2*deltat-s(i))
          mua(i)=mu0(nt1-1)+(mu0(nt1)-mu0(nt1-1))/deltat*(2*deltat-s(i))
C23456
          vb(i)=(x-x1)/(4*sqrt(pi)*s(i)**(3.0d0/2))*
     1                        exp(-(x-x1)*(x-x1)/(4*s(i)))*mub(i)
          va(i)=(x-x0)/(4*sqrt(pi)*s(i)**(3.0d0/2))*
     1                        exp(-(x-x0)*(x-x0)/(4*s(i)))*mua(i)
          dn=dn+w(i)*vb(i)-w(i)*va(i)
        enddo
        close(12)
C---------------------------------
        end subroutine



C-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
C-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.



        subroutine eval_dl(x,nt,deltat,a,b,mua,mub,dl)
C----------------------------
C   a function for the evaluation of DL, the local part of double layer potential
C   Input:
C   x: location where dl is to be evaluated, scalar,double
C   nt: integer, the index of current time in mua and mub
C   deltat: length of one time step
C   mua: density function mu(a,t), as a vector {mu(a,t_j)}_{j<=nt}
C   mub: density function mu(b,t), as a vector {mu(b,t_j)}_{j<=nt}
C   Output:
C   dl: the  evaluated local double layer potential
C------------------------------
        implicit real*8 (a-h,o-z)
        real*8 x,deltat,a,b,mua(1),mub(1),dl
        integer nt

        done=1.0d0        
        pi=4*atan(done)

        xa=x-a
        xb=x-b
        
C        A_1a=0.5d0*(1.0d0-erf(abs(xa)/sqrt(4*deltat)))*sign(1.0d0,xa)
C        A_1b=0.5d0*(1.0d0-erf(abs(xb)/sqrt(4*deltat)))*sign(1.0d0,xb)
C------------------------
        if(abs(xa).lt.1e-15) then
          A_1a=0.0d0
        else 
          A_1a=0.5d0*(1.0d0-erf(abs(xa)/sqrt(4*deltat)))*sign(1.0d0,xa)
        endif
C---------------
        if(abs(xb).lt.1e-15) then
          A_1b=0.0d0
        else
          A_1b=0.5d0*(1.0d0-erf(abs(xb)/sqrt(4*deltat)))*sign(1.0d0,xb)
        endif
C------------------------
C23456
        W_1a=-xa**2/(2*deltat)*A_1a+xa/(2*sqrt(pi*deltat))
     1                        *exp(-xa**2/(4*deltat))
        W_1b=-xb**2/(2*deltat)*A_1b+xb/(2*sqrt(pi*deltat))
     1                        *exp(-xb**2/(4*deltat))

        W_2a=(1.0d0+xa**2/(2*deltat))*A_1a-xa/(2*sqrt(pi*deltat))
     1                           *exp(-xa**2/(4*deltat))
        W_2b=(1.0d0+xb**2/(2*deltat))*A_1b-xb/(2*sqrt(pi*deltat))
     1                           *exp(-xb**2/(4*deltat))


        dl=W_1b*mub(nt-1)+W_2b*mub(nt)-W_1a*mua(nt-1)-W_2a*mua(nt)

        end subroutine



C-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
C-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.


        subroutine form_loc(x0,x1,mu0,mu1,nt,deltat,amat,rhsa)
        implicit real*8 (a-h,o-z)
        real*8 x0,x1,mu0(1),mu1(1),deltat
        real*8 amat(2,2),rhsa(2)
        integer nt
C        discretize the local part, solve A*[mu0;mu1]=2*(u_H+D_N+g)+rhsa for mu at current time step
        done=1.0d0
        pi=4*atan(done)

        xd=x0-x1
        a1=0.5d0*(1-erf(abs(xd)/sqrt(4*deltat)))*sign(1.0d0,xd)
C23456
        w1=-xd**2/(2*deltat)*a1+xd/(2*sqrt(pi*deltat))
     1                         *exp(-xd**2/(4*deltat))
        w2=(1+xd**2/(2*deltat))*a1-xd/(2*sqrt(pi*deltat))
     1                         *exp(-xd**2/(4*deltat))
C  A=[1.0,-2*w2;-2*w2,1.0];
        amat(1,1)=1.0d0
        amat(2,2)=1.0d0
        amat(1,2)=-2.0d0*w2
        amat(2,1)=-2.0d0*w2
C  rhsa=[2*w1*mu1(nt-1);2*w1*mu0(nt-1)];
        rhsa(1)=2.0d0*w1*mu1(nt-1)
        rhsa(2)=2.0d0*w1*mu0(nt-1)

        end subroutine















