        subroutine heat_int1d(xx,nx,tf,x0,x1,a,b,deltat,u)

c---------------------------------------------------
c  subroutine for solving 1D heat equation.
c    fixed boundary version
c 
c  input:
c       xx: a spatial grid where the result is to be returned
c       nx: length of xx
c       tf: final time 
c       x0, x1: boundary points
c       a, b: extended boundary points
c       deltat: time step    
c output:
c       u: u(xx,tf), solution at spatial grid xx and time T   
c IC/BC:
c         initial data: u_init from u_icbc.f
c         boundary data: u_bc0, u_bc1 from u_icbc.f
c---------------------------------------------------

        implicit real*8 (a-h,o-z)
        integer nx
        real*8 xx(nx),tf,x0,x1,a,b,deltat
        real*8 u(nx)
        real*8,external::u_bc0, u_bc1

c  local variables
        integer ntmax,ntf,nt,q,maxlf
        parameter(ntmax=100000)
        parameter(maxlf=2**12)
        parameter(q=16)
        real*8 mu0(ntmax),mu1(ntmax),ep,t
c  local-var: the bintree structure
        integer len_nds, maxst, odr, len_chnds
        integer lnds_idst(2,maxlf)
        real*8 lnds_ch(q+1,maxlf),lnds_cr(2,maxlf)
        real*8 chnodes((q+1)*maxlf)
c local-var: split the bintree
        integer ind_dir(maxlf),ind_her(maxlf)
        integer len_dir, len_her

c local-var: layer potentials and such and such...
        real*8 If0, If1
c         initial potentials at x0 and x1
        real*8 dn0, dn1
c         dlp_near, at x0 and x1
        real*8 uh0, uh1
c         history part, at x0 and x1
        real*8 uh(maxlf*(q+1)), dn(maxlf*(q+1)), uht1, dnt1
c         history part and near part of dlp, inside of the domain
c         uht1 and dnt1 for testing
        real*8 amat(2,2),rhsa(2),rhs(2),ipvt(2),z(2),rcond
c         a 2-by-2 matrix solve involved, should get rid of it soon
        real*8 uh_loc, dn_loc, fval_loc(q+1), ch_loc(q+1)
        real*8 xch(q+1),uch(q+1,q+1),vch(q+1,q+1),whts(q+1)
c         for updating the chebyshev representation
        real*8 uhf(nx), dnf(nx), dlf(nx)
c         for recovering solution at final step

        done=1.0d0
        pi=datan(done)*4


c----------------------------------------------------
c set parameters...
        ntf=ceiling(tf/deltat)+1
        ep=1d-14
        maxst=10000
	odr=q

c----------------------------------------------------
c Step 0: t=0=>mu(x_i,t)=0, since BC and IC are consistent
        nt=1
        t=(nt-1)*deltat
        write(*,*) 't= ', t
        mu0(1)=0.0d0
        mu1(1)=0.0d0

c----------------------------------------------------
c Step 1: t=deltat
        nt=nt+1
        t=(nt-1)*deltat
        write(*,*) 't= ', t
        
c 1.1: resolve u_init on a binary tree by piecewise Chebyshev polynomials
c by calling create_bintree

        call  create_bintree(a,b,q,ep,maxlf,maxst,lnds_ch,
     1                    lnds_cr,lnds_idst,chnodes,len_nds)

        len_chnds=len_nds*(q+1)
c         length of all chebyshev nodes inside the domain
        write(*,*) 'number of points: ', len_chnds


c 1.2 compute the Gauss transform of u_init, a.k.a. initial potential  
        
        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,4.0d0*deltat,ep,x0,1,If0) 
        If0=If0/dsqrt(4*pi*deltat)

        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,4.0d0*deltat,ep,x1,1,If1)   
        If1=If1/dsqrt(4*pi*deltat)

c-----------------------------------------------------------------
c 1.3 discretize DL(mu), and solve for mu(x_i,t)

        rhs(1)=2.0d0*If0-2.0d0*u_bc0(t)
        rhs(2)=2.0d0*If1-2.0d0*u_bc1(t)

        call form_loc(x0,x1,mu0,mu1,nt,deltat,amat,rhsa)

        rhs(1)=rhsa(1)+rhs(1)
        rhs(2)=rhsa(2)+rhs(2)

c       x=amat\rhs...
        call dgeco(amat,2,2,ipvt,rcond,z)
        call dgesl(amat,2,2,ipvt,rhs,0)

        mu0(nt)=rhs(1)
        mu1(nt)=rhs(2)

c	write(*,*) mu0(1), mu0(2)
c	write(*,*) mu1(1), mu1(2)
c 	  the answer seems right up to now...

c----------------------------------------------------
c step 2: t=deltat*2

        nt=nt+1
        t=(nt-1)*deltat
        write(*,*) 't= ', t

c 2.1 compute the near part DN from values of previous step

        Nquad=10       
        call dnear(x0,x0,x1,deltat,mu0,mu1,nt-1,Nquad,dn0)
        call dnear(x1,x0,x1,deltat,mu0,mu1,nt-1,Nquad,dn1)

c 2.2 compute the Gauss transform of u_init, a.k.a. initial potential G_{2*deltat}(u_init)
c     it is also the history part at this time step

        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,8.0d0*deltat,ep,x0,1,uh0) 
        uh0=uh0/dsqrt(8*pi*deltat)

        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,8.0d0*deltat,ep,x1,1,uh1)   
        uh1=uh1/dsqrt(8*pi*deltat)
        
c 2.3 discretize DL(mu), and solve for mu(x_i,t) 

        rhs(1)=2.0d0*uh0+2.0d0*dn0-2.0d0*u_bc0(t)
        rhs(2)=2.0d0*uh1+2.0d0*dn1-2.0d0*u_bc1(t)

        call form_loc(x0,x1,mu0,mu1,nt,deltat,amat,rhsa)

        rhs(1)=rhsa(1)+rhs(1)
        rhs(2)=rhsa(2)+rhs(2)

c       x=amat\rhs...
        call dgeco(amat,2,2,ipvt,rcond,z)
        call dgesl(amat,2,2,ipvt,rhs,0)

        mu0(nt)=rhs(1)
        mu1(nt)=rhs(2)

c        write(*,*) 'mu0(2deltat), mu1(2deltat)=:'
c        write(*,*) mu0(nt), mu1(nt)


c 2.4  evaluate DN and UH for off boundary points
c      to be bootstraped to the next step...

        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1  chnodes,odr,a,b,8.0d0*deltat,ep,chnodes,len_chnds,uh)

        do k=1,len_chnds
          uh(k)=uh(k)/dsqrt(8.0d0*pi*deltat)
        enddo

c----testing-------
c        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
c     1  chnodes,odr,a,b,8.0d0*deltat,ep,a,1,uht1)

c        write(*,*) uht1

c        uht1=uht1/dsqrt(8.0d0*pi*deltat)

c        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
c     1  chnodes,odr,a,b,8.0d0*deltat,ep,b,1,uht2)

c        write(*,*) uht2

c        uht2=uht2/dsqrt(8.0d0*pi*deltat)

c        a dsqrt missing at first...told you
c        this could be a silly bug, nothing disastrous!
c------------------

        do k=1,len_chnds
          call dnear(chnodes(k),x0,x1,deltat,mu0,mu1,nt-1,
     1               Nquad,dn(k))
        enddo

c        write(*,*) 'test step 2.4:'
c        write(*,*) chnodes(1), chnodes(len_chnds)
c        write(*,*) uh(1), uh(len_chnds), uht1, uht2
c        write(*,*) dn(1), dn(len_chnds)
               	
c----------------------------------------------------
c step 3: t=deltat*k, k>=3
c similar to step 2, except some modification in 2.2 and 2.4
c especially 2.2

        do while(t.lt.tf)
          nt=nt+1
          t=(nt-1)*deltat
          write(*,*) 't= ', t

c 3.1 compute the near part DN from values of previous step
c     (same as step 2.1)
          Nquad=10       
          call dnear(x0,x0,x1,deltat,mu0,mu1,nt-1,Nquad,dn0)
          call dnear(x1,x0,x1,deltat,mu0,mu1,nt-1,Nquad,dn1)
c          write(*,*) dn0, dn1

c 3.2 bootstrap uh: uh=G_deltat{uh(t-deltat)+dn(mu)(t-deltat)} (bdry points)
c to do this, we should give the piecewise chebyshev representation
c of uh+dn(mu), on an binary tree
c here without coarsening of the grid, we use the same 
c binary tree as that of the initial condition u_init
c---------------------
          call chebexps(2,q+1,xch,uch,vch,whts)

          do i=1,len_nds
            do k=1,(q+1)
              uh_loc=uh((i-1)*(q+1)+k)
              dn_loc=dn((i-1)*(q+1)+k)
              fval_loc(k)=uh_loc+dn_loc
            enddo 
c               values of uh+dn, at the local chebyshev nodes

            do j=1,(q+1)
              ch_loc(j)=0.0d0
              do k=1,(q+1)
                ch_loc(j)=ch_loc(j)+uch(j,k)*fval_loc(k)
              enddo
c               convert them to chebyshev coefficients: ch_loc=uch*fval_loc
              lnds_ch(j,i)=ch_loc(j)
c               update the chebyshev coefficients
            enddo
          enddo

c---------------------
          call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,4.0d0*deltat,ep,x0,1,uh0) 
          uh0=uh0/dsqrt(4*pi*deltat)

          call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1            chnodes,odr,a,b,4.0d0*deltat,ep,x1,1,uh1)   
          uh1=uh1/dsqrt(4*pi*deltat)

c          write(*,*) 't=3deltat, uh0 and uh1:'
c          write(*,*) uh0, uh1
c               right up to here

c 3.3 discretize DL(mu), and solve for mu(x_i,t)
c     same as step 2.3

          rhs(1)=2.0d0*uh0+2.0d0*dn0-2.0d0*u_bc0(t)
          rhs(2)=2.0d0*uh1+2.0d0*dn1-2.0d0*u_bc1(t)

          call form_loc(x0,x1,mu0,mu1,nt,deltat,amat,rhsa)

          rhs(1)=rhsa(1)+rhs(1)
          rhs(2)=rhsa(2)+rhs(2)

c         x=amat\rhs...
          call dgeco(amat,2,2,ipvt,rcond,z)
          call dgesl(amat,2,2,ipvt,rhs,0)

          mu0(nt)=rhs(1)
          mu1(nt)=rhs(2)

c          write(*,*) 'mu0(3deltat), mu1(3deltat)=:'
c          write(*,*) mu0(nt), mu1(nt)


c 3.4 bootstrap uh and eval the new dn for off bdry points
c     prepare for the next step

          call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1    chnodes,odr,a,b,4.0d0*deltat,ep,chnodes,
     2    len_chnds,uh)

          do k=1,len_chnds
            uh(k)=uh(k)/dsqrt(4.0d0*pi*deltat)
          enddo   
c           the new history part (off-bdry part)

          do k=1,len_chnds
            call dnear(chnodes(k),x0,x1,deltat,mu0,mu1,
     1                 nt-1,Nquad,dn(k))
c           the new near part (off-bdry part)
          enddo     
        

c end of the do while(t<tf) loop
c------------
        enddo

c test
c          write(*,*) 'mu0(4deltat), mu1(4deltat)=:'
c          write(*,*) nt, mu0(nt), mu1(nt)


c----------------------------------------------------
c final step: recover solution u(x,t)

        call cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1    chnodes,odr,a,b,4.0d0*deltat,ep,xx,nx,uhf)        
        
        do k=1,nx
          uhf(k)=uhf(k)/dsqrt(4.0d0*pi*deltat)
        enddo 
c         the history part, evaluated at the given grid points
c---------------
        do k=1,nx
          call dnear(xx(k),x0,x1,deltat,mu0,mu1,
     1               nt-1,Nquad,dnf(k))
c         the near part, evaluated at the given grid points
        enddo
c----------------
        do k=1,nx
          call eval_dl(xx(k),nt,deltat,x0,x1,mu0,mu1,dlf(k))
        enddo
c----------------
        do k=1,nx
          u(k)=uhf(k)+dnf(k)+dlf(k)
        enddo

c done.
c----------------

        write(*,*) 'testing output: ' 
        write(*,*) u(1:5)

        
        
        
c----------------------------------------------------        
        end subroutine







