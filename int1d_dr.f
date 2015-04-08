        program int1d_dr

c-------------------------------------------------
c       driver code for heat_int1D
c-------------------------------------------------
        implicit real*8 (a-h,o-z)
	integer nmax
	parameter(nmax=300)
	real*8 xx(nmax), u(nmax), u_ext(nmax)
        real*8 err(nmax), err_abs, err_rel
	real*8 x0,x1,a,b,deltat
        real*8,external::norm_inf

        real*8 tf(30)

	deltat=1.5625d-04
c              deltat=0.0003125d0
c              deltat=1.5625d-04
c              deltat=7.8125d-5
        open(unit=13,file='errors.dat')
c
c----------------------------------------------------------------------------------------
c         1.4e-8 and 2.7e-8 when dt=0.0003125 (lvf=2)
c               still good
c
c         7e-9 and 1.4e-8 when dt=1.5625d-04 (lvf=2)
c           grows after that
c       
c          force it to refine spatial grid 
c            one more time and see what happens (keep dt, lvf=3)
c          6.4e-9 and 1.3e-8 
c
c          force it to refine spatial grid 
c            one more time and see what happens (keep dt, lvf=4)
c          6.4e-9 and 1.3e-8
c                 no longer the problem of create_bintree
c
c
c       Problems:
c
c        spatial resolution:
c         0. different performance on different computers? :(
c            unable to achieve tail<1e-15 on a cims desktop
c            works on my laptop though.
c            but there's no big difference in spatial grid
c            from the tail<1e-14 case
c            won't change the result
c
c         1. criterion for 'resolving the initial data by a piecewise Chebyshev representation.
c            tail indicates error, but...
c
c         2. study the accuracy and stability of cfgt_herm and cfgt_dir
c----------------------------------------------------------------------------------------

	dx=0.01d0

	x0=0.0d0
	x1=1.0d0
	a=-0.1d0
	b=1.1d0

        nx = nint((x1-x0)/dx)-1
	do i=1,nx
	  xx(i)=i*dx
	enddo



c	tf=0.3d0
c         test
        do i=1,30
c-----------------------
          tf(i)=i*deltat       

	  call heat_int1d(xx,nx,tf(i),x0,x1,a,b,deltat,u)

          do k=1,nx
            call u_exact(xx(k),tf(i),u_ext(k))
            err(k)=u(k)-u_ext(k)
          enddo

          write(*,*) 'tf=', tf(i)
          write(*,*) 'absolute error: ', norm_inf(err,nx)
          write(*,*) 'relative error:', 
     1           norm_inf(err,nx)/norm_inf(u_ext,nx)

          write(13, *) norm_inf(err,nx),
     1         norm_inf(err,nx)/norm_inf(u_ext,nx)

c-----------------------
        enddo





        end program







c---------------------------------
        subroutine u_exact(x,t,res)
        implicit real*8 (a-h,o-z)
        real*8 x,t,res
        done=1.0d0
        pi=datan(done)*4

        res=dexp(-(x-0.4d0)**2/(4*(0.01d0+t)))
     1          /dsqrt(4*pi*(0.01d0+t))

        end subroutine

C-----------------------------------


        function norm_inf(v,n)
        implicit real*8 (a-h,o-z)
        integer n
        real*8 v(1)
        real*8 norm_inf
        
        norm_inf=dabs(v(1))
        do j=1,n
          if(dabs(v(j)) .gt. norm_inf) then
            norm_inf=dabs(v(j))
          endif
        enddo

        end function


