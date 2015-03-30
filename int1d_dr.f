        program int1d_dr
        implicit real*8 (a-h,o-z)
	integer nmax
	parameter(nmax=300)
	real*8 xx(nmax), u(nmax), u_ext(nmax)
        real*8 err(nmax), err_abs, err_rel
	real*8 tf,x0,x1,a,b,deltat
        real*8,external::norm_inf

	deltat=0.0003125d0
c         2e-7 and 4e-7 when dt=0.0003125
c         and stops when dt continue to shrink
	dx=0.01d0

	x0=0.0d0
	x1=1.0d0
	a=-0.1d0
	b=1.1d0

        nx = nint((x1-x0)/dx)-1
	do i=1,nx
	  xx(i)=i*dx
	enddo



	tf=0.3d0
c         test

	call heat_int1d(xx,nx,tf,x0,x1,a,b,deltat,u)

        do k=1,nx
          call u_exact(xx(k),tf,u_ext(k))
          err(k)=u(k)-u_ext(k)
        enddo

        write(*,*) 'absolute error: ', norm_inf(err,nx)
        write(*,*) 'relative error:', 
     1         norm_inf(err,nx)/norm_inf(u_ext,nx)





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


