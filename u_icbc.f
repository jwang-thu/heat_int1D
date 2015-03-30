        subroutine u_init(x,res)
        implicit real*8 (a-h,o-z)
        real*8 x,res
        done=1.0d0
        pi=4*atan(done)

        res=dexp(-(x-0.4d0)**2/(4*0.01d0))/dsqrt(4*pi*0.01d0)
        end subroutine



        function u_bc0(t)
        implicit real*8 (a-h,o-z)
        real*8 t
        real*8 u_bc0
        done=1.0d0
        pi=datan(done)*4

        u_bc0=dexp(-0.4d0**2/(4.0d0*(t+0.01d0)))
     1       /dsqrt(4.0d0*pi*(t+0.01d0))
        end function



        function u_bc1(t)
        implicit real*8 (a-h,o-z)
        real*8 t
        real*8 u_bc1
        done=1.0d0
        pi=datan(done)*4

        u_bc1=dexp(-0.6d0**2/(4.0d0*(t+0.01d0)))
     1       /dsqrt(4.0d0*pi*(t+0.01d0))
        end function


