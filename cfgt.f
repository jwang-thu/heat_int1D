        
        subroutine cfgt(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1                  chnodes,odr,a,b,d,tol,xtarg,ntarg,gx)
C---------------------------------
C       hybrid chebyshev fast gauss transform : int_a^b f(y)exp(-(x-y)^2/d) dy,
c         using direct evaluation, if the sub-interval is large compared to the supp of gaussian
c               and hermite expansion fast algorithm, if the sub-interval is small
C input: 
C     lnds_ch, lnds_cr, lnds_idst: real/int arrays, binary tree leaf nodes structure
C                                  the output of create_bintree
C     len_nds: integer, length of leaf_nodes, output of create_bintree
C     odr: integer, order of local Chebyshev coefficients, q=odr+1
C        ---pay attention to this minor inconsistency of q and odr
C     a,b: real, end points of the interval of computation
C     d: real, delta in the Gauss Transform
c     xtarg: real array, target points where the cfgt is to be computed
c     ntarg: number of target points
C
C output:
C     gx: real array of length ntarg, the values of the gauss transform
C---------------------------------

        implicit real*8 (a-h,o-z)
        integer len_nds,nmax,odr,q,ntarg
        parameter(nmax=2**12)
        integer lnds_idst(2,len_nds)
        real*8 lnds_ch(odr+1,len_nds),lnds_cr(2,len_nds)        
        real*8 a,b,d,tol
        real*8 chnodes((odr+1)*len_nds)
        real*8 xtarg(ntarg),gx(ntarg),gx_dir(ntarg)


        integer ind_dir(len_nds),ind_her(len_nds)
        integer len_dir, len_her
c         for splitting of bintree

        real*8 ends(2),mx(odr+1)
        integer ind_tmp(2),len_ind,ind(len_nds)
        logical flg1(2),flg2(2),flg(2)
c         for finding leaf nodes that cover the supp of gaussian

C---------------------------------
        done=1.0d0
        pi=4*atan(done)

        q=odr+1
        rg=6.0d0
        sd=sqrt(d)
        rs=rg*sd
C         some parameters...see the cfgt paper for reference
c        gx=0.0d0

        ends(1)=x-rs
        ends(2)=x+rs
C         support of the Gaussian
C----------------------------------

c       split the leaf nodes into two parts
c       according to their radii and rs
        call split_bintree(lnds_cr,len_nds,rs,ind_dir,len_dir
     1                          ,ind_her,len_her)

C----------------------------------

c       hermite expansion part
c         assemble all the chebyshev nodes into xat
c         compute str
c         and call gausst....
        call cfgt_herm(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1                       ind_her,len_her,chnodes,
     2                       odr,a,b,d,tol,xtarg,ntarg,gx)



C----------------------------------

c       the direct evaluation part
c         as usual...
c         consider the large leaf nodes only
        do i=1,ntarg
          call cfgt_dir(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1    ind_dir,len_dir,odr,nmax,a,b,xtarg(i),d,gx_dir(i))
        enddo

C----------------------------------

        do i=1,ntarg
          gx(i)=gx(i)+gx_dir(i)
        enddo
c         add up two parts

     
        end subroutine



C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- 
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


       

        subroutine cfgt_herm(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1                       ind_her,len_her,chnodes,
     2                       odr,a,b,d,tol,xtarg,ntarg,gx)
c------------------------------------------
c Hermite expansion part of CFGT
c   used when the supp of gaussian is large
c   compared to the size of the leaf node
c
c input: ind_her: local (within all leaf nodes) indices
c                 of hermite expansion leaves
c                 output of split_bintree
c------------------------------------------

        implicit real*8 (a-h,o-z)
        integer len_nds, len_her, odr, ntarg
        real*8 lnds_ch(odr+1,len_nds), lnds_cr(2,len_nds)
        integer lnds_idst(2,len_nds), ind_her(len_nds)
        real*8 chnodes((odr+1)*len_nds)
        real*8 a,b,d,tol,xtarg(ntarg),gx(ntarg)
c------------------------------------------
c       local variables
        integer natoms
        integer maxwrk, maxatm
        parameter(maxwrk=100000000)
        parameter(maxatm=1000000)

        integer iwork(maxwrk), lenw
        integer iout1, iout2,inform(6)

        real*8 xat((odr+1)*len_her), str((odr+1)*len_her)
        real*8 xats((odr+1)*len_her), xtargs(ntarg), ds
c              xat, xtarg, and d scaled to [0,1]
        real*8 work(maxwrk)

        real*8 x(odr+1),u(odr+1,odr+1),v(odr+1,odr+1)
        real*8 whts(odr+1)
c              for chebexps

        natoms=(odr+1)*len_her
c              number of sources (atoms, as they call them)
        lenw=maxwrk
c------------------------------------------
c       first, go through the ind_her list, assemble their chebnodes
c              and compute str(i)=f(y_i)*w_i
c              as a discretization step
        call chebexps(2,odr+1,x,u,v,whts)
c            call this guy once...

        do ip=1,len_her
c-------------------------
c          might have to consume a lot of indices, so... the name 'ip'
          k=ind_her(ip)
c          write(*,*) 'k= ',k
c          write(*,*) 'r=', lnds_cr(2,k)

c          the ip-th hermite node, which is the k-th leaf node within all leaf nodes

          do j=1,(odr+1)
            xat((ip-1)*(odr+1)+j)=chnodes((k-1)*(odr+1)+j)
c           assemble the chebyshev nodes...
c           looks right....
          enddo

c         convert the chebyshev coefficients 
c         to corresponding function values

c         f=v*chebcoeff
          do j=1,(odr+1)
            str((ip-1)*(odr+1)+j)=0.0d0
            do i=1,(odr+1)
              str((ip-1)*(odr+1)+j)=str((ip-1)*(odr+1)+j)+
     1            v(j,i)*lnds_ch(i,k)
            enddo
c             str should be the values of f at chebyshev nodes now...
c            let's test this part alone...seems right...

            str((ip-1)*(odr+1)+j)=str((ip-1)*(odr+1)+j)*whts(j)
     1                           *lnds_cr(2,k)
c       attention: don't miss the factor of lnds_cr(2,k) here!!
c       the weights are for the interval [-1,1]
c       here the interval is [c-r,c+r]

          enddo
c             str should be the real str now...
c-------------------------
        enddo


c------------------------------------------
c       second, call gausst (in gauss1d.f) to compute 
c       the discrete fgt for this part...
c       remember to scale everything to [0,1]...

c..........scale everything to [0,1]:
        ds=d/((b-a)**2)

        do i=1,natoms
          xats(i)=(xat(i)-a)/(b-a)
        enddo
        
        do i=1,ntarg
          xtargs(i)=(xtarg(i)-a)/(b-a)
        enddo
c...............................

c       call gausst to evaluate the fast gauss transform
        iout1=0
        iout2=0

        call gausst(xats,natoms,str,ds,gx,iwork,
     1              work,lenw,xtargs,ntarg,tol,
     2              iout1,iout2,inform)


c------------------------------------------
        end subroutine






C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- 





        subroutine cfgt_dir(lnds_ch,lnds_cr,lnds_idst,len_nds,
     1                    ind_dir,len_dir,odr,nmax,a,b,x,d,gx)
C---------------------------------
C       the direct evaluation part of cfgt : int_a^b f(y)exp(-(x-y)^2/d) dy
C input: 
C     lnds_ch, lnds_cr, lnds_idst: real/int arrays, binary tree leaf nodes structure
C                                  the output of create_bintree
C     len_nds: integer, length of leaf_nodes, output of create_bintree
c     len_dir: integer, length of ind_dir, the direct part of leaf_nodes,
c                       output of split_bintree
C     odr: integer, order of local Chebyshev coefficients, q=odr+1
C        ---pay attention to this minor inconsistency of q and odr
C     a,b: real, end points of the interval of computation
C     x: real, the target points where the cfgt is to be computed
C     d: real, delta in the Gauss Transform
C output:
C     val: the value of the gauss transform
C---------------------------------
        implicit real*8 (a-h,o-z)
        integer len_nds,odr,q,nmax,len_dir
        real*8 lnds_ch(odr+1,nmax),lnds_cr(2,nmax)        
        real*8 a,b,x,d,gx
        integer lnds_idst(2,nmax)
        integer ind_dir(len_dir)

        real*8 ends(2),mx(odr+1)
        integer ind_tmp(2),len_ind,ind(len_nds),n1,n2
        logical flg1(2),flg2(2),flg(2)
C---------------------------------
        done=1.0d0
        pi=4*atan(done)

        q=odr+1
        rg=6.0d0
        sd=sqrt(d)
        rs=rg*sd
C         some parameters...see the cfgt paper for reference
        gx=0.0d0
C----------------------------------
        ends(1)=x-rs
        ends(2)=x+rs
C         support of the Gaussian
C         now find the leaf nodes that span the supp 
C         pay attention to special cases...

        do k=1,2
          call find_leaf(lnds_cr,len_nds,a,b,ends(k),ind_tmp(k))
          flg1(k)=(ind_tmp(k).gt.0)
          flg2(k)=(ind_tmp(k).le. len_nds)
        enddo

C        post processing of the results
        flg=flg1 .and. flg2
        if (flg(1).and.flg(2)) then
c           both lie in [a,b]
          n1=ind_tmp(1)
          n2=ind_tmp(2)
          call sift_dir(n1,n2,len_dir,ind_dir,len_ind,ind)
        elseif (.not.(flg1(1).or.flg1(2))) then
C           both left of [a,b]
          len_ind=0
        elseif (.not.(flg2(1).or.flg2(2))) then
C           both right of [a,b]
          len_ind=0
        elseif ((.not. flg1(1)).and.(flg2(2))) then
C           left end out, right end in
          n1=1
          n2=ind_tmp(2)
          call sift_dir(n1,n2,len_dir,ind_dir,len_ind,ind)
        elseif ((flg1(1)).and.(.not. flg2(2))) then
C           left end in, right end out
          n1=ind_tmp(1)
          n2=len_nds
          call sift_dir(n1,n2,len_dir,ind_dir,len_ind,ind)
        else
C         left end left, right end right
          n1=1
          n2=len_nds
          call sift_dir(n1,n2,len_dir,ind_dir,len_ind,ind)
        endif

C       ind is the local indices (within lnds arrays) of the big (dir) leaf nodes that 
c       lie in the supp of the gaussian

C----------------------------------
C       testing...
C        write(*,*) len_ind
C        write(*,*) ind_tmp
C        write(*,*) flg1
C        write(*,*) flg2
C        if (len_ind.gt.0) then
C          write(*,*) ind(1:len_ind)
C        endif 
C----------------------------------
        if (len_ind .gt. 0) then
          do k=1, len_ind
            c=lnds_cr(1,ind(k))
            r=lnds_cr(2,ind(k))
            xmin = max(-rg, (c-r-x)/sd)
            xmax = min(rg, (c+r-x)/sd)
            call m2gauss(q,xmin,xmax,sd/r,(x-c)/r,mx)
c            write(*,*) sd/r
C               mx here should be an array of length q
            gxk=0.0d0
            do j=1,q
              gxk=gxk+mx(j)*lnds_ch(j,ind(k))
            enddo
C            write(*,*) ind(k)
C            write(*,*) lnds_ch(1:q,ind(k))
            gx=gx+gxk
          enddo
          gx=gx*sd
        endif
C----------------------------------
     
        end subroutine




C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- 
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- 


        function chebypol(n,x)
C         chebyshev polynomials
        implicit real*8 (a-h,o-z)
        integer n
        real*8 x,sig
        real*8 chebypol

c        chebypol=real(cos(n*acos(x)))
c         handle the case with abs(x)>1
c
c        attention: dacosh might not be defined
c                   for any compiler...
        if (abs(x) .le. 1.0d0) then
          chebypol=dcos(n*dacos(x))
        elseif (x .gt. 1.0d0) then
c          write(*,*) 'attention-chebypol: x>1.0 !', 'x=', x
cccccccccccc
          chebypol=dcosh(n*dacosh(x))
        else
c          write(*,*) 'attention-chebypol: x<-1.0 !', 'x=', x
cccccccccccc
c.............
          if (mod(n,2) .eq. 0) then
            sig=1.0d0
          else
            sig=-1.0d0
          endif
c..............
         chebypol=sig*dcosh(n*dacosh(-x))
        endif
        end function


C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-


        function const_2gauss(n, a, b, d, e)
C        some constant used in m2gauss
        implicit real*8 (a-h,o-z)
        integer n
        real*8 a,b,d,e
        real*8 const_2gauss

        const_2gauss=0.5d0/d*(exp(-b**2)*(chebypol(n+1,d*b+e)/(n+1)
     1              -chebypol(n-1, d*b+e)/(n-1)) -
     2    exp(-a**2)*(chebypol(n+1, d*a+e)/(n+1) - 
     3    chebypol (n-1, d*a+e)/(n-1)))

        end function



C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.- 
 
        

        subroutine m2gauss(q, a, b, d, e, mx)
C         precompute the 'moments'
c               --add some stablization anyway...just to make it robust...
        implicit real*8 (a-h,o-z)
        integer q
        real*8 a,b,d,e,mx(q)
        real*8,external::chebypol,const_2gauss

        real*8 expa, expb

        done=1.0d0
        pi=4*atan(done)

        expa=exp(-a**2)
        expb=exp(-b**2)
        
        mx(1)=0.5d0*sqrt(pi)*(erf(b)-erf(a))
        mx(2)=0.5d0*d*(expa-expb)+e*mx(1)
        mx(3)=0.5d0*expb*(-2*b*d**2-4*e*d)+0.5d0*expa*(2*a*d**2+4*d*e)+
     1       (-1.0d0+d**2+2*e**2)*mx(1)

        mx(4)=-d*(expb*chebypol(2,d*b+e)-expa*chebypol(2,d*a+e))
     1        +2*e*mx(3)-(1.0d0-4*d**2)*mx(2)

C       --------base cases-------
        if (q .gt. 5) then
          do k=5,q
            n=k-1
            mx(k)=2*e*mx(k-1)+2*((n-1)*d**2+1.0d0/(n-3))*mx(k-2)
     1          -2.0d0*e*(n-1)/(n-3)*mx(k-3)+(n-1)*1.0d0/(n-3)*mx(k-4)
     2          -2*d**2*(n-1)*const_2gauss(n-2,a,b,d,e)
          enddo
        endif
C    --------recurrence relation------
c    attention: this recurrence is unstable when d>1
c               avoid using this when d>1
        
        end subroutine




C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-     




        subroutine find_leaf(lnds_cr,leng,a,b,x,ind_loc)
C         find the leaf node that x lies in
C         if x<a return 0
C         if x>b return L+1, L=length(leaf_nodes)
C         return -1 if not found, and that should never happen
        implicit real*8 (a-h,o-z)
        integer leng,ind_loc
        real*8 lnds_cr(2,1),a,b,x

        real*8 right(leng)
C------------------------------------------
        l=leng
        do j=1,l
          right(j)=lnds_cr(1,j)+lnds_cr(2,j)
C          right end-point of each sub-interval
        enddo
        ind_loc=-1
C-----------------------------
        if (x .lt. a) then
          ind_loc=0
        elseif (x .gt. b) then
          ind_loc=l+1
        else
          ipt=1
          do while(ipt .le. l)
            if (right(ipt) .ge. x) then
              ind_loc=ipt
              ipt=l+1
C              the leaf node found, exit loop by set ipt large
            endif

            ipt=ipt+1
C              continue the loop
          enddo

        endif
C-----------------------------
        if (ind_loc .lt. 0) then
          write(*,*) 'error finding the leaf node that x lies in!'
        endif
C------------------------------

        end subroutine




C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
C.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-



        subroutine sift_dir(n1,n2,len_dir,ind_dir,len_ind,ind)

        implicit real*8 (a-h,o-z)
        integer n1, n2, len_dir, len_ind
        integer ind_dir(len_dir), ind(len_ind)

        integer nt1, nt2, pt1, pt2

c       sift ind_dir, find all those nodes that is in [n1, n2]
c       output: len_ind, ind
        pt1=1
        pt2=len_dir

        nt1=ind_dir(pt1)
        nt2=ind_dir(pt2)


        if ((nt1.gt.n2) .or. (nt2.lt.n1)) then
          len_ind=0
c         both left or right out of [n1, n2]        
        else
          do while(nt1 .lt. n1)
            pt1=pt1+1
            if (pt1 .le. len_dir) then
              nt1=ind_dir(pt1)
            endif
          enddo
c           find the first ind_dir that is >= n1

          do while(nt2 .gt. n2)
            pt2=pt2-1
            if (pt2 .ge. 1) then
              nt2=ind_dir(pt2)
            endif
          enddo   
c           find the left most ind_dir that is <= n2

          if ( (pt1 .le. len_dir).and.(pt2 .ge. 1)
     1        .and.(pt1 .le. pt2) ) then
            len_ind=pt2-pt1+1
            do j=1,len_ind
              ind(j)=ind_dir(pt1+j-1)
            enddo
          else
            len_ind=0
          endif 
c           assemble all the eligible nodes into ind  
        endif  
        
        
        end subroutine


























