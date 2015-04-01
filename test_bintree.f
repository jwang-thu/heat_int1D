        program test_bintree
        implicit real*8 (a-h,o-z)
        real*8 a,b,eps
        integer q,nmax,len_nds,max_st,max_lv, lvf
        integer,allocatable::lnds_idst(:,:)
        real*8,allocatable::lnds_ch(:,:),lnds_cr(:,:)
        real*8,allocatable::chnodes(:)
c------------------------------------------------------
c        test the create_bintree subroutine
c         eps should >= 1d-14

C------------------------------------------------------
        a=-0.1d0
        b=1.1d0
        q=16
        eps=1d-14
        lvf=2

        nmax=2**13
c        max_lv=50
c         shouldn't be deeper into the tree
c         otherwise the endpoints of the smallest interval
c         become indistinguishable
c       try 12 levels first
        max_lv=12
c-----------------------------------

c        max_st=10000
c       this isn't a good criterion anyway...
c       use max_lv instead...
c        max_st=20

        allocate(lnds_ch(q+1,nmax))
        allocate(lnds_cr(2,nmax))
        allocate(lnds_idst(3,nmax))
        allocate(chnodes((q+1)*nmax))
        call create_bintree(a,b,q,eps,nmax,max_lv,lvf,
     1             lnds_ch,lnds_cr,lnds_idst,chnodes,len_nds)
C------------------------------------------------------
c        write(*,*) len_nds
c        write(*,*) lnds_idst(1,1:len_nds)
c        write(*,*) lnds_cr(1,1:len_nds)
c        write(*,*) 
c	write(*,*) lnds_cr(2,1:len_nds)  
c	write(*,*)  
C        write(*,*) chnodes(1:(q+1)*len_nds)
c        write(*,*) lnds_ch(1:q+1,1)

        write(*,*)  
        write(*,*) 'total points in space:', len_nds*(q+1)
        write(*,*)  
        do i=1,len_nds
          ea=lnds_cr(1,i)-lnds_cr(2,i)
          eb=lnds_cr(1,i)+lnds_cr(2,i)
          write(*,*) 'nodes', i, 
     1               '[',ea,',',eb,']',',',
     1           lnds_ch(1:q+1,i), lnds_idst(2,i), lnds_idst(3,i)
          write(*,*)  
        enddo
        
        end program
