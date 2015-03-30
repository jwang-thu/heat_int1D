# heat_int1D


Heat equation solver (1D Dirichlet problem with possible moving boundary)
(Fortran version)


1. Resolve the initial condition on a binary tree, with piecewise Chebyshev representation


files: 


create_bintree.f u_init.f chebexps.f prini.f


  subroutines: 


    create_bintree.f: (done)


        1.create_bintree (call: subdivide)


        2.subdivide (call: u_init, cheb_coeff) %%% changed cheb_coeff to chebexps


        3.cheb_coeff (call: chebexps.f/chebexps)  


  
    chebexps.f (prini.f): 
        chebexps (subroutine for chebyshev quadrature and related stuff)

    
    u_init.f: (done)
        u_init




2. Gauss transform (Chebyshev representation as input)


files: 
cfgt.f(gauss1d.f) gauss1d.f 


  subroutines: 


    cfgt.f: 


      1.find_leaf: find the leaf nodes that span the supp of Gaussian, centered at a certain point x (done)


      2.m2gauss: compute the 'moments' by recurrence relation, which is stable only for sharp gaussians (which is okay, since it is called only by direct evaluation, which is performed only when the gaussian is sharp)  (done)


      3.direct_cgt: (find_leaf, m2gauss): Direct evaluation of the Gauss Transform  (done)


      4.cfgt: (direct_cgt, gauss1d) :combination of the direct and fast parts



3. Compute the near part of double layer potential DN, given values of mu at the previous time step
files: dpotentials.f legquad.dat


  subroutines:


    dpotentials.f:


        dnear (open: legquad.dat) (done)



4. Discretize and solve/eval local part of double layer potential DL


files: dpotentials.f


  subroutines:


    dpotentials.f:


        1.form_loc (replace this with fixed point iteration) (done)


        2.eval_dl  (done)



5. A driver code that assembles these stuff together...


   heat_int1D.f


   int1D_dr.f
