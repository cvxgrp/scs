.. _r_example:

R
=

.. code:: r

  library("scs")


Second-order cone programming
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
  \begin{array}{rr}
  \underset{x}{\text{maximize}}
    & x + y  \\
  \text{subject to}
    & \sqrt{x^2 + y^2} \leq \sqrt{2} \\
    & x, y \geq 0
  \end{array}


.. code:: r
  
  obj <- c(-1, -1)
  A <- rbind(c(-1,  0),  # x >= 0
             c( 0, -1),  # y >= 0
             c( 0,  0),
             c(-1,  0),
             c( 0, -1))
  b <- c(0, 0, sqrt(2), 0, 0)
  cone <- list(z = 0, l = 2, q = 3)

  solution <- scs(A, b, obj, cone = cone)
  str(solution)
  #R> List of 4
  #R>  $ x   : num [1:2] 1 1
  #R>  $ y   : num [1:5] 0 0 1.41 -1 -1
  #R>  $ s   : num [1:5] 1 1 1.41 1 1
  #R>  $ info:List of 21
  #R>   ..$ iter                : int 50
  #R>   ..$ status              : chr "solved"
  #R>   ..$ status_val          : int 1
  #R>   ..$ scale_updates       : int 0
  #R>   ..$ pobj                : num -2
  #R>   ..$ dobj                : num -2
  #R>   ..$ res_pri             : num 4.86e-05
  #R>   ..$ res_dual            : num 1.98e-06
  #R>   ..$ gap                 : num 0.000141
  #R>   ..$ res_infeas          : num NaN
  #R>   ..$ res_unbdd_a         : num 0.707
  #R>   ..$ res_unbdd_p         : num 0
  #R>   ..$ setup_time          : num 0.0236
  #R>   ..$ solve_time          : num 0.0237
  #R>   ..$ scale               : num 0.1
  #R>   ..$ comp_slack          : num 1.34e-14
  #R>   ..$ rejected_accel_steps: int 0
  #R>   ..$ accepted_accel_steps: int 0
  #R>   ..$ lin_sys_time        : num 0.00668
  #R>   ..$ cone_time           : num 0.00439
  #R>   ..$ accel_time          : num 0



Primal exponential cone
^^^^^^^^^^^^^^^^^^^^^^^

.. math::
  \begin{array}{rr}
  \underset{x}{\text{maximize}}  & x + y + z \\
  \text{subject to} & y   e^{\frac{x}{y}} \leq z \\
  & x \geq 0, y > 0, z \in [0, e]
  \end{array}


.. code:: r

  obj <- c(-1, -1, -1)
  A <- rbind(-diag(3),   # x, z >= 0, y >= 1e-12
             c(0, 0, 1), # z <= e
             -diag(3))   # K_exp
  b <- c(0, -1e-12, 0, exp(1), double(3))
  cone <- list(z = 0L, l = 4L, ep = 1L)
  solution <- scs(A, b, obj, cone = cone)
  str(solution)
  #R> List of 4
  #R>  $ x   : num [1:3] 8.46e-05 2.72 2.72
  #R>  $ y   : num [1:7] 1.17e-05 0.00 0.00 2.00 -1.00 ...
  #R>  $ s   : num [1:7] 0.00 2.72 2.72 0.00 5.68e-05 ...
  #R>  $ info:List of 21
  #R>   ..$ iter                : int 75
  #R>   ..$ status              : chr "solved"
  #R>   ..$ status_val          : int 1
  #R>   ..$ scale_updates       : int 0
  #R>   ..$ pobj                : num -5.44
  #R>   ..$ dobj                : num -5.44
  #R>   ..$ res_pri             : num 0.000273
  #R>   ..$ res_dual            : num 2.02e-05
  #R>   ..$ gap                 : num 0.000583
  #R>   ..$ res_infeas          : num NaN
  #R>   ..$ res_unbdd_a         : num 0.5
  #R>   ..$ res_unbdd_p         : num 0
  #R>   ..$ setup_time          : num 0.0273
  #R>   ..$ solve_time          : num 0.273
  #R>   ..$ scale               : num 0.1
  #R>   ..$ comp_slack          : num 8.54e-09
  #R>   ..$ rejected_accel_steps: int 0
  #R>   ..$ accepted_accel_steps: int 0
  #R>   ..$ lin_sys_time        : num 0.0118
  #R>   ..$ cone_time           : num 0.242
  #R>   ..$ accel_time          : num 0
