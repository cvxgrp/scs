.. _linear_solver:

Linear System Solver
====================
At each iteration SCS solves a system

.. math::

  \begin{bmatrix} 
  I + P & A^T  \\ 
  A & -I \end{bmatrix}
  \begin{bmatrix} 
  x \\ y
  \end{bmatrix} = 
  \begin{bmatrix} \mu_x \\ -\mu_y \end{bmatrix}

which is quasidefinite if we negate the bottom row:

.. math::

  \begin{bmatrix}
  r_x \\
  r_y
  \end{bmatrix}
  =
  \begin{bmatrix} 
  R_x + P  &  A^\top \\
  A &  -R_y   \\
  \end{bmatrix}
  \begin{bmatrix}
  q_x \\
  -q_y
  \end{bmatrix}

A direct method factorizes the above matrix.
An indirect method can solve via:

.. math::

  \begin{align}
  (R_x + P + A^\top R_y^{-1} A) r_x & = q_x - A^\top R_y^{-1} q_y \\
                            r_y & = R_y^{-1}(A z_x + q_y).
  \end{align}


Implementing your own linear solver
-----------------------------------
TODO
