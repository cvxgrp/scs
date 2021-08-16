.. _linear_solver:

Linear System Solver
====================
TODO

At each iteration SCS solves a system
.. math::
  \begin{bmatrix} I + P & A^T  \\ A & -I \end{bmatrix}\begin{bmatrix} x \\ y
  \end{bmatrix} = \begin{bmatrix} \mu_x \\ -\mu_y \end{bmatrix}



Implementing your own linear solver
-----------------------------------
TODO
