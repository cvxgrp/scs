.. _matlab_example:

MATLAB
======

.. include:: qp.prob

Matlab code to solve this is below.

.. literalinclude:: qp.m
   :language: matlab

After following the matlab :ref:`install instructions <matlab_install>`, we can
run the code yielding output:

.. /Applications/MATLAB_R2021a.app/bin/matlab -nodesktop -nojvm -nosplash \
   -nodisplay -r "run ~/git/scs/docs/src/examples/qp.m;exit" > qp.m.out

.. literalinclude:: qp.m.out
   :language: none

Workspace reuse example
-----------------------

When solving many problems that share the same :code:`A` and :code:`P` but have
different :code:`b` or :code:`c`, the workspace API avoids re-factorizing:

.. code:: matlab

  % Set up the initial problem
  data.P = sparse([3., -1.; -1., 2.]);
  data.A = sparse([-1., 1.; 1., 0.; 0., 1.]);
  data.b = [-1; 0.3; -0.5];
  data.c = [-1.; -1.];
  cone.z = 1;
  cone.l = 2;
  settings = struct('verbose', 0);

  % Initialize workspace (factorizes once)
  work = scs_init(data, cone, settings);
  [x, y, s, info] = scs_solve(work);

  % Solve again with a different b
  b_new = [-0.5; 0.1; -0.3];
  scs_update(work, b_new, []);
  [x2, y2, s2, info2] = scs_solve(work);

  % Warm-start the next solve
  warm.x = x2; warm.y = y2; warm.s = s2;
  [x3, y3, s3, info3] = scs_solve(work, warm);

  % Free workspace
  scs_finish(work);
