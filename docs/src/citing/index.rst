.. _citing :

Citing SCS
===========

If you wish to cite SCS, please use any of the following:

.. glossary::

  Original paper
      Main algorithm description, derivation, and initial numerical results `paper <https://web.stanford.edu/~boyd/papers/scs.html>`__.

      .. code:: latex

        @article{ocpb:16,
            author       = {Brendan O'Donoghue and Eric Chu and Neal Parikh and Stephen Boyd},
            title        = {Conic Optimization via Operator Splitting and Homogeneous Self-Dual Embedding},
            journal      = {Journal of Optimization Theory and Applications},
            month        = {June},
            year         = {2016},
            volume       = {169},
            number       = {3},
            pages        = {1042-1068},
            url          = {http://stanford.edu/~boyd/papers/scs.html},
        }

  Latest extension
      The paper that derived the extension to quadratics and describes the latest version of the algorithm is available `here <https://arxiv.org/abs/2004.02177>`__.

      .. code:: latex

        @article{odonoghue:21,
            author       = {Brendan O'Donoghue},
            title        = {Operator Splitting for a Homogeneous Embedding of the Linear Complementarity Problem},
            journal      = {{SIAM} Journal on Optimization},
            month        = {August},
            year         = {2021},
            volume       = {31},
            issue        = {3},
            pages        = {1999-2023},
        }

  Software
     If you need to cite a particular version of the SCS software (e.g., for replication purposes) the latest version can be cited as:


      .. code:: latex

        @misc{scs,
            author       = {Brendan O'Donoghue and Eric Chu and Neal Parikh and Stephen Boyd},
            title        = {{SCS}: Splitting Conic Solver, version 3.2.9},
            howpublished = {\url{https://github.com/cvxgrp/scs}},
            month        = nov,
            year         = 2023
        }

  Anderson Acceleration
     The acceleration scheme we use is described in the `paper <https://web.stanford.edu/~boyd/papers/nonexp_global_aa1.html>`__.

      .. code:: latex

        @article{aa2020,
          title={Globally Convergent {type--I} {A}nderson Acceleration for Non-Smooth Fixed-Point Iterations},
          author={Junzi Zhang and Brendan O'Donoghue and Stephen Boyd},
          journal={{SIAM} Journal on Optimization},
          volume={30},
          number={4},
          pages={3170--3197},
          year={2020}
        }

