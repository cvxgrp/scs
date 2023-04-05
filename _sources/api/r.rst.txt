.. _r_interface:

R
=

R interface source code available `here <https://github.com/FlorianSchwendinger/scs>`_.
After :ref:`installing <r_install>` you can load SCS using

.. code:: r

  library("scs")


Usage
-----

Use function `scs` to solve a given optimization problem. Additional information
about the arguments can be found in the `R-Manual <https://cran.r-project.org/package=scs/scs.pdf>`_
or the :ref:`matrices <matrices>` of the homepage.

.. code-block:: r

    scs(A,
        b,
        obj,
        P = NULL,
        cone,
        initial = NULL,
        control = scs_control())


The `scs_control`  is used to define additional settings, additional information
can be found in the R-Manual or the :ref:`settings <settings>` section of the
homepage.

.. code-block:: r

    scs_control(max_iters = 100000L,
                eps_rel = 1e-04,
                eps_abs = 1e-04,
                eps_infeas = 1e-07,
                alpha = 1.5,
                rho_x = 1e-06,
                scale = 0.1,
                verbose = FALSE,
                normalize = TRUE,
                warm_start = FALSE,
                acceleration_lookback = 0L,
                acceleration_interval = 1L,
                adaptive_scale = TRUE,
                write_data_filename = NULL,
                log_csv_filename = NULL,
                time_limit_secs = 0)

