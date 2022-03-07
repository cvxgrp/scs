.. _c_example:

C/C++
=====

.. include:: qp.prob

C code to solve this is below.

.. literalinclude:: qp.c
   :language: c

After following the CMake :ref:`install instructions <c_install>`, we can
compile the code (assuming the library was installed in :code:`/usr/local/` and
the :code:`gcc` compiler is available) using:

.. code::

    gcc -I/usr/local/include/scs -L/usr/local/lib/ qp.c -o qp.out -lscsdir

.. ./qp.out > qp.c.out

Then running the binary yields output:

.. literalinclude:: qp.c.out
   :language: none

