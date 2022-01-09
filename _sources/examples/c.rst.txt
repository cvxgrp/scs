.. _c_example:

C/C++
=====

.. include:: qp.prob

C code to solve this is below.

.. literalinclude:: qp.c
   :language: c

After following the CMake :ref:`install instructions <c_install>`, we can
compile the code using:

.. code::

	  gcc -I/usr/local/include/scs -L/usr/local/lib/ -lscsdir qp.c -o qp.out

.. ./qp.out > qp.c.out

Then run the code yielding output

.. literalinclude:: qp.c.out
   :language: none

