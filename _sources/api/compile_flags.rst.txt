.. _compile_flags:

Compile flags
-------------

Typically the user will not need to interact with these flags, but we document them
here for completeness.  Note this list is not exhaustive, others are not exposed
but defined in :code:`include/glbopts.h`.  The list below are defined in
:code:`scs.mk` and can be overridden via the command line when compiling by
executing, e.g., :code:`make DLONG=1`, to set the :code:`DLONG` flag to True.


.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Name
     - Description
     - Permitted values
     - Default
   * - :code:`DLONG`
     - If True use 64 bit integers, else 32 bit
     - True/False
     - 0
   * - :code:`SFLOAT`
     - If True use 32 bit floats, else 64 bit (WARNING: currently broken)
     - True/False
     - 0
   * - :code:`CTRLC`
     - Listen to CTRL-C interruptions
     - True/False
     - 1
   * - :code:`NO_TIMER`
     - Disables code timing
     - True/False
     - 0
   * - :code:`NO_VALIDATE`
     - Disables data validation
     - True/False
     - 0
   * - :code:`NO_PRINTING`
     - Disables all printing in compiled binary
     - True/False
     - 0
   * - :code:`NO_READ_WRITE`
     - Disables the read/write code
     - True/False
     - 0
   * - :code:`COPYAMATRIX`
     - Make a copy of A in memory
     - True/False
     - 1
   * - :code:`GPU_TRANSPOSE_MAT`
     - If on GPU store A transpose in memory
     - True/False
     - 1
   * - :code:`VALIDATE`
     - Whether to perform problem validation or not
     - True/False
     - 1
   * - :code:`VERBOSITY`
     - Verbosity level (for debugging)
     - :math:`\mathbf{N}`
     - 0
   * - :code:`USE_LAPACK`
     - Whether to link in :ref:`BLAS/LAPACK <blas_lapack>`
     - True/False
     - 1
   * - :code:`USE_OPENMP`
     - Use openmp to parallelize some computation
     - True/False
     - 0
   * - :code:`BLAS64`
     - The BLAS library is 64 bits
     - True/False
     - 0
   * - :code:`NOBLASSUFFIX`
     - The BLAS library has no function name suffix
     - True/False
     - 0
   * - :code:`BLASSUFFIX`
     - The BLAS library uses this suffix
     - String
     - '_'
   * - :code:`MATLAB_MEX_FILE`
     - If compiling for use in MATLAB
     - True/False
     - 0
   * - :code:`PYTHON`
     - If compiling for use in python
     - True/False
     - 0
   * - :code:`USING_R`
     - If compiling for use in R
     - True/False
     - 0
