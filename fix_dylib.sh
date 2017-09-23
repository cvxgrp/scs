PYTHON_VERSION="3.5"
SO="scs-1.3.0-py3.5-macosx-10.6-x86_64.egg/_scs_indirect.cpython-35m-darwin.so"

install_name_tool -change libiomp5.dylib $CONDA_PREFIX/lib/libiomp5.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_intel_lp64.dylib $CONDA_PREFIX/lib/libmkl_intel_lp64.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_intel_thread.dylib $CONDA_PREFIX/lib/libmkl_intel_thread.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_core.dylib $CONDA_PREFIX/lib/libmkl_core.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

SO="scs-1.3.0-py3.5-macosx-10.6-x86_64.egg/_scs_direct.cpython-35m-darwin.so"

install_name_tool -change libiomp5.dylib $CONDA_PREFIX/lib/libiomp5.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_intel_lp64.dylib $CONDA_PREFIX/lib/libmkl_intel_lp64.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_intel_thread.dylib $CONDA_PREFIX/lib/libmkl_intel_thread.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

install_name_tool -change libmkl_core.dylib $CONDA_PREFIX/lib/libmkl_core.dylib $CONDA_PREFIX/lib/python${PYTHON_VERSION}/site-packages/${SO}

