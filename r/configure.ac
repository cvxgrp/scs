dnl R interface to SCS

AC_INIT([scs], 1.1.8)

SCS_FLAGS="-DLAPACK_LIB_FOUND -DCOPYAMATRIX"
SCS_LIBS="-lm"
SCS_OUT_LIB="out/libscsdir.a"

PKG_CFLAGS="-I../../include"
PKG_LIBS="-L../../out -lscsdir"

dnl Find the R home directory.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "Could not determine R_HOME."
  exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`

LAPACK_LIBS=`"${R_HOME}/bin/R" CMD config LAPACK_LIBS`
BLAS_LIBS=`"${R_HOME}/bin/R" CMD config BLAS_LIBS`
dnl FLIBS=`"${R_HOME}/bin/R" CMD config FLIBS`

AC_PROG_CC

PKG_LIBS="${PKG_LIBS}" 

CFLAGS="${CFLAGS} ${PKG_CFLAGS}"

dnl substitute in src/Makevars.in
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(BLAS_LIBS)
AC_SUBST(LAPACK_LIBS)
AC_SUBST(SCS_FLAGS)
AC_SUBST(SCS_LIBS)
AC_SUBST(SCS_OUT_LIB)
AC_OUTPUT(src/Makevars)

exit 0

