AC_INIT([MNI SurfReg],[0.3])
AC_CONFIG_SRCDIR([surflib/SurfaceMapBase.hpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADERS([config.h])

AC_PREREQ(2.58)

dnl Add extra include or link directories
smr_WITH_BUILD_PATH

# Require CGAL.
smr_CGAL_MAKEFILE


#------------------------------------------------------------------------
# Checks for compiler toolset.
#------------------------------------------------------------------------
AC_PROG_CC
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET

AC_PROG_LIBTOOL


#------------------------------------------------------------------------
# Look for other tools.
#------------------------------------------------------------------------

AC_PATH_PROGS(PERL, [perl5 perl], no)

if test "$PERL" = no; then
    AC_MSG_WARN(No perl program found!)
    PERL=
dnl else check perl version ....
fi

AC_PATH_PROGS(PYTHON, [python])
AC_PATH_PROGS(CREATE_TETRA_PATH, [create_tetra])
AC_PATH_PROGS(BICOBJ2VTK_PATH, [bicobj2vtk])

AC_CHECK_PROG(HAVE_DOT,dot,YES,NO)

AC_CHECK_PROG(DOXYGEN,doxygen,doxygen)
AM_CONDITIONAL(HAVE_DOXYGEN, test "$DOXYGEN" = doxygen)


#------------------------------------------------------------------------
# Check for libraries.
#------------------------------------------------------------------------

# Locate CppUnit for testing
AM_PATH_CPPUNIT(1.8.0)
AM_CONDITIONAL(HAVE_CPPUNIT, test "$CPPUNIT_LIBS")

smr_LIBS=

mni_REQUIRE_VOLUMEIO
mni_REQUIRE_BICPL
bic_LIBS=$LIBS
LIBS=$smr_LIBS

AC_PATH_XTRA
LIBS=$smr_LIBS

# GSL 1.3 has a bug in the simplex optimizer initialization.
AM_PATH_GSL([1.4])

AC_SUBST([bic_LIBS])


#------------------------------------------------------------------------
# Check that we have some required header files (and macros in them)
#------------------------------------------------------------------------
AC_HEADER_STDC
AC_CHECK_HEADERS(float.h limits.h malloc.h unistd.h)


#------------------------------------------------------------------------
# Checks for typedefs, structures, and compiler characteristics.
#------------------------------------------------------------------------
AC_C_INLINE
AC_C_CONST
AC_TYPE_SIZE_T


#------------------------------------------------------------------------
# Checks for library functions.
#------------------------------------------------------------------------
AC_FUNC_VPRINTF
AC_CHECK_FUNCS(strtod strtol)

AC_CONFIG_FILES([Makefile 
surflib/Makefile
surflib/doc/Makefile surflib/doc/Doxyfile
surflib/test/Makefile
registration/Makefile
tools/Makefile
])
AC_OUTPUT
