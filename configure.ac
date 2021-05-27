#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Smearing_radius], [0.1], [nelson.lachini@ed.ac.uk])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR([Smear_radius_singleLoad.cpp])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
m4_ifdef([AM_SILENT_RULES],[AM_SILENT_RULES([yes])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_LANG([C++])
AC_ARG_WITH([hadrons],
    [AS_HELP_STRING([--with-hadrons=<prefix>],
    [try this for a non-standard install prefix of Hadrons])],
    [PATH="$PATH$PATH_SEPARATOR$with_hadrons/bin"]
    [CXXFLAGS="-I$with_hadrons/include $CXXFLAGS"]
    [LDFLAGS="-L$with_hadrons/lib $LDFLAGS"])
AC_CHECK_PROG([HADRONSCONF],[hadrons-config],[yes])
if test x"$HADRONSCONF" != x"yes" ; then
    AC_MSG_ERROR([hadrons-config not found])
fi
AC_ARG_WITH([latan],
    [AS_HELP_STRING([--with-latan=<prefix>],
    [try this for a non-standard install prefix of LatAnalyze])],
    [PATH="$PATH$PATH_SEPARATOR$with_latan/bin"]
    [CXXFLAGS="-I$with_latan/include $CXXFLAGS"]
    [LDFLAGS="-L$with_latan/lib $LDFLAGS"])
AC_CHECK_PROG([LATANCONF],[latan-config],[yes])
if test x"$LATANCONF" != x"yes" ; then
    AC_MSG_ERROR([latan-config not found])
fi
CXXFLAGS="$CXXFLAGS `hadrons-config --cxxflags` `latan-config --cxxflags`"
LDFLAGS="$LDFLAGS `hadrons-config --ldflags` `latan-config --ldflags`"
CXXFLAGS="$AM_CXXFLAGS $CXXFLAGS -O0 -g"
LDFLAGS="$AM_LDFLAGS $LDFLAGS"
LIBS="$LIBS `hadrons-config --libs` `latan-config --libs` -lHadrons -lGrid -lLatAnalyze "
AC_MSG_CHECKING([that a minimal Grid/Hadrons program compiles]);
AC_LINK_IFELSE(
        [AC_LANG_SOURCE([[
    #include <Grid/Grid.h>
    #include <Hadrons/Application.hpp>

    using namespace Grid;
    using namespace Hadrons;

    int main(int argc, char *argv[])
    {
        Grid_init(&argc, &argv);
        Application application;
        Grid_finalize();

        return 0;
    }

    ]])],
        [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])]
    [AC_MSG_ERROR([Could not compile a minimal Grid/Hadrons program])])

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

