dnl Process this file with autoconf to produce a configure script

dnl ------------------------------------------------
dnl Initialization and Versioning
dnl ------------------------------------------------

AC_INIT([libtheora],[1.2.0],[theora-dev@xiph.org])

AC_CANONICAL_HOST

AC_CONFIG_SRCDIR([lib/fdct.c])
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE([enable])

dnl we use doc_DATA in doc/Makefile.am which requires autoconf >= 2.60
dnl to define docdir for us.
AC_PREREQ([2.71])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

################################################################################
# Set the shared versioning info, according to section 6.3 of the libtool info #
# pages. CURRENT:REVISION:AGE must be updated immediately before each release: #
#                                                                              #
#   * If the library source code has changed at all since the last             #
#     update, then increment TH*_LIB_REVISION (`C:R:A' becomes `C:r+1:A').     #
#                                                                              #
#   * If any interfaces have been added, removed, or changed since the         #
#     last update, increment TH*_LIB_CURRENT, and set TH*_LIB_REVISION to 0.   #
#                                                                              #
#   * If any interfaces have been added since the last public release,         #
#     then increment TH*_LIB_AGE.                                              #
#                                                                              #
#   * If any interfaces have been removed since the last public release,       #
#     then set TH*_LIB_AGE to 0.                                               #
#                                                                              #
################################################################################
TH_LIB_CURRENT=5
TH_LIB_REVISION=1
TH_LIB_AGE=4
AC_SUBST(TH_LIB_CURRENT)
AC_SUBST(TH_LIB_REVISION)
AC_SUBST(TH_LIB_AGE)

THDEC_LIB_CURRENT=3
THDEC_LIB_REVISION=1
THDEC_LIB_AGE=1
AC_SUBST(THDEC_LIB_CURRENT)
AC_SUBST(THDEC_LIB_REVISION)
AC_SUBST(THDEC_LIB_AGE)

THENC_LIB_CURRENT=4
THENC_LIB_REVISION=1
THENC_LIB_AGE=2
AC_SUBST(THENC_LIB_CURRENT)
AC_SUBST(THENC_LIB_REVISION)
AC_SUBST(THENC_LIB_AGE)

dnl Extra linker options (for version script)
THEORA_LDFLAGS=""

dnl --------------------------------------------------  
dnl Check for programs
dnl --------------------------------------------------  

AM_PROG_AS

dnl save $CFLAGS since AC_PROG_CC likes to insert "-g -O2"
dnl if $CFLAGS is blank
cflags_save="$CFLAGS"
AC_PROG_CC
AC_PROG_CPP
CFLAGS="$cflags_save"

AM_PROG_CC_C_O
LT_INIT([win32-dll])

dnl Add parameters for aclocal
AC_CONFIG_MACRO_DIR([m4])

dnl Check for doxygen
AC_ARG_ENABLE([doc],
	AS_HELP_STRING([--disable-doc], [Do not build API documentation]),
	[ac_enable_doc=$enableval], [ac_enable_doc=auto])

if test "x$ac_enable_doc" != "xno"; then
	AC_CHECK_PROG(HAVE_DOXYGEN, doxygen, true, false)

	if test "x$HAVE_DOXYGEN" = "xfalse" -a "x$ac_enable_doc" = "xyes"; then
		AC_MSG_ERROR([*** API documentation explicitly requested but Doxygen not found])
	fi
else
	HAVE_DOXYGEN=false
fi
AM_CONDITIONAL(HAVE_DOXYGEN,$HAVE_DOXYGEN)
if test $HAVE_DOXYGEN = "false"; then
        AC_MSG_WARN([*** doxygen not found, API documentation will not be built])
fi

dnl Check for tools used to build the format specification
BUILD_SPEC="false"
ac_build_spec=yes
AC_ARG_ENABLE(spec,
     AS_HELP_STRING([--disable-spec], [Do not build the specification]),
     [
       if test "x$enableval" = "xno"; then
         ac_build_spec=$enableval
       fi
     ], [
       ac_build_spec=yes
     ] )
if test "x$ac_build_spec" = "xyes"; then
  AC_CHECK_PROG(HAVE_PDFLATEX, pdflatex, yes)
  AC_CHECK_PROG(HAVE_BIBTEX, bibtex, yes)
  AC_CHECK_PROG(HAVE_TRANSFIG, fig2dev, yes)
  AC_MSG_NOTICE([Checking for packages in ${srcdir}/doc/spec/spec.tex...])
  if test -r ${srcdir}/doc/spec/spec.tex; then
   if test "x$HAVE_PDFLATEX" = "xyes"; then
    if test "x$HAVE_BIBTEX" = "xyes"; then
      if test "x$HAVE_TRANSFIG" = "xyes"; then
        tex_pkg_list=`fgrep usepackage ${srcdir}/doc/spec/spec.tex | grep \{ | grep -v ltablex`
        tex_pkg_ok="yes"
        for pkg_line in $tex_pkg_list; do
	  pkg_name=`echo $pkg_line | sed -e 's/.*{\(.*\)}.*/\1/'`
	  AC_MSG_CHECKING([for Tex package $pkg_name])
          cat >conftest.tex <<_ACEOF
\\documentclass{book}
$pkg_line
\\begin{document}
Hello World.
\\end{document}
_ACEOF
	  if pdflatex -interaction batchmode -halt-on-error conftest < /dev/null > /dev/null 2>&1; then
            AC_MSG_RESULT([ok])
          else
	    tex_pkg_ok="no"
	    AC_MSG_RESULT([no])
          fi
        done
        if test -w conftest.tex; then rm conftest.tex; fi
        if test -w conftest.tex; then rm conftest.aux; fi
        if test -w conftest.pdf; then rm conftest.pdf; fi
        if test "x$tex_pkg_ok" = "xyes" && test x$cross_compiling = xno; then
          BUILD_SPEC="true"
        fi
      fi
    fi
   fi
  fi
fi
AM_CONDITIONAL(BUILD_SPEC, $BUILD_SPEC)
if test $BUILD_SPEC = "false"; then
	AC_MSG_WARN([*** Format Specification will not built.])
fi

# Check for valgrind
AC_ARG_ENABLE(valgrind-testing,
     AS_HELP_STRING([--enable-valgrind-testing],[enable running of tests inside Valgrind (default disabled)]))
if test "x$enable_valgrind_testing" = "xyes"
then
  AC_PATH_PROG(VALGRIND,valgrind,)
  if test "x$VALGRIND" != "x"
  then
    VALGRIND="$VALGRIND -q --error-exitcode=99 --leak-check=full --show-reachable=yes --num-callers=50"
    TESTS_INFO="Test suite will be run under:
    ${VALGRIND}"
  else
    TESTS_INFO="Type 'make check' to run test suite (Valgrind not found)"
  fi
else
  TESTS_INFO="Type 'make check' to run test suite (Valgrind testing not enabled)"
fi

dnl --------------------------------------------------
dnl Set build flags based on environment
dnl --------------------------------------------------

dnl Set some host options

cc_compiler=unknown
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#if ! __clang__
                                    #error
                                    #endif
                                    ]])],
                  cc_compiler=clang,
                  [AC_COMPILE_IFELSE([AC_LANG_SOURCE([[#if ! __GNUC__
                                                       #error
                                                       #endif
                                                       ]])],
                                     cc_compiler=gcc, [])])

cflags_save="$CFLAGS"
if test $cc_compiler != "gcc" ; then
        case $host in 
        *)
                DEBUG="-g -DDEBUG"
                CFLAGS="-O"
                PROFILE="-g -p -DDEBUG" ;;
        esac
else
        case $host in 
        *)
                DEBUG="-g -Wall -Werror=uninitialized -Winit-self -Wno-parentheses -DDEBUG -D__NO_MATH_INLINES"
                CFLAGS="-Wall -Werror=uninitialized -Winit-self -Wno-parentheses -O3 -fomit-frame-pointer -finline-functions -funroll-loops"
                PROFILE="-Wall -Werror=uninitialized -Winit-self -Wno-parentheses -pg -g -O3 -fno-inline-functions -DDEBUG";;
        esac
fi
CFLAGS="$CFLAGS $cflags_save"

AC_ARG_ENABLE(gcc-sanitizers,
    AS_HELP_STRING([--enable-gcc-sanitizers], [Enable GCC sanitizers]),
    [ ac_enable_gcc_sanitizers=$enableval ], [ ac_enable_gcc_sanitizers=no] )

if test $cc_compiler = "gcc" && test "x${ac_enable_gcc_sanitizers}" = xyes; then
  CFLAGS="${CFLAGS} -fsanitize=address -fsanitize=undefined -g"
  LDFLAGS="${CFLAGS} -fsanitize=address"

  TEST_ENV="env UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1"
  AC_SUBST(TEST_ENV)
fi

cpu_x86_64=no
cpu_x86_32=no
cpu_arm=no
cpu_c64x=no
AC_ARG_ENABLE(asm,
    AS_HELP_STRING([--disable-asm], [Disable assembly optimizations]),
    [ ac_enable_asm=$enableval ], [ ac_enable_asm=yes] )

if test "x${ac_enable_asm}" = xyes; then
  cpu_optimization="no optimization for your platform, please send a patch"
  case $host_cpu in
  i[[3456]]86)
    cpu_x86_32=yes
    cpu_optimization="32 bit x86"
    AC_DEFINE([OC_X86_ASM], [],  [make use of x86 asm optimization])
    if test "x$host_vendor" = "xapple"; then
      THEORA_LDFLAGS="$THEORA_LDFLAGS  -Wl,-read_only_relocs,suppress"
    fi
      ;;
  x86_64)
    cpu_x86_64=yes
    cpu_optimization="64 bit x86"
    AC_DEFINE([OC_X86_ASM], [],  [make use of x86 asm optimization])
    AC_DEFINE([OC_X86_64_ASM], [],  [make use of x86_64 asm optimization])
    ;;
  arm*|aarch*)
    cpu_arm=yes
    cpu_optimization="ARM"
    AC_DEFINE([OC_ARM_ASM], [],  [make use of arm asm optimization])
    AC_CHECK_DECLS([vabal_u8], [HAVE_ARM_ASM_NEON=1], [HAVE_ARM_ASM_NEON=0], [[#include <arm_neon.h>]])

    if test x$HAVE_ARM_ASM_NEON = x1 ; then
      AC_DEFINE(OC_ARM_ASM_NEON, 1,
        [Define if compiler supports NEON instructions])
      cpu_optimization="$cpu_optimization (NEON)"
    fi
    AC_SUBST(HAVE_ARM_ASM_NEON)
    ;;
  tic6x)
    cpu_c64x=yes
    cpu_optimization="TI C64x+"
    AC_DEFINE([OC_C64X_ASM], [],  [make use of c64x+ asm optimization])
    ;;
  wasm*)
    cpu_wasm=yes
    cpu_optimization="WASM"
    AC_ARG_ENABLE(wasm-simd128,
      AS_HELP_STRING([--enable-wasm-simd128], [Enable WASM SIMD128 decoder optimizations]),
      [ ac_enable_wasm_simd128=$enableval ], [ ac_enable_wasm_simd128=no ])

    if test x$ac_enable_wasm_simd128 = xyes ; then
      CFLAGS="$CFLAGS -msimd128"
      AC_DEFINE(OC_WASM_SIMD128, 1,
        [Define if targeting WASM with SIMD128 support])
      cpu_optimization="$cpu_optimization (SIMD128)"
    fi
    AC_SUBST(OC_WASM_SIMD128)
    ;;
  esac
else
  cpu_optimization="disabled"
fi
AM_CONDITIONAL([CPU_x86_64], [test x$cpu_x86_64 = xyes])
AM_CONDITIONAL([CPU_x86_32], [test x$cpu_x86_32 = xyes])
AM_CONDITIONAL([CPU_arm], [test x$cpu_arm = xyes])
AM_CONDITIONAL([CPU_c64x], [test x$cpu_c64x = xyes])
AM_CONDITIONAL([CPU_wasm], [test x$cpu_wasm = xyes])

# Test whenever ld supports -version-script
LT_PATH_LD
AC_PROG_LD_GNU
AC_MSG_CHECKING([how to control symbol export])

THDEC_VERSION_ARG=""
THENC_VERSION_ARG=""
TH_VERSION_ARG=""
if test "x$lt_cv_prog_gnu_ld" = "xyes"; then
 case "$host_os" in
  *mingw*)
    THEORA_LDFLAGS="$THEORA_LDFLAGS -no-undefined"
    THDEC_VERSION_ARG="-export-symbols \$(top_srcdir)/win32/xmingw32/libtheoradec-all.def"
    THENC_VERSION_ARG="-export-symbols \$(top_srcdir)/win32/xmingw32/libtheoraenc-all.def"
    THENC_VERSION_ARG="$THENC_VERSION_ARG -ltheoradec"
    THC_VERSION_ARG="-export-symbols \$(top_srcdir)/win32/libtheora.def"
    AC_MSG_RESULT([-export-symbols])
    ;;
  linux* | solaris* | gnu* | k*bsd*-gnu)
    THDEC_VERSION_ARG='-Wl,--version-script=$(srcdir)/Version_script-dec'
    THENC_VERSION_ARG='-Wl,--version-script=$(srcdir)/Version_script-enc'
    TH_VERSION_ARG='-Wl,--version-script=$(srcdir)/Version_script'
    AC_MSG_RESULT([--version-script])
    ;;
  *)
    # build without versioning
    AC_MSG_RESULT([no])
    ;;
 esac
else
 case "$host_os" in
  darwin*)
    THDEC_VERSION_ARG='-Wl,-exported_symbols_list,$(srcdir)/theoradec.exp'
    THENC_VERSION_ARG='-Wl,-exported_symbols_list,$(srcdir)/theoraenc.exp'
    TH_VERSION_ARG='-Wl,-exported_symbols_list,$(srcdir)/theora.exp'
    AC_MSG_RESULT([-exported_symbols_list])
    ;;
  *)
    # build without versioning
    AC_MSG_RESULT([no])
    ;;
 esac
fi
 
THEORADEC_LDFLAGS="$THEORA_LDFLAGS $THDEC_VERSION_ARG"
THEORAENC_LDFLAGS="$THEORA_LDFLAGS $THENC_VERSION_ARG"
THEORA_LDFLAGS="$THEORA_LDFLAGS $TH_VERSION_ARG"
AC_SUBST(THEORADEC_LDFLAGS)
AC_SUBST(THEORAENC_LDFLAGS)
AC_SUBST(THEORA_LDFLAGS)

dnl --------------------------------------------------
dnl Checks for support libraries and headers
dnl --------------------------------------------------

dnl check for Ogg
HAVE_OGG=no

dnl first check through pkg-config since it's more flexible

dnl check for pkg-config itself so we don't try the m4 macro without pkg-config
AC_CHECK_PROG(HAVE_PKG_CONFIG, pkg-config, yes)

THEORA_LIBOGG_REQ_VERSION=1.3.4
AC_SUBST(THEORA_LIBOGG_REQ_VERSION)
if test "x$HAVE_PKG_CONFIG" = "xyes"
then
  PKG_CHECK_MODULES(OGG, ogg >= $THEORA_LIBOGG_REQ_VERSION, HAVE_OGG=yes, HAVE_OGG=no)
fi
if test "x$HAVE_OGG" = "xno"
then
  dnl fall back to the old school test
  XIPH_PATH_OGG(, AC_MSG_ERROR([
    libogg is required to build this package!
    please see https://www.xiph.org/ for how to
    obtain a copy.
  ]))
  cflags_save=$CFLAGS
  libs_save=$LIBS
  CFLAGS="$CFLAGS $OGG_CFLAGS"
  LIBS="$LIBS $OGG_LIBS"
  AC_CHECK_FUNC(oggpackB_read, , [
    AC_MSG_ERROR([newer libogg version ($THEORA_LIBOGG_REQ_VERSION or later) required])
  ])
  CFLAGS=$cflags_save
  LIBS=$libs_save
fi


dnl check for Vorbis
HAVE_VORBIS=no

dnl first check through pkg-config since it's more flexible

if test "x$HAVE_PKG_CONFIG" = "xyes"
then
  PKG_CHECK_MODULES(VORBIS, vorbis >= 1.0.1, HAVE_VORBIS=yes, HAVE_VORBIS=no)
  dnl also set VORBISENC_LIBS since an examples needs it
  dnl the old .m4 sets this to a value to use on top of VORBIS_LIBS,
  dnl so we do the same here.
  VORBISENC_LIBS="-lvorbisenc"
  AC_SUBST(VORBISENC_LIBS)
fi
if test "x$HAVE_VORBIS" = "xno"
then
  dnl fall back to the old school test
  XIPH_PATH_VORBIS(HAVE_VORBIS=yes, HAVE_VORBIS=no)
fi

dnl check for SDL
HAVE_SDL=no
if test "x$HAVE_PKG_CONFIG" = "xyes"
then
  PKG_CHECK_MODULES(SDL, sdl, HAVE_SDL=yes, [
    HAVE_SDL=no
    AC_MSG_WARN([*** Unable to find SDL -- Not compiling example players ***])])
fi

dnl check for OSS
HAVE_OSS=no
AC_CHECK_HEADERS([sys/soundcard.h soundcard.h machine/soundcard.h],[
  HAVE_OSS=yes
  break
])
if test x$HAVE_OSS != xyes; then
  AC_MSG_WARN([OSS audio support not found -- not compiling player_example])
fi

dnl OpenBSD needs -lossaudio to use the oss interface
OSS_LIBS=
case "$host_os" in
  openbsd*)
    OSS_LIBS='-lossaudio'
    ;;
esac
AC_SUBST(OSS_LIBS)

dnl check for libpng
HAVE_PNG=no
if test "x$HAVE_PKG_CONFIG" = "xyes"
then
  PKG_CHECK_MODULES(PNG, libpng, HAVE_PNG=yes, HAVE_PNG=no)
fi
AC_SUBST(PNG_CFLAGS)
AC_SUBST(PNG_LIBS)

dnl check for libtiff
HAVE_TIFF=no
TIFF_CFLAGS=''
TIFF_LIBS=''
AC_CHECK_LIB([tiff], [TIFFReadRGBAImage], [
  TIFF_LIBS='-ltiff'
  AC_CHECK_HEADER([tiffio.h], [
    HAVE_TIFF=yes
  ])
])
AC_SUBST(HAVE_TIFF)
AC_SUBST(TIFF_CFLAGS)
AC_SUBST(TIFF_LIBS)

dnl check for libcairo
HAVE_CAIRO=no
AC_ARG_ENABLE(telemetry,
    AS_HELP_STRING([--enable-telemetry], [Enable debugging output controls]),
    [ ac_enable_telemetry=$enableval ], [ ac_enable_telemetry=no] )

if test "x${ac_enable_telemetry}" = xyes; then
   if test "x$HAVE_PKG_CONFIG" = "xyes"
   then
     PKG_CHECK_MODULES(CAIRO, cairo, HAVE_CAIRO=yes, HAVE_CAIRO=no)
     AC_DEFINE([HAVE_CAIRO], [],  [libcairo is available for visual debugging output])  
   fi
   if test x$HAVE_CAIRO != xyes; then
     AC_MSG_WARN([libcairo not found -- not compiling telemetry output support ])
   fi
   AC_SUBST(CAIRO_CFLAGS)
   AC_SUBST(CAIRO_LIBS)
fi

dnl add safety rails for fuzzing
AC_ARG_ENABLE(mem-constraint,
    AS_HELP_STRING([--enable-mem-constraint], [Abort if size exceeds 16384x16384 (for fuzzing only)]),
    [ ac_enable_mem_constraint=$enableval ], [ ac_enable_mem_constraint=no] )

if test "x${ac_enable_mem_constraint}" = xyes; then
   AC_DEFINE([HAVE_MEMORY_CONSTRAINT], [],  [Abort if size exceeds 16384x16384 (for fuzzing only)])
fi

dnl --------------------------------------------------
dnl Overall build configuration options
dnl --------------------------------------------------

dnl Configuration option for building of encoding support.

ac_enable_encode=yes
AC_ARG_ENABLE(encode,
     AS_HELP_STRING([--disable-encode], [Disable encoding support]),
     [ ac_enable_encode=$enableval ], [ ac_enable_encode=yes] )

if test "x${ac_enable_encode}" != xyes ; then
    AC_DEFINE([THEORA_DISABLE_ENCODE], [],
  [Define to exclude encode support from the build])
else
    if test x$HAVE_VORBIS = xyes; then
      BUILDABLE_EXAMPLES="$BUILDABLE_EXAMPLES encoder_example\$(EXEEXT)"
    else
      AC_MSG_NOTICE([Vorbis missing, cannot build example encoder])
    fi
fi
AM_CONDITIONAL(THEORA_DISABLE_ENCODE, [test "x${ac_enable_encode}" != xyes])

dnl Configuration option for examples

ac_enable_examples=no
AC_ARG_ENABLE(examples,
     AS_HELP_STRING([--enable-examples], [Enable examples]),
     [ ac_enable_examples=$enableval ], [ ac_enable_examples=no] )
AM_CONDITIONAL(THEORA_ENABLE_EXAMPLES, [test "x${ac_enable_examples}" != xno])

dnl --------------------------------------------------
dnl Check for headers
dnl --------------------------------------------------

dnl none here

dnl --------------------------------------------------
dnl Check for typedefs, structures, etc
dnl --------------------------------------------------

dnl none

dnl --------------------------------------------------
dnl Check for library functions
dnl --------------------------------------------------

# The dump_video example requires either clock_gettime or ftime.
# clock_gettime is used only if time.h defines CLOCK_REALTIME and the
# function is available in the standard library; on platforms such as
# glibc < 2.17 where -lrt or another library would be required, ftime
# will be used.
AC_MSG_CHECKING([for clock_gettime])
AC_LINK_IFELSE([
  AC_LANG_PROGRAM([[#include <time.h>]], [[
    struct timespec ts;
    return clock_gettime(CLOCK_REALTIME, &ts);
  ]])
], [
  AC_MSG_RESULT([yes])
  AC_DEFINE([OP_HAVE_CLOCK_GETTIME], [1],
    [Enable use of clock_gettime function])
], [
  AC_MSG_RESULT([no])
  AC_SEARCH_LIBS(ftime, [compat], , )
])

dnl substitute the included getopt if the system doesn't support long options
AC_CHECK_FUNC(getopt_long,
              [GETOPT_OBJS=''],
              [GETOPT_OBJS='getopt.$(OBJEXT) getopt1.$(OBJEXT)'])
AC_SUBST(GETOPT_OBJS)

if test x$HAVE_SDL = xyes -a x$HAVE_OSS = xyes -a x$HAVE_VORBIS = xyes; then
  BUILDABLE_EXAMPLES="$BUILDABLE_EXAMPLES player_example\$(EXEEXT)"
fi
if test x$HAVE_PNG = xyes; then
  BUILDABLE_EXAMPLES="$BUILDABLE_EXAMPLES png2theora\$(EXEEXT)"
fi
if test x$HAVE_TIFF = xyes; then
  BUILDABLE_EXAMPLES="$BUILDABLE_EXAMPLES tiff2theora\$(EXEEXT)"
fi
AC_SUBST(BUILDABLE_EXAMPLES)

dnl --------------------------------------------------
dnl Do substitutions
dnl --------------------------------------------------

AC_SUBST(DEBUG)
AC_SUBST(PROFILE)

AC_CONFIG_FILES([
  Makefile 
  lib/Makefile
  include/Makefile include/theora/Makefile
  examples/Makefile
  doc/Makefile doc/Doxyfile doc/spec/Makefile
  tests/Makefile
  m4/Makefile
  libtheora.spec
  theora.pc
  theora-uninstalled.pc
  theoradec.pc
  theoradec-uninstalled.pc
  theoraenc.pc
  theoraenc-uninstalled.pc
])
AC_CONFIG_HEADERS([config.h])

AC_OUTPUT

AS_AC_EXPAND(LIBDIR, ${libdir})
AS_AC_EXPAND(INCLUDEDIR, ${includedir})
AS_AC_EXPAND(BINDIR, ${bindir})
AS_AC_EXPAND(DOCDIR, ${docdir})

if test $HAVE_DOXYGEN = "false"; then
  doc_build="no"
else
  doc_build="yes"
fi
if test $BUILD_SPEC = "false"; then
  spec_build="no"
else
  spec_build="yes"
fi

AC_MSG_RESULT([
------------------------------------------------------------------------
  $PACKAGE $VERSION:  Automatic configuration OK.

  General configuration:

    Encoding support: ........... ${ac_enable_encode}
    Assembly optimization: ...... ${cpu_optimization}
    Debugging telemetry: ........ ${ac_enable_telemetry}
    Abort on huge files: ........ ${ac_enable_mem_constraint}
    Build example code: ......... ${ac_enable_examples}
    API Documentation: .......... ${doc_build}
    Format Documentation: ....... ${spec_build}

  Installation paths:

    libtheora: ................... ${LIBDIR}
    C header files: .............. ${INCLUDEDIR}/theora
    Documentation: ............... ${DOCDIR}

  Building:

    Type 'make' to compile $PACKAGE.

    Type 'make install' to install $PACKAGE.

    ${TESTS_INFO}

  Example programs will be built but not installed.
------------------------------------------------------------------------
])

