/* common.h -- Process this file with configure to produce common.h
   Copyright (C) 2000, 2006 Gary V. Vaughan

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA.
*/
/** @start 1 */
#ifndef PSI_COMMON_H
#define PSI_COMMON_H 1

#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#if STDC_HEADERS
#  include <stdlib.h>
#  include <string.h>
#elif HAVE_STRINGS_H
#  include <strings.h>
#endif /*STDC_HEADERS*/

#if HAVE_UNISTD_H
#  include <unistd.h>
#endif

#if HAVE_ASSERT_H
#  include <assert.h>
#endif

#if HAVE_MATH_H
#  include <math.h>
#endif

#if HAVE_FLOAT_H
#  include <float.h>
#endif

#if HAVE_STDARG_H
#  include <stdarg.h>
#endif

#if HAVE_LIBGEN_H
#  include <libgen.h>
#endif

#if HAVE_TIME_H
#  include <time.h>
#endif

#if HAVE_GETOPT_H
#  include <getopt.h>
#endif

/** @end 1 */
#if HAVE_SYS_WAIT_H
#  include <sys/wait.h>
#endif
#ifndef WIFEXITED
#  define WIFEXITED(stat)       (((stat) & 0xff) == 0)
#endif
#ifndef WEXITSTATUS
#  define WEXITSTATUS(stat)     ((unsigned)(stat) >> 8)
#endif
#ifndef WIFSTOPPED
#  define WIFSTOPPED(stat)      (((stat) & 0xff) == 0x7f)
#endif
#ifndef WSTOPSIG
#  define WSTOPSIG(stat)        WEXITSTATUS(stat)
#endif
#ifndef WIFSIGNALED
#  define WIFSIGNALED(stat)     (!WIFEXITED(stat) && !WIFSTOPPED(stat))
#endif
#ifndef WTERMSIG
#  define WTERMSIG(stat)        ((stat) & 0x7f)
#endif

/** @start 1 */
#if HAVE_ERRNO_H
#  include <errno.h>
#endif /*HAVE_ERRNO_H*/
#ifndef errno
/* Some systems #define this! */
extern int errno;
#endif

/** @end 1 */
/** @start 4 */
#ifdef __cplusplus
#  define BEGIN_C_DECLS         extern "C" {
#  define END_C_DECLS           }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

/** @end 4 */
#ifdef __GNUC__
#  ifndef const
#    define const       __const
#  endif
#  ifndef signed
#    define signed      __signed
#  endif
#  ifndef volatile
#    define volatile    __volatile
#  endif
#else
#  ifdef __STDC__
#    undef  signed
#    define signed
#    undef  volatile
#    define volatile
#  endif
#endif

#ifdef __STDC__
#define STR(x)          #x
#define CONC(x, y)      x##y
#else
#define STR(x)          "x"
#define CONC(x, y)      x/**/y
#endif

/** @start 2 */
#define EXIT_OTHERS   2
#ifndef EXIT_SUCCESS
#  define EXIT_SUCCESS  0
#  define EXIT_FAILURE  1
#endif
/** @end 2 */
#if !HAVE_BZERO && HAVE_MEMSET
# define bzero(buf, bytes)      ((void) memset (buf, 0, bytes))
#endif

#if !HAVE_STRCHR
#  define strchr index
#endif

#if !HAVE_STRRCHR
#  define strrchr rindex
#endif

#if !HAVE_ISFINITE
#  define isfinite finite
#endif

/** @start 4 */
#define XCALLOC(type, num)                                  \
        ((type *) xcalloc ((num), sizeof(type)))
#define XMALLOC(type, num)                                  \
        ((type *) xmalloc ((num) * sizeof(type)))
#define XREALLOC(type, p, num)                              \
        ((type *) xrealloc ((p), (num) * sizeof(type)))
#define XFREE(stale)                            do {        \
        if (stale) { free (stale);  stale = 0; }            \
                                                } while (0)

BEGIN_C_DECLS

extern void *xcalloc    (size_t num, size_t size);
extern void *xmalloc    (size_t num);
extern void *xrealloc   (void *p, size_t num);
extern char *xstrdup    (const char *string);
extern char *xstrerror  (int errnum);

/** @end 4 */
#if !HAVE_BASENAME
extern char *basename   (const char *path);
#endif

#if !HAVE_STRCSPN
extern size_t strcspn   (const char *string, const char *accept);
#endif

#if !HAVE_STRERROR
extern char *strerror   (int err);
#endif

#if !HAVE_STRSIGNAL
extern char *strsignal  (int signo);
#endif

#if !HAVE_STRSPN
extern size_t strspn    (const char *string, const char *reject);
#endif

#if !HAVE_WAITPID
extern pid_t waitpid    (pid_t pid, int *pstatus, int options);
#endif
/** @start 4 */

#define GMEL_N_MIN 0.000001
#define GMEL_N_MAX 0.999999

enum { 
  GEN_INIT, 
  GEN_NEXT 
};

enum {
  PSI_SAMPLE_A,
  PSI_SAMPLE_C,
  PSI_SAMPLE_G,
  PSI_SAMPLE_T,
  PSI_SAMPLE_S,
  PSI_SAMPLE_P,
  PSI_NUM_SAMPLE
};

enum {
  PSI_SAMPLE_K = 6,
  PSI_SAMPLE_W
};



/*!\brief 

  For a grid point, the following five things are stored:
  S: \f$-2E_s(i)\f$
  P: \f$-2E_p(i)\f$
  A: number of nucleotide A's
  C: number of nucleotide C's
  G: number of nucleotide G's
 */
typedef struct gibbsPart{
   double S, P;
   int A, C, G;
} GibbsPart;

/*!\brief grid point 

   This is a C-structure representation of a grid point.
 */
typedef struct Gridpoint {
   int s;
   int p;
   int a;
   int c;
   int g;
} gridpoint;

/*!\brief parameter

   This is a C-structure representation of a parameter \f$theta\f$. It is a
   vector having a list of parameters.
 */
typedef struct Parameter {
   double s;
   double p;
   double a;
   double c;
   double g;
   double t;
   int w;
   gridpoint gp;
   double max_energy;
   double gibbs_sum;
} parameter;

END_C_DECLS
/** @end 4 */
/** @start 1 */
#endif /* !PSI_COMMON_H */
/** @end 1 */
