# polyfit

Latest version: polyfit-0.13, Sep 17, 2015

INTRODUCTION
------------

Polyfit is a polynomial chain fitting tool developed to pre-process eclipsing
binary (EB) light curves and get them ready for the artificial intelligence
based engine EBAI. The code is released under GNU General Public License,
which allows you to freely use, share and modify the code. For the details
on the computational logic of the code please refer to Prsa et al. (2008),
ApJ 687, 542.

COMPILING
---------

To compile polyfit, you must be running a GNU-compatible build system. This
means that polyfit will build readily on any linux platform; the code should
be readily portable to other platforms, but we have not tested this yet. Any
success reports would be greatly appreciated!

Polyfit depends on the GNU Scientific Library (GSL), although this dependency
will be removed in the upcoming version. Any C compiler, such as gcc, should
compile the code without any problems.

To build polyfit, modify Makefile to your liking and issue `make`. That will
create the polyfit executable.

To install polyfit, copy the executable to /usr/local/bin or any alternative
directory in the bin path.

USAGE
-----

Polyfit is run on light curves. The input file can either contain 1 column
(equidistant fluxes), 2 columns (phase and flux), or three columns (phase,
flux and standard deviation). You can additionally pass the following switches
to polyfit:

  -o order          ..  fitting polynomial order (default: 2)
  -i iters          ..  number of iterations (default: 10000)
  -s step           ..  step for random knot displacement (default: 0.01)
  -k k1 k2 ... kN   ..  explicit list of knots
  -n vertices       ..  number of equidistant vertices in the computed fit
  --find-knots      ..  attempt to find knots automatically
  --find-step       ..  attempt to find step automatically
  --chain-length    ..  minimum chain length for automatic knot search
  --apply-pshift    ..  shift phase so that the polyfit minimum is at phase 0
  --residuals       ..  output residuals instead of the fit itself
  --csv-input       ..  comma-separated value file input
  --csv-output      ..  comma-separated value output
  --summary-output  ..  output a single summary line from polyfit

The best way to get started is to use the demo input file 'lc.dat'. Invoke
polyfit with:

    polyfit --find-knots -s 0.005 -n 500 lc.dat > lc.out
    polyfit --find-knots --find-step --csv-input -c 1 2 --csv-output \
            lc.csv > pf.csv

This will compute a theoretical curve that consists of 4 segments, each fit
with a quadratic function. You may want to display it along with the data.
In gnuplot, for example, you might do:

    plot 'lc.dat', 'lc.out' with lines
    set datafile separator ","
    plot 'lc.csv' u 2:3 w p, 'pf.csv' w l

WHAT IS NEW IN VERSION 0.13
---------------------------

Comma-separated-value files (csv) can now be passed to polyfit by using
the --csv-input and --csv-output switches.

WHAT IS NEW IN VERSION 0.12
---------------------------

When eclipses are very narrow, the test knots could overlap and cause the
segfault if --find-step is *not* used. This is now fixed by always sorting
the array of test knots.

The 'dist' and 'install' rules added to Makefile.
