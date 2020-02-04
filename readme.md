## Skycover

This code is for simulating the probes that collect light for the
Aquisition Guide Star and Wave Front Sensing (AGWS) system of the
Giant Magellan Telescope (GMT). The goal of this project is to
determine the probabilities of finding sufficiently bright stars for
doing wavefront sensing.

## Dependencies

This code depends on the Geometry and Foreach libraries from the Boost package found at boost.org. These are "header only" libaries so only the include files need to be pointed to in the Makefile.

## Running

I wrote this code on OS X 10.11.4 (El Capitan). This is what I get
when I run `g++ --version`.

    Configured with: --prefix=/Applications/Xcode.app/Contents/Developer/usr --with-gxx-include-dir=/usr/
    include/c++/4.2.1                                                                                   
    Apple LLVM version 7.3.0 (clang-703.0.31)
    Target: x86_64-apple-darwin15.6.0
    Thread model: posix
    InstalledDir: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin

To build the project, run `make` in the top level directory. This will
compile an executable named <b>skycov</b>.

The program takes eight command line arguments, all mandatory. Here is
the usage message:

    usage: ./skycov <--4probe | --phasing> <--gclef | --m3 | --dgnf | --dgwf> <--track | --notrack> <--print | --noprint> <wfsmag> <gdrmag> <phmag> <nfiles>

And brief explanations of the arguments:

    arg 1: regular simulation or phasing
    arg 2: obscuration type. '--dgnf' for no obscuration
    arg 3: tracking or notracking
    arg 4: print configurations or not
    arg 5: wfsmag
    arg 6: gdrmag
    arg 7: phmag -- phasing star magnitude
    arg 8: number of files to test

The program will read in star field data from the star catalogues in
the <b>Bes</b> directory, and return the probability of finding the
given wave-front/guide-star magnitude pair in any one of the
files. The probability formula is

    (# of valid files) / (# of files tested)

The validity of a file depends on which test is being run. There are
two test types, <b>4probe</b> and <b>phasing</b>.

For the <b>4probe</b> test, a valid file contains a configuration of
stars such that three probes can reach a star of at least the given
wave-front magnitude, and one probe can reach a star of at least the
given guide-star magnitude. A valid configuration must also be one
where none of the probes are colliding with each other, or are blocked
by an obscuration.

For the <b>phasing</b> test, only three probes need to be able to
reach a star of at least the given wave-front magnitude. The
guide-star magnitude is ignored in this case. The phasing test also
takes into account probe collisions and obscuration by M3.

The second command line argument determines which obscuration, if any,
will be used. Passing in '--gclef' will cause an obscuration polygon
to be read out of the file 'gclef_obsc.txt'. Passing in '--m3' will
cause an obscuration polygon to be read out of 'm3_obsc.txt'. And
passing in '--dgnf' will cause no obscuration to be used, as there is
no obscuration to take into account in the direct gregorian narrow
field configuration.  With '--dgwf', the AGWS sees no obscuration but instead creates shadows to the science focal plane below.  The size of the shadowed area must be less than 15% of the 20 arcminute diameter DGWF field of view.

The format of the obscuration files must be one ordered pair per line,
with x and y coordinates separated by a single tab.

Ex. gclef_obsc.txt

    -306.63	1388.34
    330.15	-149
    -128.09	-338.81
    -764.88	1198.53

The order of the points given in the obscuration file should follow
clockwise or counterclockwise order around the polygon, as the
geometry algorithms in the program use this format to reason about
polygons.

The third command line argument toggles the tracking test. The
tracking configuration tests whether all four probes are able to
follow a bright enough star as the field of view rotates for 60
degrees above the stationary mirror. This means for the <b>4probe</b>
configuration, three probes must follow a star at least as bright as
the given wave-front magnitude, and one probe must follow a star at
least as bright as the given guide-star magnitude. A caveat of the
tracking configuration is that any probe may 'backtrack' to another
star of sufficient magnitude if the star it is currently following
leaves its range. This backtracking ability is simulated in the
program, and as usual probe collisions and obscuration are taken into
account.

The fourth command line argument tells the program whether or not you
want the coordinates of each probe and obscuration to be printed out
at every valid configuration. For the non-tracking test, this results
in a single frame of animation per valid configuration. For the
tracking test, this results in 60 frames of animation per valid
configuration.

When using the --print option, it's a good idea to redirect output to
some output file. Example

     ./skycov --4probe --gclef --track --print 14 15 13.5 10 > sky.out

Since the resulting probability will be printed to stderr, you get
your answer and the printed coordinates in separate locations.

## Running a full test

I've included two scripts, <b>run4probe.sh</b> and
<b>runphasing.sh</b>, that will test a whole range of magnitudes at
once. <b>run4probe.sh</b> tests all possible wfs/gdr magnitude pairs
between 13 and 19, inclusive. <b>runphasing.sh</b> tests all wfs
magnitudes between 13 and 19. The 4probe script takes the same
arguments as the skycov executable, minus the simulation type and
given magnitudes.

    ./run4probe.sh <--gclef | --m3 | --dgnf | --dgwf> <--notrack | --track> <nfiles>

The phasing script takes only the number of files to test, as the
phasing configuration must use M3 and does not do tracking.

    ./runphasing <nfiles>

## Probe geometry

The four probes are identical in shape, and are described by the
polygon files 'probe_slider_body.txt', 'probe_slider_shaft.txt', and
'probe_baffle_tube.txt'. Each file contains one of three polygons that
make up the probe's geometry. The probe described as these polygons is
rotated about the origin to describe the other three probes.

These three files are generated by the script agwsprobe.py.  This
script takes one argument which is the buffer size in mm by which to
increase the effective probe size so as to set a minimum distance
between the physical probes.  The minimum distance between probes is
twice the buffer size.  The script produces a figure probe.png showing
the probe geometry.  The default buffer size is 15mm.  To change the
geometry, run the agwsprobe.py script before running skycov.


There are two reasons the probes are described as multiple
polygons. First, the baffle tube extends and contracts based on the
radial distance of the probe's star from the origin. Second, the
polygon intersection algorithm used by the program cannot handle
concave shapes.

## Results

The output of the run4probe.sh script will look something like this

    grid	config	wfsmag	gdrmag	prob
    ----	------	------	------	----
    grid	dgnf	13	13	0.01
    grid	dgnf	13	14	0.07
    grid	dgnf	13	15	0.16
    grid	dgnf	13	16	0.19
    grid	dgnf	13	17	0.19
    grid	dgnf	13	18	0.19
    grid	dgnf	13	19	0.19
    grid	dgnf	14	13	0.17
    grid	dgnf	14	14	0.18
    grid	dgnf	14	15	0.53
    grid	dgnf	14	16	0.66
    grid	dgnf	14	17	0.68
    grid	dgnf	14	18	0.68
    grid	dgnf	14	19	0.68
    grid	dgnf	15	13	0.57
    grid	dgnf	15	14	0.69
    grid	dgnf	15	15	0.7
    grid	dgnf	15	16	0.98
    grid	dgnf	15	17	0.99
    grid	dgnf	15	18	0.99
    grid	dgnf	15	19	0.99
    grid	dgnf	16	13	0.77
    grid	dgnf	16	14	0.98
    grid	dgnf	16	15	0.99
    grid	dgnf	16	16	0.99
    grid	dgnf	16	17	1
    grid	dgnf	16	18	1
    grid	dgnf	16	19	1
    grid	dgnf	17	13	0.78
    grid	dgnf	17	14	0.99
    grid	dgnf	17	15	1
    grid	dgnf	17	16	1
    grid	dgnf	17	17	1
    grid	dgnf	17	18	1
    grid	dgnf	17	19	1
    grid	dgnf	18	13	0.78
    grid	dgnf	18	14	0.99
    grid	dgnf	18	15	1
    grid	dgnf	18	16	1
    grid	dgnf	18	17	1
    grid	dgnf	18	18	1
    grid	dgnf	18	19	1
    grid	dgnf	19	13	0.78
    grid	dgnf	19	14	0.99
    grid	dgnf	19	15	1
    grid	dgnf	19	16	1
    grid	dgnf	19	17	1
    grid	dgnf	19	18	1
    grid	dgnf	19	19	1

This is starbase table format with two header rows and columns
separated by a single tab. Included in the repo is a python script
called <b>plotgrid.py</b>, which can be used to create graphs using
output in the above format. The plotgrid script takes three command
line arguments.

    python plotgrid.py <input file> <output file> <plot title>

The input file should be in the starbase format above. The output of
run4probe.sh. The output file will be a .png image. And the plot title
would be something like "M3 Obscuration - 60 degree Tracking (25mm
buffer)."

![DGNF - No Tracking (10mm buffer)](plots/4probe_100dgnf_notrack_10mm.png)

![M3 - No Tracking (10mm buffer)](plots/4probe_100m3_notrack_10mm.png)

![M3 - Tracking (10mm buffer)](plots/4probe_100m3_track_10mm.png)

![GCLEF - Tracking (10mm buffer)](plots/4probe_100gclef_track_10mm.png)

## Validate probe positions

If instead of computing the sky coverage probability you simply want to validate a number of probe configurations, use the program <b>agwsvalid</b>.

    ./agwsvalid <--m3 | --gclef | --dgnf | --dgwf> <--plot | --bool> < inputfile

The input file consists of 8 numbers per line: four x y pairs of star positions expressed in degrees.  Lines beginning with # will be ignored.

The obscuration configurations are the same as above in <b>skycover</b>.  The --bool flag will cause the program to print 1 or 0 for each line processed.   Alternately, the flag --plot will generate output suitable for use in the visual simulation described in the next section.  Text diagnostics are printed to standard error.

An example text file <b>slist.txt</b> is provided:

    ./agwsvalid --m3 --bool <slist.txt

In python you can call agwsvalid as follows:

   import agwsvalid
   mode = 'dgnf'
   dgnf = agwsvalid.validator(mode)
   result = dgnf.check(0 , .12, -.12, 0, 0, -.12, .12, 0)

where mode should be 'dgnf', 'dgwf', or 'm3', and result is a boolean.
See also agwsvalid-example.py

   

## Visual Simulation

The '--print' option enables viewing the valid configurations found by
the program. When printing is turned on, the program will output the
coordinates of all the polygons in the system whenever a valid
configuration is found. This means for the non-tracking setting, a
single set of coordinates will be printed out for each polygon in the
system. For the tracking setting, a set of coordinates will be printed
out for every polygon in the system, for every one of 60 positions
along the 60 degree tracking path.

With the print option, the program will also attempt to print the
valid starfields to files in a directory named 'starfiles'. So if you
would like to use the print option to see the visual simulation, you
should create a folder in the toplevel called 'starfiles'.

The coordinates of a polygon are printed out in ordered pairs in a
single line. The m3 obscuration is printed out like this

    (829.100000, 918.240000) (80.740000, 66.650000) (10.100000, 104.200000) (297.630000, 1200.830000) 

I've included matlab scripts that will plot the configurations output
by the simulation program. But since the coordinates are output as
text it's possible to write a parser/plotter in any language.

This is a frame displaying one valid configuration output by the
command

    ./skycov --4probe --m3 --notrack --print 15 15 13.5 10 > sky.out

![Valid config. example](/plots/valid_config_ex.png)

To run the visual simulations, open matlab and cd to the matlab
directory. To run a non-tracking simulation, run
`iterprobes('../starfiles', '../sky.out')`. The arguments to
iterprobes are the directory containing the starfields and the file
containing the polygon coordinates printed out by the simulation. To
run a tracking simulation, run `iterprobes_tracking('../starfiles',
'../sky.out')`. The arguments are the same.

Beware there are a few limitations of the matlab scripts that I
haven't found a good way to abstract away. If you would like to
simulate a dgnf test, you must go into the file, either iterprobes.m
or iterprobes_tracking.m, and change the 'dims' variable at the top of
the script from [4 4 4 4] to [4 4 4]. This is because the dgnf doesn't
have an obscuration polygon to read out of the file. The 'dims'
variable gives the number of point on the polygons described in the
file.

Another thing you must do if you would like to run a dgnf visual
simulation is go into the file, iterprobes.m or iterprobes_tracking.m,
and change the 'npolygons' variable from 13 to 12. This tells the
script to display 12 polygons (3 for each probe) per frame, instead of
13 (3 for each probe plus the obscuration).

The blue dots in the picture represent stars. These stars are read
from conventionally named text files in the 'starfiles'
directory. Each run of the simulation overwrites these files to
contain the starfields from valid files found during that
simulation. So the simulation might look off if you are viewing
polygons from a few simulations ago.

## Making a Movie of the Visual Simulation

The iterprobes_tracking matlab function records each frame it displays
in a matrix that it returns to the caller. If you are just interested
in watching the simulation once, there is no need to capture this
value. But returning the captured frames means it is possible to
record a visual simulation as a movie file.

To capture the recorded frames from an iterprobes run, assign a return
variable like this

    M = iterprobes_tracking('../starfiles', '../sky.out')

Then you can call the included iterprobes_movie function like so,
passing in a filename where the movie will be stored.

    iterprobes_movie(M, 'iterprobes.avi')

Be careful about how many configurations you record in the movie. File
size will increase rapidly. A movie showing four valid tracking
configurations can run to almost 20M.