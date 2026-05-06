KiT (Kinetochore Tracking) software
===================================

KiT is a MATLAB program for tracking kinetochores in microscopy images and movies
in 2D or 3D, and in multiple channels.

This software is licensed for use under the GNU GPL v3 license. See LICENSE file
for details.

Installation
------------

No installation is necessary, simply copy the source code to the directory of your
choice. Enter that directory within MATLAB (use directory toolbar at top, or 'cd'
command) and then enter commands.

Usage
-----

Detailed instructions for single-timepoint tracking for Mad2 pseudo-timeline creation
are provided in the accompanying document 'KiT how to v2.4.1.pdf'. For live-cell
tracking, please see 'KiT2manual_2.1.11.pdf' for a legacy document on how to do this.

A tracking experiment is organized in KiT as a 'jobset' stored in a mat-file. To
create a jobset:

	kitRun

This will display a GUI that will enable you to load in or create new jobsets. If you
choose to create a new jobset, you will have the options of 'Tracking' for a movie,
'Single timepoint' for fixed cell images, or 'Chromatic shift' (should be use on
single-timepoint images).

This will display a GUI to allow modification of tracking parameters. Click 'Select
directory' to choose a directory to search for image files. Select the images to
analyse by highlighting them in the list and click 'Add' or 'Add and Crop'. If 'Add
and Crop' was selected, for each image draw one or more ROI(s) around the cell(s) to
be tracked. Double click in each ROI after drawing to save it. After ROI selection,
configure tracking parameters. Typically the default tracking parameters will be fine,
but may need to be changed on a case-by-case basis. Choose a name for the jobset.
Validate your data, then click 'Save jobset' to return to the main KiT window. You
can click 'Run job' to run tracking immediately, or 'Close' to come back to it later.

After tracking is complete you will find a file named something like
'kittracking001_expname_moviefile.mat', where expname is replaced with the name of
jobset mat-file and moviefile is replaced with the name of the movie this result
is associated with.

To load a previously created jobset mat-file use:

   jobset = kitLoadJobset;

To save a jobset to its mat-file after making changes:

   kitSaveJobset(jobset);

To load a job's results from the 'kittracking*.mat' file, e.g. movie #2 of a
jobset:

   job = kitLoadJob(jobset,2);


Many functions within KiT have optional parameters. The 'help' command is your
friend here...
 

Bugs
----

Send bug or crash reports and feature requests to catriona.c.conway(at)warwick.ac.uk.

Include:
  - Brief description of steps to reproduce error.
  - Copy of the error text.
  - Jobset mat-file (if available).

Movies will be useful but should be sent via e.g. Dropbox or Files@Warwick.


Credits
-------

- KiT evolved from both Maki (Jaqaman, et al., JCB 188:665-679 (2010)) and
  u-Track (Jaqaman, et al., Nat. Meth. 5:695-702 (2008)). Large portions of the
  core tracking code is lifted directly from Maki and u-Track.

- E. Harry modified the u-Track mixture-model fitting to operate on 3D data,
  plus various other tweaks throughout.

- E. Vladimirou made modifications to the coordinate system fitting for
  non-congressed kinetochores, plus various other tweaks throughout.

- C. A. Smith contributed some bugfixes and additions to some of the plotting code.

- J. W. Armond organised the code into its current form as KiT. He added the
  adaptive wavelet spot finder, image-moment based coordinate system, spot intensity
  measurement, multiple channel tracking, and movie rendering. He also made a
  number of improvements toward code cleanliness and efficiency.

- C. C. Conway contributed bugfixes and additions for fixed cell pseudo-timeline
  creation and graph plotting.

- A. Daniyan contributed bugfixes.


The KiT package contains code from various other sources. See file headers to
identity sources and licenses.


MATLAB requirements
--------------------

MATLAB R2014b-2024b (Intel installation)
MATLAB Image Processing Toolbox 9.1
MATLAB Statistics Toolbox 9.1

Optional toolboxes
------------------

For parallel execution of tracking
    MATLAB Parallel Computing Toolbox 6.5
For adaptive threshold particle detection
    MATLAB Global Optimization Toolbox 3.3
For correction of photobleaching in intensity data and enhancing
adaptive particle finding algorithm
    MATLAB Curve Fitting Toolbox 3.5
