*Along*-tract statistics
========================
1. Overview
2. Install
3. Support
4. License

Overview
--------
This is a set of tools for conducting an along-tract analysis of white matter fiber tracts derived from diffusion MRI data. It allows you to analyze a scalar metric (e.g. fractional anisotropy; FA) parameterized along a tract, rather than the more typical method of collapsing the variability in these measures into single tract-averaged mean estimates.

This package contains 3 things:

1. **Tools** - A modular set of MATLAB functions to perform individual tasks on the tract groups (load, plot, save, interpolate, etc.).
2. **Experiment wrappers** - Several examples of how you can link these tools together to perform a full analysis. These will be a good starting point to adapt to your own needs.
3. **Example data** - So you can get a feel for running these tools and verify that they are behaving the same way as in the example documentation.

Install
-------
### Requirements

- [MATLAB](http://www.mathworks.com/products/matlab)
- [MATLAB Curve Fitting Toolbox](http://www.mathworks.com/products/curvefitting) (or Spline Toolbox prior to r2010b)
- [FSL](http://www.fmrib.ox.ac.uk/fsl) (mainly for `read_avw`)

### Optional
- [MATLAB Statistics Toolbox](http://www.mathworks.com/products/statistics): Nicely handle Excel-type data tables
- [MATLAB Parallel Computing Toolbox](http://www.mathworks.com/products/parallel-computing): Speed up streamline processing
- [TrackVis](http://www.trackvis.org): 3D visualization
- [R](http://www.r-project.org): Statistical inference and graphics
 
### Install
1. Click the [Downloads](http://github.com/johncolby/along-tract-stats/archives/master) link towards the top-right on Github. Download and extract either the .zip or .tar.gz versions.
2. Add the `along-tract-stats` directory to your matlab path. (Adjust these paths according to your setup)
        addpath('/path/to/along-tract-stats')
3. Add the FSL `matlab` directory to your matlab path.
        addpath('/usr/local/fsl/etc/matlab')
4. Save these changes with `savepath`.
5. Set the `$FSLDIR` environment variable in MATLAB.
        setenv('FSLDIR', '/usr/local/fsl')
6. To save this setting for future sessions, consider putting it in your [`startup.m`](http://www.mathworks.com/help/techdoc/ref/startup.html) file.

Support and Usage
-----------------
- [Github Wiki](http://github.com/johncolby/along-tract-stats/wiki) - The main source of online documentation. Information on basic usage and tutorials.
- [Colbyimaging Wiki](http://www.colbyimaging.com/wiki/neuroimaging/along-tract-stats) - Video tutorials (can't be embedded on Github).
- [Google groups](http://groups.google.com/group/along-tract-stats) - Questions, comments, feedback, bugs, etc.
- Offline documentation - There are several sources of offline documentation included with this package:
    - MATLAB help - Type `help <command_name>` at the MATLAB command prompt (e.g. `help trk_plot`) to get specific info on usage, inputs, and outputs for each function.
    - MATLAB demo - Open the index.html file in the `html` directory to see the demo output that is published by `trk_demo.m`
    - Inline code comments - The code is written in plain text with comments throughout. 

License
-------
Copyright 2010, John Colby

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.