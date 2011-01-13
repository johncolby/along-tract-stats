*Along*-tract statistics
========================
1. Overview
2. Install
3. Support
4. Contact
5. License

Overview
--------
This is a set of tools for conducting an along-tract analysis of white matter fiber tracts derived from diffusion MRI data. It allows you to analyze a scalar metric (e.g. fractional ansiotropy; FA) parameterized along a tract, rather than the more typical method of collapsing the variability in these measures into single tract-averaged mean estimates.

This package contains 3 things:

1. **Tools** - A modular set of MATLAB functions to perform individual tasks on the tract groups (load, plot, save, interpolate, etc.).
2. **Experiment wrappers** - Several examples of how you can link these tools together to perform a full analysis. These will be a good starting point to adapt to your own needs.
3. **Example data** - So you can get a feel for running these tools and verify that they are behaving the same way as in the example documentation.

Install
-------
### Requirements

- MATLAB
- MATLAB Spline toolbox (or Curve Fitting Toolbox as of r2010b)
- FSL (mainly for `read_avw`)
 
### Install
1. Click the Downloads link towards the top-right on Github. Download and extract either the .zip or .tar.gz versions.
2. Add the `along-tract-stats` directory to your matlab path.
        addpath('/path/to/along-tract-stats')
3. Add the FSL `matlab` directory to your matlab path.
        addpath('/path/to/fsl/etc/matlab')

Support and Usage
-----------------
- [Github Wiki](http://github.com/johncolby/along-tract-stats/wiki/wiki) - The main source of online documentation. Information on basic usage and tutorials.
- [Colbyimaging Wiki](http://www.colbyimaging.com/wiki/neuroimaging/along-tract-stats) - Video tutorials (can't be embedded on Github).
- [Google groups](http://groups.google.com/group/along-tract-stats) - Questions, comments, feedback, bugs, etc.
- Offline documentataion - There are several sources of offline documentation inluded with this package:
    - MATLAB help - Type `help <command_name>` at the MATLAB command prompt (e.g. `help trk_plot`) to get specific info on usage, inputs, and outputs for each function.
    - MATLAB demo - Open the index.html file in the `html` directory to see the demo output that is published by `trk_demo.m`
    - Inline code comments - The code is written in plain text with comments throughout.
    
Contact
-------
Support requests should be posted on the Google Groups discussion forum. For personal contact, use:

John Colby

Email: <a href="&#x6D;&#x61;&#105;&#x6C;&#116;&#x6F;:&#106;&#111;&#104;&#x6E;&#x63;&#x6F;&#108;&#98;&#121;&#64;&#x75;c&#x6C;&#97;&#x2E;&#101;&#x64;&#117;">&#106;&#111;&#104;&#x6E;&#x63;&#x6F;&#108;&#98;&#121;&#64;&#x75;c&#x6C;&#97;&#x2E;&#101;&#x64;&#117;</a>  
UCLA Developmental Cognitive Neuroimaging Laboratory  
Department of Neurology  
David Geffen School of Medicine at UCLA  

License
-------
Copyright 2010, John Colby

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.