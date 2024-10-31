This repository contains the code for a model of signal processing in the macaque retina.

The code runs in Matlab. Ensure all .m and .mat files are on the Matlab path.

Run runCol. This should plot a map of cone locations.

Now edit runCol. The main switch statement controls the analysis. You just ran case array.x.y: try some other cases.

All cases set metadata used by stream, which executes the analyses. See Guide to stream and streamTute.m for short tutorials on stream.

The file setCol contains code for setting model parameters and for running the ratio-of-Gaussians (ROG) model. After unzipping the data file Colour.zip, run setCol to fit the ROG to empirical data.
