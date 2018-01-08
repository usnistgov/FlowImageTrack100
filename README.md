# FlowImageTrack100
Counting of 100% of particles passing through a flow imaging system

8 JAN 2018  
Dean Ripple  

General Notes  
This program is an implementation of the algorithm discussed in Dean C. Ripple & Paul C. DeRose, “Primary Determination of Particle Concentration 
with Light Obscuration and Dynamic Imaging Particle Counters“ (2018) J. Res. Nat. Inst. Stds. Tech.

The actual program used for the paper was compiled FORTRAN code.  That code used libraries that cannot be distributed, so a Python equivalent is given here.  Because the routine requires many least-squares fits with a small number of data points per fit, Python is not a particularly fast language for this application.  Translation into a compiled language such as C or FORTRAN will give large increases in processing speed.

Files:  
FlowImageTrack100.py:  the Python code    
Track100_Sample_In.csv:  a relatively short input file created by truncating one of the files from Ripple & DeRose.  
Track100_Sample_Out.csv:  sample output file using as program parameters 0.301 second time window, 5 pixel x variation, 500 pixel y variation, 0.3 diameter threshold.  

Notes for FlowImageTrack100.py  
    VERSION: 1.1  
    DATE:  08 JAN 2018  
    PURPOSE:  Sorting of flow imaging exported data files for analyzing apparent  
    particle densities by sedimentation properties and for performing 100%  
    particle counts  
    INPUT & OUTPUT FILES:  Both are .csv format  
    PYTHON VERSION:  Tested on Python 3.6  
    INPUT FILE: specified number of header rows; succeeding rows have:  
      particle image id, area (um^2), x & y corner location (pixels), diameter  
      (um), elapsed time (seconds), image height & width (pixels) 
    OUTPUT FILE:   
      Input file name echoed  
      Allowable change in x & y between images (pixels), allowable fractional  
      change in dia., time to look ahead for next matching image (seconds)  
      Rows reporting: Part. #, ID 1st image, # images in track, ave. time  
      (seconds), ave. x (pixels), ave. diameter (um), slope of linear fit  
      (um per second), rms dev. from fit  
    Optional output section containing:   
      ID, diameter (um), area-based diameter (ABD) derived from the reported particle area (um), area (um^2), center x (pixels), 
      time (seconds), y (pixels), x_corner (pixels), y_corner  
      (pixels), track number  

      For 100% counting, the number of particles found = the total number of  
      tracks found that have at least the specified n_min number of tracks.  
      Since the track number will increment even if there is a single particle  
      in the track, the total particle number is not equal to the highest track  
      number.  The number of valid tracks is reported at the end of the first output section.  
