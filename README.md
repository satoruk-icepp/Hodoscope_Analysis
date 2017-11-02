# Hodoscope_Analysis
ROOT(C++) Sources for the analysis of HBU calibration with two hodoscopes

# Directory Structure
Parent Directory/

+src/--this file

+data/--raw EASIROC file

+rootfile/--output of these programs

At this point, you can analyze the data from two EASIROC modules.

To extract the raw data, run "root "tEASIROC(<runnum>)"".
Then you can find two files are created in rootfile/.

Tasks
- make the geometry the right position(modify channel assign)
- make output file that contains the reconstruction
- database of MPPC gain, MPPC pedestal, fiber gain
- convert LCIO file to more easier format(for us)
- make program to solve ambiguity
- make program to determine the track
