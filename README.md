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
1. convert LCIO file, which is from the data collector, to more easier format(ROOT is most desirable for us)
2. correct channel assignment
3. make output file that contains the reconstruction
4. database of MPPC gain, MPPC pedestal, fiber gain
5. make a program to solve ambiguity
6. make a program to determine the track

As for reconstruction, essential information is the angle because the number of reconstruction points is not constant(five or six).

So the number of most essential parameters from the EASIROC side is four(two directions times two counters)
These parameters are eventually compared with parameters from the HBU side.
