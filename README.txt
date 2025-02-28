Authors: dalton durant --- ddurant@utexas.edu, thedaltondurant@gmail.com

Run the respective MAIN script. It will command a configuration class Config
that contains all the properties and methods required for the filter. In addition,
there are the filter functions which are separate files.

Scenario - two Targets cross at the exact same point in time and space. 
They are being tracked by a radar sensor, placed at the origin, which 
produces noisy range, azimuth, and elevation measurements, recorded at a 
rate of dt. The targets are assumed trackable at all times.

Note - this is a MTT example, but:
    . no births and deaths of true targets at different times
    . no spawning
    . no measurement gating
    . no labeling 

Enjoy :)