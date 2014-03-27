
Command file for recon3d


Number of Scenes (-)                                        ::  9 
Master Scene (-)                                            ::  1

Radar Wavelength (in preferred output distance units)       ::  0.24
Time converstion factor (applied to all scenes equally)     ::  365 
Distance converstion factor (applied to all scenes equally) ::  none

Output folder                       :: . 
Output file prefix                  :: synthetic_test_1

Output diag(G^T*W*G)^-1)            :: y
Output offdiag(G^T*W*G)^-1)         :: y
Output diag(G^T*G)^-1               :: y
Output offdiag(G^T*G)^-1            :: y
Output GDOP sqrt(tr((G^TG)^-1))     :: y
Output sqrt(tr((G^T*W*G)^-1))       :: y 
Output mean square error estimate   :: y

Output velocity magnitude           :: y
Output average correlation          :: y

Output number of scenes             :: y
Output rank(G) at each pixel        :: y
 
Output null value                   :: -100
Prior model covariance weighting factor :: 1.e-4

Scene :: 1
   Unwrapped or LOS displacement file           :: Synth_32005_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 223
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_32005_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_32005_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.1400447158
   Upper left corner Longitude (deg)         :: -19.3893619348
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000



Scene :: 2

   Unwrapped or LOS displacement file           :: Synth_02001_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 161
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_02001_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_02001_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.0579829482
   Upper left corner Longitude (deg)         :: -19.0748096276
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000


Scene :: 3

   Unwrapped or LOS displacement file           :: Synth_07500_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 170
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0
 
   Correlation file                          :: Synth_07500_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_07500_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 64.910140541
   Upper left corner Longitude (deg)         :: -19.2250624235
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555600017833

   Number of looks                           :: 90000


Scene :: 4

   Unwrapped or LOS displacement file           :: Synth_07502_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 170
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_07502_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_07502_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.0449290542
   Upper left corner Longitude (deg)         :: -19.2942839201
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555600017833

   Number of looks                           :: 90000



Scene :: 5
   Unwrapped or LOS displacement file           :: Synth_14009_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 222
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_14009_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_14009_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.0936503676
   Upper left corner Longitude (deg)         :: -19.75213583
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000



Scene :: 6
   Unwrapped or LOS displacement file           :: Synth_20001_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 165
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_20001_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_20001_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.173987892
   Upper left corner Longitude (deg)         :: -19.2723632623
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000



Scene :: 7
   Unwrapped or LOS displacement file           :: Synth_20003_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 164
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_20003_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_20003_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.1623759536
   Upper left corner Longitude (deg)         :: -19.5280275877
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555600017833

   Number of looks                           :: 90000




Scene :: 8
   Unwrapped or LOS displacement file           :: Synth_32006_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 223
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_32006_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_32006_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.2287745742
   Upper left corner Longitude (deg)         :: -19.1414752771
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000



Scene :: 9
   Unwrapped or LOS displacement file           :: Synth_25501_12039_12040_HH.unw.grd 
   Samples in Unwrapped or LOS displacement (-) :: 170
   Type of unwrapped or LOS displacement file (enter phase or disp) :: phase
   Phase or displacement null value          :: 0

   Correlation file                          :: Synth_25501_12039_12040_HH.cor.grd
   Correlation null value                    :: 0

   LOS file (must be by-pixel interleave in order ENU) :: Synth_25501_12039_12040_HH.los.grd
   LOS null value                            :: -100

   Upper left corner Latitude (deg)          :: 65.1738810805
   Upper left corner Longitude (deg)         :: -19.3073306848
   Latitude Spacing (deg/pixel)              :: -0.00277800008917
   Longitude Spacing (deg/pixel)             :: 0.00555500009796

   Number of looks                           :: 90000







End of command file





