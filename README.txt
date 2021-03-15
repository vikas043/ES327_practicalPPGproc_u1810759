This repository contains all the work done as a part of ES327 Technical Project at University of Warwick.
It contains 14 data files from two datasets labelled with the subject IDs from their originating database 
and in the same style in their originating database. Files starting with 'DATA' were obtained from IEEE SP Cup 2015.
Other data files from Uni. of Marche database on PPG. The task was to create a signal processing workflow
and test it across these datasets. For this 5 function files are provided, each output a struct with performance metrics
and figures showing signal charactersitics at each stage along with some metrics and filter properties. Some of the key figures
have also been stored in the sub-folder labelled figures. There are four files meant for datasets,
ones ending in SP for SP cup.
ones ending in Mar for Marche datasets.
The user can simply comment and uncomment as directed within to change subjects.
Because of lack of expected functionality, a test was also created on known sinusoids and noise.
This is described by the AdapTest function and again produces figures at each stage of processing.
User will require access to MATLAB to run these functions and view the figure outputs,
ones not already provided as images.
