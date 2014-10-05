% THIS SCRIPT COMPILES THE BISECTION-SOURCE
[Bispathstrmain, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(Bispathstrmain),'-begin');
compileBisectionSource;
cd(Bispathstrmain);
