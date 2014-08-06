% THIS SCRIPT COMPILES THE BISECTION-SOURCE
[pathstrmain, name, ext] = fileparts(mfilename('fullpath'));
addpath(genpath(pathstrmain),'-begin');
compileBisectionSource;
cd(pathstrmain);
