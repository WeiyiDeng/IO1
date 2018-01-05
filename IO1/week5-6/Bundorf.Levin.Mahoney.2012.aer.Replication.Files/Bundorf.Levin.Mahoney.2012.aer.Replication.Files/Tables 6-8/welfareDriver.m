clear all; clc;

%generates draws of epsilon from conditional distribution; generates costs;
structureWelfareData;

%calculates net social surplus over grid space of incremental contributions;
welfareOpt;

%determines max net social surplus under various non descrimination
%constraints;
welfareMax;

%value of choice;
welfareNew;

%generates statistics for output
welfareStats;

%writes output
welfareOut;

