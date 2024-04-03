% usrpth=userpath;
% usrpth(end)='';
addpath(genpath('C:\\Fiji.app\\scripts'));
usrpth=fileparts(mfilename('fullpath')); %#ok<GPFST>
cd(usrpth);
addpath(usrpth);
addpath(genpath(fullfile(usrpth,'utrack','u-track-2.1.3','software')));
addpath(genpath(fullfile(usrpth,'u-track-2.1.3','software')));
addpath(genpath(fullfile(usrpth,'LuosScript')));
rmpath(genpath(fullfile(usrpth,'Luosscript','Carmine')));
% addpath(genpath(fullfile(usrpth,'AutoUtrack')));
% addpath(genpath(fullfile(usrpth,'AutoTM')));
% addpath(genpath(fullfile(usrpth,'Carmine')));
% addpath(genpath(fullfile(usrpth,'partdrift')));
warning off MATLAB:dispatcher:nameConflict;
clear all;
close all;
fclose all;
clear java;
clear classes;
% if ispref('Luos'),rmpref('Luos');end
% if ispref('LuosTrack'),rmpref('LuosTrack');end
if pref('is','Luos'),pref('rm','Luos');end
if pref('is','LuosTrack'),pref('rm','LuosTrack');end
if pref('is','inpdlg'),pref('rm','inpdlg');end
global slash usrpth;
slash=filesep;
pref('set', 'Luos', 'slash', slash); 
usrpth=fileparts(mfilename('fullpath'));
pref('set', 'Luos', 'usrpth', usrpth); 

