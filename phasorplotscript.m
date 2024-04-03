% Script to automatize the use of the function "phasor_plots.m" using the
% same parameters for many files.
% Starts asking multiple files also in multiple directories, then creates
% the image+phasor plots in linear and logarithmic scale.
% Define the principal points (satp, with coordinates along columns and
% different points along rows) and eventually the center point (if
% satp is a cell array of length 2, {[principal points], [center point]}, a
% suffix for the file name (popt.suff), the width of the gaussian used for
% smoothing images and phasors (popt.sgauss), the two saturation parameters
% (popt.spar), the minimum for the logarithmic scale (logI) and how to draw
% the principal and central points in the phasor plot (popt.drawsatp - see 
MultiDir2filelist;
if isempty(filelist)
    %%
    filelist.name='';
    filelist.folder='';
end
%%
warnsit=warning('off','LUIN:wmean1');
% satp=[0.68,0.433;0.992,0.078;0.24,0.265];%triangle DMF
% popt.suff='triangle';
% satp=[0.45,0.362;0.54,0.335;0.525,0.195];%cells
% popt.suff='cells';
satp={[0,0;0.07696,0.07792;0.24,0.265;0.35823,0.34331;0.68,0.433;0.992,0.078],[0.4,0.06]};%NPs
popt.suff='NPs';
popt.sgauss=1;
popt.spar=[0.85,0.95];
popt.drawsatp=-3;
logI=0.01;
% logI=0.005;
for cf=1:length(filelist)
    if iscell(filelist(cf).name)
        filelist(cf).name=filelist(cf).name{1};
    end
    display([fullfile(filelist(cf).folder,filelist(cf).name) ' ' datestr(now)]);
    if mod(cf,25)==0 && length(filelist)-cf > 20
        close all;
    end
    popt.logI=0;
    %[hs,rgb,limxy,hsim,rgbim,h2, Filename,Pathname,Ithvec]=...
    [~,~,~,~,~,~,filelist(cf).name,filelist(cf).folder]=phasor_plots([],filelist(cf).name,filelist(cf).folder,2,0,satp,[],popt);
    popt.logI=logI;
    phasor_plots([],filelist(cf).name,filelist(cf).folder,2,0,satp,[],popt);
    %phasor_plots([],'','',2,1,satp);(c,Filename,Pathname,doplot,doscatter,satp,Ith,sgauss)
end
warning(warnsit);