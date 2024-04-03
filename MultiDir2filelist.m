%MultiDir2filelist script.
%Ask for multiple files possibly in multiple directories, and for an output
%directory, and return:
%filelist: list of files organized as in output from dir;
%startdir: directory of the first selected group of files;
%rootdir: output directory
%[fnm,pnm,fim]: output as from MultiDirOpen

options.Resize='on';
options.WindowStyle='normal';
options.Interpreter='none';
format compact
format short g
if ~exist('titl','var'), titl='Select files (Cancel: files in the dirtree from output folder)'; end
if ~exist('strtfl','var'), strtfl='*.csv'; end

[fnm,pnm,fim]=MultiDirOpen(fileext('text+'),strtfl,titl,true);
if ~isempty(fim)
    fnm(~cellfun(@iscell,fnm))=cellstr(fnm(~cellfun(@iscell,fnm)));
    pnmt=cellfun(@(p,f) repmat({p},size(f)),pnm,fnm,'unif',false);
    pnmt=[pnmt{:}];
    fnmt=[fnm{:}];
    filelist=struct('name',fnmt,'folder',pnmt,'isdir',0);
    startdir=pnmt{1};
%     clear fnmt pnmt
else
    filelist=struct('name',{},'folder',{},'isdir',0);
    startdir=pwd;
end

% rootdir = uigetdir(startdir,'Output folder');
[strtfl,rootdir,fi]=uiputfile(fileext('text',true),'Output folder; input type from .ext',fullfile(startdir,strtfl));
if ~rootdir, return; end
if ~exist('outfile','var') || pref('get', 'Luos', 'raman.','singleoutfile',false)
    outfile=strtfl;
end
[~, ~, strtfl]=fileparts(strtfl);
if fi==2 || isempty(strtfl)
    strtfl='*.*';
else
    strtfl=['*' strtfl];
end
if ~exist('filelist','var') || isempty(filelist)
    %%
    filelist = dir(fullfile(rootdir, '**', strtfl));  %get list of files and folders in any subfolder
    filelist=filelist(~arrayfun(@(s) s.name(1)=='.',filelist));
%     dirlist = filelist([filelist.isdir]);
    filelist = filelist(~[filelist.isdir]);  %remove folders from list
    fnm={filelist.name};
    pnm={filelist.folder};
    fim=repmat({1},size(pnm));
end