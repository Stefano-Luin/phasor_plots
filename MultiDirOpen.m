function [FileName,PathName,FilterIndex]=MultiDirOpen(f,startdir,titl,isfile,savepp)
l=1;
k=1;
if nargin<5 || isempty(savepp)
    savepp=false;
end
if nargin<4 || isempty(isfile)
    isfile=false;
end
if nargin<3 || isempty(titl)
    titl='Open... (cancel to stop)';
else
    titl=[titl '... (Cancel to stop)'];
end
try
    nmax=30;
    maxchar=255;
    FileName={};
    PathName={};
    FilterIndex={};
    
    if nargin<1||isempty(f)
        f=fileext('traj');
    end
    if nargin<2 || isempty(startdir)
        startdir=pwd;
    end
    ppath=startdir;
    if isfile
        suffix='';
    else
        suffix='\';
    end
    while 1
        stp='No';
        fname={};
        while strcmpi(stp,'No')
            drawnow; pause(0.01);
            if isfile
                DefaultName=fullfile(ppath,startdir);
            else
                DefaultName=[ppath suffix];
            end
            [fn,pname,findex] = uigetfile(f(:,1:2),titl,DefaultName,'MultiSelect','on');
            if isequal(fn,0)
                stp=questdlg('Stop?','Multi Dir Open','Yes','No','Yes and separate','Yes');
                %if strcmpi(stp,'Cancel'), stp='No'; end
            else
                ppath=pname;
                if size(f,2)>2, findex=f{findex,3}; end
                fname=[fname,fn];
                display(['first file: ' fname{1}]);
                display(['last file: ' fname{end}]);
                display([int2str(length(fname)) ' files; ']);
                if length(fname)>nmax || any(cellfun(@length,fname)>maxchar)
                    altt=input('Ok? (Y/n) ','s');
                    if (isempty(altt) || altt(1)=='y' || altt(1)=='Y')
                        break;
                    else
                        %fname(end)=[];
                        %display('select remaining files, including "last file" even if correct');
                        display('select remaining files');
                    end
                else
                    break;
                end
            end
        end
        
        if (stp(1)=='y' || stp(1)=='Y')
            break;
        elseif findex==99
            for kk=1:length(fname)
                temp=load(fullfile(pname,fname{kk}));
                if isfield(temp,'FileName')
                    fn='FileName';
                    pn='PathName';
                    fi='FilterIndex';
                elseif isfield(temp,'pn')
                    fn='fn';
                    pn='pn';
                    fi='fi';
                else
                    disp(['Something wrong in file ' fname{kk} ' in folder ' pname]);
                    continue;
                end
                l=length(temp.(fi));
                for ii=1:l
                    fprintf('%d %s \n',ii,temp.(pn){ii});
                end
                %%
                if l>1
                    a=input('Start from file number ');
                    if ~isscalar(a)||~isnumeric(a)||~isfinite(a)
                        a=1;
                    end
                    if a<=0,a=a+1;end
                    a=mod(int16(a),l);
                    if ~a, a=l; end
                    ind=[a:l,1:(a-1)];
                else
                    ind=1;
                end
                %%
                if l
                    FileName=[FileName,temp.(fn)(ind)];
                    PathName=[PathName,temp.(pn)(ind)];
                    FilterIndex=[FilterIndex,temp.(fi)(ind)];
                    k=k+l;
                end
            end
        else
            FileName{k}=fname; %#ok<*AGROW>
            PathName{k}=pname;
            FilterIndex{k}=findex;
            k=k+1;
        end
    end
catch ME
    disp(getReport(ME));
%     keyboard;
end

if strcmp('Yes and separate',stp)
    for k=1:length(FilterIndex)
        if ~iscell(FileName{k})
            FileName{k}=FileName(k);
        end
        PathName{k}=repmat(PathName(k),size(FileName{k}));
        FilterIndex{k}=repmat(FilterIndex(k),size(FileName{k}));
    end
    PathName=[PathName{:}];
    FilterIndex=[FilterIndex{:}];
    FileName=[FileName{:}];
    FileName=cellfun(@(c) {c},FileName,'Unif',false);
end

drawnow; pause(0.02);
if savepp || k>l+1 && size(f,2)>2 && any(cellfun(@(c) c==99,f(:,3)))
    [fout,pout]=uiputfile('*.mat','Save file list?',['pp', datestr(now,'yymmddHHMMSS'),'.mat']);
    if ~isequal(fout,0)
        save(fullfile(pout,fout),'FileName','PathName','FilterIndex');
    end
    drawnow; pause(0.05);
else
    fout=0;
    pout=0;
end
pref('set','LuosTrack','pout',pout);
pref('set','LuosTrack','fout',fout);
