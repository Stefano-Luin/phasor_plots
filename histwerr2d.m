function res = histwerr2d(xdata,dmin,dmax,dx,varargin)

[cax,args,~] = axescheck(varargin{:});

ndata = size(xdata,1);
switch size(xdata,2)
    case 1
        xdata(:,[2,3,4,5])=repmat([0,0,0,1],ndata,1);
    case 2
        xdata=[xdata(:,1),zeros(ndata,1),xdata(:,2),zeros(ndata,1),ones(ndata,1)];
    case 3
        xdata=[xdata(:,1),zeros(ndata,1),xdata(:,2),zeros(ndata,1),xdata(:,3)];
    case 4
        xdata(:,5)=ones(ndata,1);
end
    
xdata=xdata((~any(isnan(xdata),2)&all(isfinite(xdata(:,[2,3])),2)),:);
ndata = size(xdata,1);

% Process input parameter name/value pairs
pnames = {'nbins','ctrs','edges','log1','log2','makeplot','halfhalf'};
dflts =  { [],     [],       [],    0,     0,      0,         0};
if verLessThan('matlab', '8')
    [errid,errmsg,nbins,ctrs,edges,islog1,islog2,makeplot,half,arg] = internal.stats.getargs(pnames, dflts, args{:});
    if ~isempty(errmsg)
        error(['Luin:histwerr2D:' errid], errmsg);
    end
else
    [nbins,ctrs,edges,islog1,islog2,makeplot,half,~,arg] = internal.stats.parseArgs(pnames, dflts, args{:});
end
if makeplot==2, plotArgs=arg; end

islog=logical([islog1,islog2]);
xdata(:,logical([islog1,0,islog2]))=log10(xdata(:,logical([islog1,0,islog2])));

if (isempty(nbins)+isempty(ctrs)+isempty(edges)) < 2
    error('LUIN:histwerr2D:AmbiguousBinSpec', 'Ambiguous specification of bins.');
end

if ndata    
    if ~isempty(nbins)
        arg = [{'nbins',nbins},arg];    
    elseif ~isempty(ctrs)
        if islog(1), ctrs{1}=log10(ctrs{1}); end
        if islog(2), ctrs{2}=log10(ctrs{2}); end
        arg = [{'ctrs',ctrs},arg];    
    elseif ~isempty(edges)
        if islog(1), edges{1}=log10(edges{1}); end
        if islog(2), edges{2}=log10(edges{2}); end
        arg = [{'edges',edges},arg];    
    else
        %calculate binning
        x=num2cell(xdata(:,[1,3]),1);
        for i=1:2

            if islog(i)
               x{i}=real(x{i}(imag(x{i})==0));
               if length(dmin)>=i
                   if dmin(i)>0,dmin(i)=log10(dmin(i));
                   else dmin(i)=-Inf; end;
               end;
               if length(dmax)>=i
                   if dmax(i)>0,dmax(i)=log10(dmax(i));
                   else dmax(i)=-Inf; end;
               end;
            end

            if (nargin<3 ||  length(dmin)<i || length(dmax)<i || dmin(i) == dmax(i) || ~isfinite(dmin(i)) || ~isfinite(dmax(i)))
                if isempty(x{i}(isfinite(x{i})))
                    dmin(i)=0;
                    dmax(i)=0;
                else
                    dmin(i)=min(x{i}(isfinite(x{i}(:,1)),1));
                    dmax(i)=max(x{i}(isfinite(x{i}(:,1)),1));
                end
                minmaxgiven=false;
            else
                minmaxgiven=true;
            end

            if dmax(i)<dmin(i)
                temp=dmax(i);
                dmax(i)=dmin(i);
                dmin(i)=temp;
                clear temp;
            end

            if (nargin<4 || length(dx)<i || dx(i) == 0 || ~isfinite(dx(i)))
                dx(i)=(dmax(i)-dmin(i))./sqrt(ndata);
            end

            if dx(i)==0
                if islog(i), dx(i)=mean(xdata(:,2*i));
                else dx(i)=mean(xdata(:,2*i))/mean(abs(10.^x{i})); end;
            end;
            if dx(i)==0, dx(i)=1; end; %Something strange in the data...

            edges{i}=(dmin(i)-dx(i)/2.0:dx(i):dmax(i)+dx(i)/2.0);
            ctrs{i}=mean([[dmin(i)-1.5*dx(i), edges{i}];[edges{i},edges{i}(end)+dx(i)]]);
            if minmaxgiven
                edges{i}=edges{i}(2:end-1);
                ctrs{i}=ctrs{i}(2:end-1);
            end
        end
        arg = [{'ctrs',ctrs},arg];    
    end
    if (makeplot && makeplot~=2), arg = [arg,{'makeplot',1}]; end;    
    arg=[{xdata},arg,{'log1',islog1,'log2',islog2,'halfhalf',half}];
    if ~isempty(cax), arg=[{cax},arg]; end;
    [nn,ctrs,edges]=hist3werr(arg{:});
else
    nn=0;ctrs={0 0};edges={[0;0] [0;0]};
end
if islog1
    ctrs{1}(:,2)=real(10.^ctrs{1}(:,1));
    edges{1}(:,2)=real(10.^edges{1}(:,1));
end;
if islog2
    ctrs{2}(:,2)=real(10.^ctrs{2}(:,1));
    edges{2}(:,2)=real(10.^edges{2}(:,1));
end;

res.rmat=nn;
res.r=[ctrs,edges];
if makeplot == 2
    drawhist2d(ctrs,edges,nn,{'x','y'},'hist2D',plotArgs);
end
end
