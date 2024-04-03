function h=drawhist2d(ctrs,edges,nn,name,fname,plotArgs,cax)
%h=drawhist2d(ctrs,edges,nn,name,fname,plotArgs,cax)

% Build xy-coords.
[xx,yy] = ndgrid(ctrs{1}(:,1),ctrs{2}(:,1));
if nargin<7,cax=[];end
if nargin<6,plotArgs={};end
cax = newplot(cax);
holdState = ishold(cax);

% Plot the surface, using any specified graphics properties to override
% defaults.
h = surfc(cax, xx, yy, nn, 'FaceColor','interp',plotArgs{:});
set(cax,'Title',text('String', fname));
set(cax,'XLabel',text('String',name{1}))
set(cax,'YLabel',text('String',name{2}))

if ~holdState
    % Set ticks for each bar if fewer than 16 and the centers/edges are
    % integers.  Otherwise, leave the default ticks alone.
    if (length(ctrs{1}(:,1))<16) && all(floor(ctrs{1}(:,1))==ctrs{1}(:,1))
        set(cax,'xtick',ctrs{1}(:,1));
    end
    if (length(ctrs{2}(:,1))<16) && all(floor(ctrs{2}(:,1))==ctrs{2}(:,1))
        set(cax,'ytick',ctrs{2}(:,1));
    end
    
    % Set the axis limits to have some space at the edges.
    dx = range(edges{1}(:,1))*.05;
    dy = range(edges{2}(:,1))*.05;
    set(cax,'xlim',[edges{1}(1,1)-dx edges{1}(end,1)+dx]);
    set(cax,'ylim',[edges{2}(1,1)-dy edges{2}(end,1)+dy]);
    
    view(cax,3);
    grid(cax,'on');
    set(get(cax,'parent'),'renderer','zbuffer');
end