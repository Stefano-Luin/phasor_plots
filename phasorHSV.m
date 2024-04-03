function [hs,rgb,satp,cntr]=phasorHSV(satp,pxy,cntr,spar,spar2,val,doplot)
%   function [hs,rgb]=phasorHSV(satp,pxy,cntr,spar,spar2,val,doplot)
% Input:
%   satp: "principal" points, where saturation has to be (close to) maximum
%   pxy: 2-column vector, points where hue and saturation must be
%        calculated. If it is a cell array, it contatins two vectors with
%        the x (g) and y (s) coordinates of a grid in the phasor plot.
%   cntr: center for angles for hue values, point where saturation=0;
%         default: average of satp.
%   spar <0.5>; %saturation parameter: minimum saturation at sat points or
%         semicircle borders
%   spar2 <0.75>; %saturation parameter: maximum saturation at sat points
%         or semicircle borders
%   val: vector with intensities/values on points pxy (must have the right
%        length; default all 1)
%   doplot: if true, make a plot in the phasor plane
%
% Output:
%   hs: n x 2 double array [hue, saturation] for each of the n "pxy"
%       points in each row;
%   rgb: n x 3 double array of RGB values for corresponding
%       [hue, saturation, val] for each of the n "pxy" points in each row;
%   satp: principal points
%   center: center point
%
%   To do: consider to make sat=sat^a, spar[2]=spar[2]^(1/a)
%   To do: correct what's happening if cntr is on (or close to) a line
%   connecting two satp or on semicircle borders, especially if both.

if nargin<1 || isempty(satp)
    %satp=[0.5,0.45;0.7,0.3;0.7,0.45;0.1,0.0];
    %satp=[0.5,0.45;0.7,0.3;0.1,0.0];
    %satp=[0.1,0.1;0.5,0.49;0.9,0.02];
    satp=[0.0,0.0;0.5,0.5;1,0.0];
end
if nargin<2 || isempty(pxy)
    % pxy=[0.4999,0.45;0.7,0.30001;0.70001,0.45;0.1,0.0001];
    [X,Y]=ndgrid(-0.1:0.02:1.1,-0.1:0.02:0.6);
    pxy=[X(:),Y(:)];
    %pxy=[0.3,0.4;0.4,0.2;0.6,0.01;0,0];
    %pxy=[0.4999,0.45;0.7,0.30001;0.70001,0.45;0.1,0.0001];
end
if iscell(pxy)
    [X,Y]=ndgrid(pxy{1}(:),pxy{2}(:));
    pxy=[X(:),Y(:)];
end
if nargin<3 || isempty(cntr)
    cntr=mean(satp);
    %cntr=[0.5,0.166];
end
szsat=size(satp,1);
lp=size(pxy,1);
if nargin<4 || isempty(spar)
    spar=0.5; %saturation parameter: minimum saturation at sat points or semicircle borders
end
if nargin<5 || isempty(spar2)
    spar2=0.75; %saturation parameter: maximum saturation at sat points or semicircle borders
end

if nargin<7 || isempty(doplot)
    doplot=false;
end

if (doplot || nargout>1) && (nargin<6 || isempty(val))
    val=ones(lp,1);
end

if doplot
    figure
    %scatter(satp(:,1),satp(:,2),'.')
    scatter(satp(:,1),satp(:,2),'SizeData',25.*(1:szsat))
    axis equal;
    ax = gca;
    ax.Color = 'k';
%     ax.XLim=[-0.15,1.15];
%     ax.YLim=[-0.15,0.65];
    ax.XAxisLocation='origin';
    ax.YAxisLocation='origin';
    ax.XColor='w';
    ax.YColor='w';
    ax.Layer = 'top';
    hold on;
    t=0:0.01:1;
    plot(cos(pi*t)/2+1/2,sin(pi*t)/2,'w');
    plot(t,zeros(size(t)),'w');
    %scatter(pxy(:,1),pxy(:,2),'r','SizeData',25.*(1:lp))
    scatter(cntr(:,1),cntr(:,2),'+w')
end



[phc,rc]=cart2pol(1/2-cntr(1),0-cntr(2));

[phs,rs]=cart2pol(satp(:,1)-cntr(1),satp(:,2)-cntr(2));
[php,rp]=cart2pol(pxy(:,1)-cntr(1),pxy(:,2)-cntr(2));
[phssrt,indsrt]=sort(phs);
cs=find(indsrt==1,1);
satp(:,:)=satp(circshift(indsrt,1-cs),:);
phsdst=diff([-pi;phssrt;pi]);
if ~all([phsdst(2:end-1);phsdst(1)+phsdst(end)])
    error('Luin:phasorHSV','phasorHSV: aligned points');
end
phpsatdif=php-[-pi;phssrt;pi]';
indp=sum(phpsatdif>=0,2); %index of the last positive difference, minimum 1; along the columns if more than one pxy, to be checked.
indp(indp>szsat+1)=szsat+1;
distph=NaN(lp,2);
for c= 1:lp
    if isnan(indp(c)) || ~indp(c)
        continue
    end
    distph(c,:)=[phpsatdif(c,indp(c)),-phpsatdif(c,indp(c)+1)];
    if indp(c)==1
        distph(c,1)=distph(c,1)+phsdst(end);
    elseif indp(c)==szsat+1
        distph(c,2)=distph(c,2)+phsdst(1);
    end
end
%%
satfirst=find(indsrt==1);
ord=indsrt(mod(satfirst,szsat)+1)<indsrt(mod(satfirst-2,szsat)+1);
clear indcol;
%indcol=mod(satfirst-1+(2*ord-1)*(0:(szsat-1)),szsat)+1;
indcol(mod(satfirst-1+(2*ord-1)*(0:(szsat-1)),szsat)+1)=1:szsat;
% indcol(end+1)=1;
indp(indp==szsat+1)=1;
hue=NaN(lp,1);
sat=hue;
rfromline=hue;
rfromcirc=hue;
rfromx=hue;
for c=1:lp
    if isnan(indp(c)) || ~indp(c)
        continue
    end
    pre=indcol(mod(indp(c)-2,szsat)+1);
    post=indcol(mod(indp(c)-1,szsat)+1);
    if pre==1 && post==szsat
        pre=szsat+1;
    elseif post==1 && pre==szsat
        post=szsat+1;
    end
    hue(c)=(distph(c,1)*(post-1)+distph(c,2)*(pre-1))/(szsat)/(distph(c,1)+distph(c,2));
    %hue(c)=(sin(distph(c,1))*(post-1)+sin(distph(c,2))*(pre-1))/(szsat)/(sin(distph(c,1))+sin(distph(c,2)));
    
    %     ppre=pre;
    %     ppost=post;
    ppre=indcol(mod(indp(c)-2,szsat)+1);
    ppost=indcol(mod(indp(c)-1,szsat)+1);
    
    rfromline(c)=rs(ppre)/(cos(distph(c,1))+(rs(ppre)/rs(ppost)-cos(distph(c,1)+distph(c,2)))*sin(distph(c,1))/sin(distph(c,1)+distph(c,2)));
    rfromcirc(c)=rc*cos(php(c)-phc)+sqrt(1/4.0-rc^2*(sin(php(c)-phc))^2);
    %rfromcircm(c)=rc*cos(php(c)-phc)-sqrt(1/4.0-rc^2*(sin(php(c)-phc))^2);
    rfromx(c)=cntr(2)/sin(-php(c));
    if ~isreal(rfromcirc(c)) || rfromcirc(c)<=0
        if ~isfinite(rfromline(c)) || rfromline(c)<0
            if isfinite(rfromx(c))&&rfromx(c)>0
                sat(c)=rp(c)/rfromx(c);
            else
                sat(c)=1-exp(-rp(c)/(-rfromline(c)));
            end
        elseif rp(c)<rfromline(c)
            sat(c)=spar2*rp(c)/rfromline(c);
        else
            sat(c)=spar2+(1-spar2)*(1-exp(-(rp(c)/rfromline(c)-1)));
        end
    elseif rfromcirc(c)<rfromx(c) || rfromx(c)<=0 % rfromcirc(c)>0
        if ~isfinite(rfromline(c)) || rfromline(c)<0
            sat(c)=spar2*rp(c)/rfromcirc(c);
        elseif rfromcirc(c)<rfromline(c)
            alpha=spar+(spar2-spar)*(rfromcirc(c)/rfromline(c))^2;
            if rp(c)<rfromcirc(c)
                %sat(c)=alpha*rp(c)/rfromcirc(c);
                sat(c)=alpha*1/sin(1)*sin(rp(c)/rfromcirc(c));
            elseif rp(c)>rfromline(c)
                sat(c)=spar2+(1-spar2)*(1-exp(-2*(rp(c)-rfromline(c))/(rfromcirc(c)+rfromline(c))));
            else
                sat(c)=alpha+(spar2-alpha)*(rp(c)-rfromcirc(c))/(rfromline(c)-rfromcirc(c));
            end
        else %rfromcirc(c)>rfromline(c)
            alpha=spar+(spar2-spar)*(rfromline(c)/rfromcirc(c))^2;
            if rp(c)<rfromline(c)
                %sat(c)=alpha*rp(c)/rfromline(c);
                sat(c)=alpha*1/sin(1)*sin(rp(c)/rfromline(c));
            elseif rp(c)>rfromcirc(c)
                sat(c)=spar2+(1-spar2)*(1-exp(-2*(rp(c)-rfromcirc(c))/(rfromcirc(c)+rfromline(c))));
            else
                sat(c)=alpha+(spar2-alpha)*(rp(c)-rfromline(c))/(rfromcirc(c)-rfromline(c));
            end
        end
    else % rfromcirc(c)>rfromx(c)>0
        if ~isfinite(rfromline(c)) || rfromline(c)<0
            sat(c)=spar2*rp(c)/rfromx(c);
        elseif rfromx(c)<rfromline(c)
            alpha=spar+(spar2-spar)*(rfromx(c)/rfromline(c))^2;
            if rp(c)<rfromx(c)
                %sat(c)=alpha*rp(c)/rfromx(c);
                sat(c)=alpha*1/sin(1)*sin(rp(c)/rfromx(c));
            elseif rp(c)>rfromline(c)
                sat(c)=spar2+(1-spar2)*(1-exp(-2*(rp(c)-rfromline(c))/(rfromx(c)+rfromline(c))));
            else
                sat(c)=alpha+(spar2-alpha)*(rp(c)-rfromx(c))/(rfromline(c)-rfromx(c));
            end
        else %rfromx(c)>rfromline(c)
            alpha=spar+(spar2-spar)*(rfromline(c)/rfromx(c))^2;
            if rp(c)<rfromline(c)
                %sat(c)=alpha*rp(c)/rfromline(c);
                sat(c)=alpha*1/sin(1)*sin(rp(c)/rfromline(c));
            elseif rp(c)>rfromx(c)
                sat(c)=spar2+(1-spar2)*(1-exp(-2*(rp(c)-rfromx(c))/(rfromx(c)+rfromline(c))));
            else
                sat(c)=alpha+(spar2-alpha)*(rp(c)-rfromline(c))/(rfromx(c)-rfromline(c));
            end
        end
    end
    if sat(c)>1, sat(c)=1; end
    if sat(c)<0 %should not be necessary
        sat(c)=0;
        warning('Luin:phasorHSV','phasorHSV: negative saturation found');
        %keyboard;
    end
end
hs=[hue,sat];
if doplot || nargout>1
    rgb=hsv2rgb([hs,val]);
end
if doplot
    scatter(pxy(:,1),pxy(:,2),[],rgb,'.');
    scatter(cntr(:,1),cntr(:,2),'+w')
end

%     if ~isreal(rfromcirc(c)) || rfromcirc(c)<0
%         if ~isfinite(rfromline(c)) || rfromline(c)<0
%             if isfinite(rfromx(c))&&rfromx(c)>0
%                 sat(c)=rp(c)/rfromx(c);
%             else
%                 sat(c)=exp(-2*rp(c));
%             end
%         elseif rp(c)<rfromline(c)
%             sat(c)=spar*rp(c)/rfromline(c);
%         else
%             sat(c)=spar+(1-spar)*exp(-(rp(c)/rfromline(c)-1));
%         end
%     elseif rfromcirc(c)<rfromx(c) || rfromx(c)<0 % rfromcirc(c)>0
%         if ~isfinite(rfromline(c)) || rfromline(c)<0
%             sat(c)=rp(c)/rfromcirc(c);
%         elseif rfromcirc(c)<rfromline(c)
%             sat(c)=rp(c)/rfromline(c);
%         elseif rp(c)<rfromline(c)
%             sat(c)=spar*rp(c)/rfromline(c);
%         else
%             sat(c)=spar+(1-spar)*(rp(c)-rfromline(c))/(rfromcirc(c)-rfromline(c));
%         end
%     else % rfromcirc(c)>rfromx(c)>0
%         if ~isfinite(rfromline(c)) || rfromline(c)<0
%             sat(c)=rp(c)/rfromx(c);
%         elseif rfromx(c)<rfromline(c)
%             sat(c)=rp(c)/rfromline(c);
%         elseif rp(c)<rfromline(c)
%             sat(c)=spar*rp(c)/rfromline(c);
%         else
%             sat(c)=spar+(1-spar)*(rp(c)-rfromline(c))/(rfromx(c)-rfromline(c));
%         end
%     end