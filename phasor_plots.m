function [hs,rgb,limxy,hsim,rgbim,h2, Filename,Pathname,Ithvec]=...
    phasor_plots(c,Filename,Pathname,doplot,doscatter,satp,Ith,ppopt)
% [hs,rgb,hsim,rgbim,Filename,Pathname,Ithvec]=...
%    phasor_plots(c,Filename,Pathname,doplot,doscatter,satp,Ith,sgauss)
% Plot image and phasor plot with hue and saturation depending on the
% position on phasor plot and intensity as for the image in the image, and
% as the sum of the intensities of points within a 2D-bin in the phasor
% plot.
%Input:
% c is an m x n x 3 (min) hypermatrix, with m x n image dimension, and
% intensity, and phase and modulation of the phasors, in the three 2D
% matrices.
% If c is not given, open one or more .R64 file(s) with given (cell of)
% Filename(s) and Pathname, or ask for them.
% If c is given, Filename and Pathname are the file name(s) and the path
% for the output figure (see below).
% doplot <true>: draw the figures with image and phasorplot. If doplot>1,
% save it.
% doscatter <true>: show the position of the phasor for each pixel.
% satp: n x 2 array, with principal points  coordinates along columns and
% different points along rows, and eventually the center point (if
% satp is a cell array of length 2: {[principal points], [center point]})
% Ith <0>: intensity threshold; NaN for asking for each file. If vector,
% the first two entries are used to determine a minimum and a maximum
% intensities (included).
% ppopt: struct with fields:
%       spar <[0.75,9]>: [spar, spar2] for phasorHSV;
%       sgauss <0>: width of gaussian for smoothing of intensity or
%           weighted with intensity
%       logI <0>: if different from 0, use a logarithmic intensity scale
%           between logI and 1.
%       drawsatp <0>: absolute value: 0: don't draw central points; 1 (or
%           any number different from 0, 2, or 3): draw segments joining
%           the principal points; 2: draw crosses on the principal points;
%           3: both; if negative: draw a cyan x on the central point.
%
% output (if applicable, they are cell arrays with the same information if
%    more than one phasor hypermatrix or file have been considered):
%
%    hs: nhist x 2 double array [hue, saturation] for each of the nhist 
%           points in the FLIM image, a row for each point;
%    rgb: n x 3 double array of RGB values for corresponding
%        [hue, saturation, val] for each of the n "pxy" points in each row;
%    limxy: limits of the 2D histogram in the phasor plot
%    hsim: m x n x 2 array with hue and saturation for the FLIM image
%    rgbim: m x n x 3 array encoding in RGB the colors of the FLIM image
%
%    h2: handle of the image containing image and phasor plot
%    Filename, Pathname: file name and path for the file containing FLIM
%       information for an image.
%    Ithvec: min accepted intensity or extremes for acceptable
%       intensities.

sgauss=0;
logI=0;
spar=0.75;
spar2=0.9;
drawsatp=0;
suff='';

if nargin>=8 && ~isempty(ppopt) && isstruct(ppopt)
    if isfield(ppopt,'sgauss') && ~isempty(ppopt.sgauss)
        sgauss=ppopt.sgauss;
    end
    if isfield(ppopt,'logI') && ~isempty(ppopt.logI)
        logI=ppopt.logI;
    end
    if isfield(ppopt,'spar') && ~isempty(ppopt.spar)
        spar=ppopt.spar(1);
        spar2=ppopt.spar(2);
    end
    if isfield(ppopt,'drawsatp') && ~isempty(ppopt.drawsatp)
        drawsatp=ppopt.drawsatp;
    end
    if isfield(ppopt,'suff') && ~isempty(ppopt.suff)
        suff=ppopt.suff;
    end
end
if nargin<2
    Filename='';
end
if nargin<3
    Pathname='';
end
if nargin<7 || isempty(Ith)
    Ith=0;
end
askTh=false;
Ithvec=Ith;
if isnan(Ith)
    askTh=true;
    Ithvec(1,2)=Inf;
end
if nargin<4 || isempty(doplot)
    doplot=true;
end
if nargin<5 || isempty(doscatter)
    doscatter=true;
end
cntr=[];
if nargin<6 || isempty(satp)
    satp=[];
elseif iscell(satp)
    cntr=satp{2};
    satp=satp{1};
end
if doplot || doscatter
%     hold on;
    t=0:0.01:1;
    %     plot(cos(pi*t)/2+1/2,sin(pi*t)/2,'k');
    %     plot(t,zeros(size(t)),'k');
end
if nargin<1||isempty(c)
    [c,Filename,Pathname] = Open_R64(Filename,Pathname); %put "Open_R64.m" in the same folder as this script before running
end
% cgr=colororder50(length(c));
if ~iscell(c)
    c={c};
end
if ~iscell(Filename)
    Filename={Filename};
end
if length(c)>1
    hsout=cell(length(c),1);
    rgbout=cell(length(c),1);
    hsimout=cell(length(c),1);
    rgbimout=cell(length(c),1);
    Ithvecout=cell(length(c),1);
    h2out=cell(length(c),1);
    limxyout=cell(length(c),1);
end
for i=1:length(c)
    if doplot>2, close all; end
    matr = c{i}; %extracts the content of each cell to a 3-dimensional matrix
    
    intensity = matr(:,:,1); %extracts the various images from the matrix: intensity
    har1phase = matr(:,:,2); %phase of the 1st harmonic
    har1mod = matr(:,:,3); %modulation of the 1st harmonic
    %     har2phase = matr(:,:,4); %phase of the 2nd harmonic
    %     har2mod= matr(:,:,5); %modulation of the 2nd harmonic
    
    if askTh
        intt=intensity;
        
        hint=figure;
        Ok=false;
        while(~Ok)
            hold off;
            imagesc(intt);
            disp(fullfile(Pathname,Filename{i}));
            Ith=input('Threshold ');
            figure(hint);hold on;
            if numel(Ith)==1
                Ith(1,2)=Inf;
            end
            image(repmat(zeros(size(intt)),[1,1,3]),'AlphaData',intt<=Ith(1)|intt>=Ith(2));
            Ok=input('Ok? (0/1) ');
        end
        Ithvec(i,:)=Ith;
    end
    
    har1phasor_g = har1mod.*(cosd(har1phase)); %obtains the g coordinate by multiplying elementwise modulation and cosine of the phase
    har1phasor_s = har1mod.*(sind(har1phase)); %obtains the s coordinate by multiplying elementwise modulation and sine of the phase
    
    [m,n] = size(har1phasor_g);
    hsim=NaN(m,n,2);
    rgbim=NaN(m,n,3);
    chos_g = har1phasor_g;
    chos_s = har1phasor_s;
    chos_I =intensity;
    
    if numel(Ith)==1
        chosen=intensity>=Ith;
    else
        chosen=intensity>=Ith(1)&intensity<=Ith(2);
    end
    chos_g(~chosen)=NaN;
    if all(isnan(chos_g(:)))
        warning('Luin:phasor_from_R64',['No intensities above threshold for ' Filename{i}]);
        %         [av_gfw(i),av_sfw(i),parw(i,1:2),mingf(i),maxgf(i),av_gf(i),av_sf(i),par(i,1:2),Dgf(i)]=deal(NaN);
        continue
    end
    
    chos_s(~chosen)=NaN;
    chos_I(~chosen)=NaN;
    if ~isempty(sgauss) && isfinite(sgauss) && sgauss>0
        chos_I=imageWA(chos_I,sgauss);
        chos_g=imageWA(chos_g,sgauss,chos_I);
        chos_s=imageWA(chos_s,sgauss,chos_I);
    end

    if doscatter || nargin<6 || isempty(satp)%&& i==1
        h=figure;
        hold on;
        %scatter(vec_gf,vec_sf,'.','MarkerEdgeColor',[0.8,0.8,0.8]);
        binscatter(chos_g(:),chos_s(:),'ShowEmptyBins','on');
        axis xy;
        plot(cos(pi*t)/2+1/2,sin(pi*t)/2,'w');
        plot(t,zeros(size(t)),'w');
        axis equal;
        ax = gca;
        ax.Color = 'k';
        ax.XLim=[-0.1,1.1];
        ax.YLim=[-0.1,0.6];
        ax.XAxisLocation='origin';
        ax.YAxisLocation='origin';
        ax.XColor='w';
        ax.YColor='w';
        ax.Layer = 'top';
    end
    limxy=[[max(min(chos_g(:)),-0.2),min(max(chos_g(:)),1.2)];...
        [max(min(chos_s(:)),-0.2),min(max(chos_s(:)),0.7)]];
    dxy=min(diff(limxy')./[201,101],[0.005,0.005]);
    %figure;
    res=histwerr2d([chos_g(:),chos_s(:),chos_I(:)],limxy(:,1),limxy(:,2),dxy,'makeplot',0,'halfhalf',0);
    if nargin<6 || isempty(satp) % choose on scatter plot
        figure(h);
        set(h,'numbertitle','off','name','Select ROI/vertices');
        if i==1
            msgbox({'Select the ROI/vertices, and the center, for determining hue and saturation.',...
                '<Esc> for default: limits of histogram for vertices, mean of vertices for center.',...
                'When done: ''A'' to add a vertex, right click to delete it, double-click to accept;'});
        end
        hroi = impoly; %permette di tracciare la roi poligonale
        wait(hroi); %permette di cambiare la roi a piacimento;
        %'A' per aggiugnere un vertice, doppio click sulla ROI per accettare.
        if ~isempty(hroi)
            satp = getPosition(hroi);
        else
            satp=[];
        end
        if size(satp,1)<2||all(isnan(satp(:)))||size(satp,2)~=2
            satp=[[res.r{1}(1),res.r{2}(1)];[res.r{1}(end),res.r{2}(end)];...
                [res.r{1}(1),res.r{2}(end)];[res.r{1}(end),res.r{2}(1)]];
        end
        set(h,'numbertitle','off','name','Select center');
        hpoint=impoint;
        wait(hpoint);
        if isempty(hpoint)
            cntr=[];
        else
            cntr=getPosition(hpoint);
        end
    end
    val=res.rmat(:)/max(max(res.rmat(2:end-1,2:end-1)));
    val(val>1)=1;
    val(val<0)=0;
    if logI>0 && logI<1
        val(val>logI)=logI+(log10(logI)-log10(val(val>logI)))*(1-logI)/log10(logI);
    end
    [hs,rgb,satp_ord,cntr]=phasorHSV(satp,res.r,cntr,spar,spar2,val,doscatter);
    rgb=reshape(rgb,length(res.r{1}),length(res.r{2}),3);
    if doplot
        h2=figure('Color','k','Position',[200,200,1350,500],...
            'PaperUnits','points','PaperSize', [1350 500],...
            'PaperPositionMode', 'auto','InvertHardcopy','off');
        subplot(1,5,[3,4,5])
        hold on;
        imagesc([res.r{1}(1),res.r{1}(end)],[res.r{2}(1),res.r{2}(end)],permute(rgb,[2 1 3]));
        axis xy;
        plot(cos(pi*t)/2+1/2,sin(pi*t)/2,'w','Linewidth',1.5);
        plot(t,zeros(size(t)),'w','Linewidth',1.5);
        if drawsatp
            switch abs(drawsatp)
                case 2
                    satp_ord(end+1,:)=satp_ord(1,:); %#ok<AGROW>
                    plot(satp_ord(:,1),satp_ord(:,2),'Color',[0.75,0.75,0.75],'Linewidth',1.5);
                case 3
                    satp_ord(end+1,:)=satp_ord(1,:); %#ok<AGROW>
                    plot(satp_ord(:,1),satp_ord(:,2),'Color',[0.75,0.75,0.75],'Linewidth',1.5);
                    scatter(satp_ord(:,1),satp_ord(:,2),80,[0.75,0.75,0.75],'+','Linewidth',1.5);
                otherwise
                    scatter(satp_ord(:,1),satp_ord(:,2),80,[0.75,0.75,0.75],'+','Linewidth',1.5);
            end
            if drawsatp<0
                scatter(cntr(1),cntr(2),80,'c','x','Linewidth',1.5);
%                 scatter(0.35823,0.34331,70,[0.75,0.75,0.75],'+','Linewidth',1.5); %temp
            end
        end
        axis equal;
        ax = gca;
        ax.Color = 'k';
        ax.XLim=[-0.09,1.09];
        ax.YLim=[-0.05,0.59];
        ax.XAxisLocation='origin';
        ax.YAxisLocation='origin';
        ax.XColor='w';
        ax.YColor='w';
        ax.Layer = 'top';
        ax.FontSize = 20;
        ax.XLabel.String='G';
        ax.XLabel.Position=[0.53,-0.13];
        ax.YLabel.String='S';
        ax.YLabel.Position=[-0.15,0.335];
        ax.LineWidth=1.5;
    end
    val=chos_I(:)/max(chos_I(:));
    if logI>0 && logI<1
        val(val>logI)=logI+(log10(logI)-log10(val(val>logI)))*(1-logI)/log10(logI);
    end
    [hsim,rgbim]=phasorHSV(satp,[chos_g(:),chos_s(:)],cntr,spar,spar2,val,doscatter);
    rgbim=reshape(rgbim,m,n,3);
    if doplot
%         figure;
%         imagesc(chos_I);
%         axis image;
        figure(h2);
        axim=subplot(1,5,[1,2]);
        imagesc(rgbim);
        axis image;
        axim.YTickLabel={};
        axim.XTickLabel={};
        axim.TickLength=[0,0];
        axim.LineWidth=1;
        axim.YColor=[1,1,1];
        axim.XColor=[1,1,1];
        axim.LineWidth=1.5;
        if doplot>1
            if isempty(Filename{i})
                Filename{i}=['fig' num2char(i)];
            end
            [fdir, fileroot, ~]=fileparts(Filename{i});
            if drawsatp
                fileroot=[fileroot, '_d_']; %#ok<AGROW>
            end
            fileroot=[fileroot, suff];
            fileroot=[fileroot 'g' num2str(sgauss) 'l' num2str(logI)]; %#ok<AGROW>
            saveas(h2,fullfile(Pathname,fdir,[fileroot '.pdf']),'pdf'); %implement try catch or check for existing/open file
            print(fullfile(Pathname,fdir,[fileroot '.png']),'-dpng');
        end
    end
    
    if length(c)>1
        h2out{i}=h2;
        hsout{i}=hs;
        rgbout{i}=rgb;
        hsimout{i}=hsim;
        rgbimout{i}=rgbim;
        Ithvecout{i}=Ithvec;
        limxyout{i}=limxy;
    end
end
if length(c)>1
    h2=h2out;
    hs=hsout;
    rgb=rgbout;
    hsim=hsimout;
    rgbim=rgbimout;
    Ithvec=Ithvecout;
    limxy=limxyout;
end

% if doplot
%     h2=figure; hold on;
%     plot(cos(pi*t)/2+1/2,sin(pi*t)/2,'k');
%     plot(t,zeros(size(t)),'k');
%     for i=1:length(c)
%         g(i)=scatter(av_gfw(i),av_sfw(i),'+','MarkerEdgeColor',cgr(i,:));
%         %scatter(av_gf(i),av_sf(i),'x','MarkerEdgeColor',cgr(i,:))
%         plot(mingf(i):(Dgf(i)/50):maxgf(i),par(i,1)+par(i,2)*(mingf(i):(Dgf(i)/50):maxgf(i)),'-','Color',cgr(i,:))
%         %plot(mingf(i):(Dgf(i)/50):maxgf(i),parw(i,1)+parw(i,2)*(mingf(i):(Dgf(i)/50):maxgf(i)),':','Color',cgr(i,:))
%     end
%     legend(g,Filename,'Location','northoutside','Interpreter','none','FontSize',8);
%     legend('boxoff');
% end
