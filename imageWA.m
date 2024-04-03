function imout=imageWA(imin,r,imw)
if nargin<1 || isempty(imin)
    imout=[];
    return;
end
if nargin<2 || isempty(r)
    r=0.5;
end
if nargin<3 || isempty(imw)
    imout=imgaussfilt(imin,r);
    return;
end
szim=size(imin);
szw=size(imw);
if any(szim>szw)
    imw(end+1:szim(1),:)=0;
    imw(:,end+1:szim(2))=0;
end
imw(isnan(imw)|isnan(imin))=0;
imin(isnan(imin))=0;

if any(szw>szim) % may be useless
    imin(end+1:szw(1),:)=0;
    imin(:,end+1:szw(2))=0;
end
szim=size(imin);

imout=NaN(szim);
n=max(floor(3*abs(r)),1);
[Rv,Cv]=meshgrid(-n:n);
gw=normpdf(sqrt(Rv.^2+Cv.^2),0,r);
for rw=1:szim(1)
    for cl=1:szim(2)
        w=zeros(2*n+1);
        rww=rw+Rv;
        clw=cl+Cv;
        isin=rww>0 & rww<=szim(1) & clw>0 & clw<=szim(2);
        indim=rww(isin)+(clw(isin)-1)*szim(1);
%         indim=sub2ind(szim,rww(isin),clw(isin));
        w(isin)=gw(isin).*imw(indim);
        imout(rw,cl)=sum(imin(indim).*w(isin))/sum(w(isin));
%         if isnan(imout(rw,cl))
%             keyboard;
%         end
%         imout(rw,cl)=wmean(imin(indim),w(isin));
    end
end
        
    