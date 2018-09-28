function visiblePtInds=HPR(p,C,param)

%visiblePtInds=HPR(p,[-2,0,0],exp(0.3));    
dim=size(p,2);
numPts=size(p,1);
% Move the points s.t. C is the origin
p=p-repmat(C,[numPts 1]);
% Calculate ||p||
normp=sqrt(dot(p,p,2));
% Sphere radius
R=repmat(max(normp)*(10^param),[numPts 1]);
%Spherical flipping
P=p+2*repmat(R-normp,[1 dim]).*p./repmat(normp,[1 dim]);
%convex hull
visiblePtInds=unique(convhulln([P;zeros(1,dim)]));
visiblePtInds(visiblePtInds==numPts+1)=[];

end

