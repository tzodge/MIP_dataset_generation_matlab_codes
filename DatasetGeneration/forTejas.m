clear;clc;
close all
warning('off');
shapename='bunny';%name of the stl file
[tF,tP]=stlread(strcat(shapename,'.stl'));
tP=tP/(max(max(tP)-min(tP)))*1;
tP=tP-mean(tP);
partial_model_switch = 1
M=unique(tP,'rows');
cam_position=[-2,0,0];

if partial_model_switch == 1
    indices_partial_model = HPR(M, cam_position, exp(0.5));
    partial_model = M(indices_partial_model,:);
    ind_not_in_partial_model = setdiff(1:size(M,1),indices_partial_model);
end


gg=M;
% vertices are repeated in stl file. So I need to select only unique ones.
tF1=tF;
for ii=1:length(tF)
    for jj=1:3
        tF1(ii,jj)=find(gg(:,1)==tP(tF(ii,jj),1) & gg(:,2)==tP(tF(ii,jj),2) & gg(:,3)==tP(tF(ii,jj),3));
    end
end

if partial_model_switch == 1
    faces_to_delete= zeros(size(tF1,1),1) ;
    j = 1;
    for i =1:size(tF1,1)
        if (ismember(tF1(i,1),ind_not_in_partial_model)||ismember(tF1(i,2),ind_not_in_partial_model)||ismember(tF1(i,3),ind_not_in_partial_model))
            faces_to_delete(j,1) = i;
            j = j+1;
        end    
    end
    
    faces_to_delete(all(~faces_to_delete,2),:) = [];
    tF1(faces_to_delete,:) = [];
    [Lia,map_full_partial_M_pts] = ismember( M , partial_model,'rows');
    
    for i =1:size(tF1,1)
        for j = 1:3 
            tF1(i,j) = map_full_partial_M_pts(tF1(i,j));
        end 
    end
    M = partial_model;
end


Nf=length(tF1);
Nm=length(M);

%%

F=zeros(Nf,Nm);

for ii=1:Nf
    F(ii,tF1(ii,:))=1;
end

%%
% I sampled 10 points from each face to form Msampled
pointPerFace=10;
Nmsampled=Nf*pointPerFace;
Msampled=zeros(Nmsampled,3);
tt=1;

for jj=1:Nf
    for ii=1:pointPerFace
        w=rand(1,3);
        w=w/sum(w);
        Msampled(tt,:)=w(1)*M(tF1(jj,1),:)+w(2)*M(tF1(jj,2),:)+w(3)*M(tF1(jj,3),:);
        tt=tt+1;
    end
end

%%
d2r=pi/180;
Xgt=[0.5*(2*rand(1,3)-1),180*d2r*(2*rand(1,3)-1)];
Xgt(4:6)=rotm2eul(eul2rotm(Xgt(4:6)));

gt=eye(4);
Rgt=eul2rotm(Xgt(4:6));
tgt=Xgt(1:3)';
gt(1:3,1:3)=Rgt;
gt(1:3,4)=tgt;

invTgt=inv(gt);
%%

Ns=6000;

outliers=0;%IF you have outliers, turn this to 1
outlierPercentage=10;
noise= 0;% If you have noise, turn this to 1
noiseSigma=5*1e-2;% also try 1e-1 for more noise and 1e-2 for less noise

partial=0;
partialPercentage=40;% also try 80 and 30

if partial==1
    Ns=floor(Ns*100/40);
end

o=zeros(Ns,1);
S=zeros(Ns,3);
C=zeros(Ns,Nm);
Cb=zeros(Ns,Nmsampled);
f1=zeros(Ns,Nf);
SigmaS=zeros(Ns,6);
B=zeros(Ns,6);
t=1;
wb=waitbar(0);
for ii=1:ceil(Ns/Nf)*Nf
    
    waitbar(ii/(ceil(Ns/Nf)*Nf),wb);
    idx=randi(Nf);
    w=rand(1,3);
    w=w/sum(w);
    C(t,tF1(idx,:))=w;
    f1(t,idx)=1;
    
    S(t,:)= w(1)*M(tF1(idx,1),:)+w(2)*M(tF1(idx,2),:)+w(3)*M(tF1(idx,3),:);
        idx1=knnsearch(Msampled((idx-1)*pointPerFace+1:idx*pointPerFace,:),S(t,:));
    Cb(t,(idx-1)*pointPerFace+idx1)=1;

    CovS=eye(3);
    SigmaS(t,:)=[CovS(1,1:3),CovS(2,2:3),CovS(3,3)];
    B(t,:)=SigmaS(t,:);
    % if the variable outliers in 1 then I add outliers. 1% of points
    if outliers==1 && rand>1-(outlierPercentage/100)
        S(t,:)=S(t,:)+(2*rand(1,3)-1)*0.5;
        C(t,:)=0;
        o(t,:)=1;
        f1(t,:)=0;
        
        Cb(t,:)=0;
        
    end
    Rn=eye(3);
    if noise==1
        nvec=cross((tP(tF(idx,2),:)-tP(tF(idx,1),:)),(tP(tF(idx,3),:)-tP(tF(idx,1),:)));
        nvec=nvec/norm(nvec);
        tmp2=cross(rand(1,3),nvec);tmp2=tmp2/norm(tmp2);
        tmp1=cross(tmp2,nvec);
        Rn(:,1)=tmp1';
        Rn(:,2)=tmp2';
        Rn(:,3)=nvec';
        S(t,:)=S(t,:)+(Rn*[normrnd(0,1e-6),normrnd(0,1e-6),normrnd(0,noiseSigma)]')';
        
        CovS=invTgt(1:3,1:3)*Rn*diag([1e-6,1e-6,noiseSigma])*Rn'*invTgt(1:3,1:3)';
        cholInvCov=chol(inv(CovS));
        B(t,:)=[cholInvCov(1,1:3),cholInvCov(2,2:3),cholInvCov(3,3)];
        SigmaS(t,:)=[CovS(1,1:3),CovS(2,2:3),CovS(3,3)];
    end
    
    
    
    
    
    
    
    
    S(t,:)=transpose(invTgt(1:3,1:3)*S(t,:)'+invTgt(1:3,4));
    
    
    
    
    if t==Ns
        break;
    end
    
    t=t+1;
end

if partial==1
    
    idx1=randi(length(S));
    idx2=knnsearch(S,S(idx1,:),'k',floor(partialPercentage/100*Ns));
    C=C(idx2,:);
    f1=f1(idx2,:);
    
    Cb=Cb(idx2,:);
    o=o(idx2,:);

    S=S(idx2,:);
        B=B(idx2,:);
    SigmaS=SigmaS(idx2,:);

    idx3=randperm(length(S));
    S=S(idx3,:);
    B=B(idx3,:);
    SigmaS=SigmaS(idx3,:);
    Ns=length(S);
end

for ii=1:Ns
St=transpose(gt(1:3,1:3)*S'+gt(1:3,4));
phi=0;
    phi=phi+sum(abs([B(ii,1),B(ii,2),B(ii,3);B(ii,2),B(ii,4),B(ii,5);B(ii,3),B(ii,5),B(ii,6)]*(gt(1:3,1:3)'*M'*C(ii,:)'-gt(1:3,1:3)'*gt(1:3,4)-S(ii,:)')));
    
end
phi=phi/Ns;
%%
% trisurf(tF,tP(:,1),tP(:,2),tP(:,3),'FaceAlpha',0.8); hold on; axis equal
% scatter3(Msampled(:,1),Msampled(:,2),Msampled(:,3),'b','fill');
trisurf(tF1,M(:,1),M(:,2),M(:,3),'FaceAlpha',0.8); hold on; axis equal
scatter3(St(:,1),St(:,2),St(:,3),'.m');
scatter3(St(1:20,1),St(1:20,2),St(1:20,3),'k','fill');
scatter3(S(:,1),S(:,2),S(:,3),'.b');
scatter3(cam_position(1),cam_position(2),cam_position(3), '.r')


%%
% 
filename="D:\CMU Thesis july 2018\MIP\datasets\datasets";
dlmwrite(strcat(filename,'\',shapename,"\F.txt"),F);
dlmwrite(strcat(filename,'\',shapename,'\M.txt'),M);
dlmwrite(strcat(filename,'\',shapename,'\S.txt'),S);
dlmwrite(strcat(filename,'\',shapename,'\St.txt'),St);
dlmwrite(strcat(filename,'\',shapename,'\gt.txt'),gt);
dlmwrite(strcat(filename,'\',shapename,'\Msampled.txt'),Msampled);
dlmwrite(strcat(filename,'\',shapename,'\SigmaS.txt'),SigmaS);
dlmwrite(strcat(filename,'\',shapename,'\B.txt'),B);
% dlmwrite(strcat(filename,'\',shapename,'\C.txt'),C);
% dlmwrite(strcat(filename,'\',shapename,'\Cb.txt'),Cb);
% dlmwrite(strcat(filename,'\',shapename,"\f.txt"),f1);
dlmwrite(strcat(filename,'\',shapename,'\o.txt'),o);
