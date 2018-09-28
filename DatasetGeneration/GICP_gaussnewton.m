function [Rk,tk] = GICP_gaussnewton(model,source,SigmaSi)

% model= Nx3
% source=Nx3
% SigmaSi=Nx6


Rk=eye(3);tk=zeros(3,1);
n=size(model,1);

iterNum=100;

for kk=1:iterNum
    
    
    f0=zeros(3*n,1);
    for ii=1:n
        f0(3*ii-2:3*ii)=model(ii,:)'-Rk*source(ii,:)'-tk;
    end
    
    J=zeros(3*n,6);
    
    for ii=1:n
        Rkp2n=Rk*source(ii,:)';
        skewRkp2n=[0,-Rkp2n(3),Rkp2n(2);...
            Rkp2n(3),0,-Rkp2n(1);...
            -Rkp2n(2),Rkp2n(1),0];
        J(3*ii-2:3*ii,1:3)=skewRkp2n;
        J(3*ii-2:3*ii,4:6)=-eye(3);
        
    end
    
    Sigma=zeros(3*n,3*n);
    
    for ii=1:n
        SigmaSiMat=[SigmaSi(ii,1),SigmaSi(ii,2),SigmaSi(ii,3);...
            SigmaSi(ii,2),SigmaSi(ii,4),SigmaSi(ii,5);...
            SigmaSi(ii,3),SigmaSi(ii,5),SigmaSi(ii,6)];
        
        Sigma(3*ii-2:3*ii,3*ii-2:3*ii)=Rk*SigmaSiMat*Rk';
    end
    
    
    
    deltaP=-(J'*(Sigma\J))\(J'*(Sigma\f0));
    deltaA=deltaP(1:3);
    deltaT=deltaP(4:6);
    
    
    theta =norm(deltaA);
    v=deltaA'/norm(deltaA);
    cth = cos(theta);
    sth = sin(theta);
    vth = (1 - cth);
    
    vx = v(1);
    vy = v(2);
    vz = v(3);
    
    Rtmp =[ vx.*vx.*vth+cth,     vy.*vx.*vth-vz.*sth, vz.*vx.*vth+vy.*sth; ...
        vx.*vy.*vth+vz.*sth, vy.*vy.*vth+cth,     vz.*vy.*vth-vx.*sth; ...
        vx.*vz.*vth-vy.*sth, vy.*vz.*vth+vx.*sth, vz.*vz.*vth+cth];
    
    
    Rk=Rtmp*Rk;
    
    tk=tk+deltaT;
    
    if norm(deltaA)<1e-4*pi/180 && norm(deltaT)<1e-4
        break;
    end
    
end



end

