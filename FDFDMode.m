function mode=FDFDMode(GuideCross,k0,kfind,ModeNum)
% GuideCross: The structures conatin the material and geometry informations
% kfind: 

% Fetch the geometry and material informations
nx_square=GuideCross.nx_square;
ny_square=GuideCross.ny_square;
nz_square=GuideCross.nz_square;
NX=GuideCross.NX;NY=GuideCross.NY;
NXY=NX*NY;
dx=GuideCross.dx;
dy=GuideCross.dy;

% %get all the coeffiecients and form a sparse matrix
% global n_clad   nx_guide    ny_guide    nz_guide    nx_square ny_square nz_square k0 lamda;
% global NX NY dx dy AAA;
%the following codes are devoted to obtain the coefficients of the coupled
%equations.we have 14 coeffiecients with different length.
SS=zeros(NXY-2*NX,1); WW1=zeros(NXY-2*NX,1);EE1=zeros(NXY-2*NX,1);NN=zeros(NXY-2*NX,1);
CSS1=zeros(NXY-2*NX,1); CNN1=zeros(NXY-2*NX,1); CWW1=zeros(NXY-2*NX,1); CEE1=zeros(NXY-2*NX,1);CCC1=zeros(NXY-2*NX,1);
CSS2=zeros(NXY-2*NX,1); CNN2=zeros(NXY-2*NX,1); CWW2=zeros(NXY-2*NX,1); CEE2=zeros(NXY-2*NX,1);CCC2=zeros(NXY-2*NX,1);
in_epsx=1./nx_square;   in_epsy=1./ny_square;   in_epsz=1./nz_square;
clear nz_square nx_square ny_square;
net_ratio=dx/dy;

%***************************************************************
q=0;
for r=2:NY-1
    for s=1:NX
            q=q+1;
            if s==1||s==NX
                WW1(q)=0;
                 EE1(q)=0;
                 NN(1)=0;
                 SS(q)=0;
                 CNN1(q)=0;
                 CNN2(q)=0;
                 CSS1(q)=0;
                 CSS2(q)=0;
                 CEE1(q)=0;
                 CEE2(q)=0;
                 CWW1(q)=0;
                 CWW2(q)=0;
                 CCC1(q)=0;
                 CCC2(q)=0;
            else
                
                WW1(q)=0.5*net_ratio*((in_epsz(r-1,s)-in_epsz(r,s))/(in_epsy(r-1,s)+in_epsy(r,s))-(in_epsy(r,s-1)-in_epsy(r-1,s-1))/(in_epsy(r,s-1)+in_epsy(r-1,s-1)));
                 EE1(q)=0.5*net_ratio*((in_epsz(r,s-1)-in_epsz(r-1,s-1))/(in_epsy(r,s-1)+in_epsy(r-1,s-1))-(in_epsy(r-1,s)-in_epsy(r,s))/(in_epsy(r,s)+in_epsy(r-1,s)));
                 NN(q)=0.5*net_ratio*((in_epsz(r,s-1)-in_epsz(r,s))/(in_epsx(r,s-1)+in_epsx(r,s))-(in_epsx(r-1,s)-in_epsx(r-1,s-1))/(in_epsx(r-1,s-1)+in_epsx(r-1,s)));
                 SS(q)=0.5*net_ratio*((in_epsz(r,s)-in_epsz(r,s-1))/(in_epsx(r,s-1)+in_epsx(r,s))-(in_epsx(r,s-1)-in_epsx(r,s))/(in_epsx(r,s-1)+in_epsx(r,s)));
           
                CSS1(q)=net_ratio*(in_epsz(r,s-1)/(in_epsy(r-1,s-1)+in_epsy(r,s-1))+in_epsz(r,s)/(in_epsy(r-1,s)+in_epsy(r,s)));
                CSS2(q)=1;  
                CNN1(q)=net_ratio*(in_epsz(r-1,s-1)/(in_epsy(r-1,s-1)+in_epsy(r,s-1))+in_epsz(r-1,s)/(in_epsy(r-1,s)+in_epsy(r,s)));
                CNN2(q)=1; 
                CWW2(q)=net_ratio*(in_epsz(r-1,s-1)/(in_epsx(r-1,s-1)+in_epsx(r-1,s))+in_epsz(r,s-1)/(in_epsx(r,s)+in_epsx(r,s-1)));
                CWW1(q)=1;   
                CEE2(q)=net_ratio*(in_epsz(r-1,s)/(in_epsx(r-1,s-1)+in_epsx(r-1,s))+in_epsz(r,s)/(in_epsx(r,s)+in_epsx(r,s-1)));
                CEE1(q)=1; 
                CCC1(q)=k0^2*(dx^2)*(1/(in_epsy(r,s)+in_epsy(r-1,s))+1/(in_epsy(r-1,s-1)+in_epsy(r,s-1)))-2-(dx^2/dy^2)*((in_epsz(r-1,s)+in_epsz(r,s))/(in_epsy(r,s)+in_epsy(r-1,s))+(in_epsz(r-1,s-1)+in_epsz(r,s-1))/(in_epsy(r-1,s-1)+in_epsy(r,s-1)));
                CCC2(q)=k0^2*(dy^2)*(1/(in_epsx(r-1,s-1)+in_epsx(r-1,s))+1/(in_epsx(r,s-1)+in_epsx(r,s)))-2-(dy^2/dx^2)*((in_epsz(r-1,s-1)+in_epsz(r-1,s))/(in_epsx(r-1,s-1)+in_epsx(r-1,s))+(in_epsz(r,s-1)+in_epsz(r,s))/(in_epsx(r,s-1)+in_epsx(r,s)));
            end
    end
end%pass
%WW1;
%************************************************************************
 %*************************************************************************

 %*************************************************************************
 %*****************End of Coefficients definition**************************
 %*************************************************************************
%to form a sparse matrix with a band width of 14.

DimS=2*NXY;
SS=SS(2:NXY-2*NX-1,1);
NN=NN(2:NXY-2*NX-1,1);
WW1=WW1(2:NXY-2*NX-1,1);
EE1=EE1(2:NXY-2*NX-1,1);
CNN1=CNN1(2:NXY-2*NX-1,1);
CNN2=CNN2(2:NXY-2*NX-1,1);
CSS1=CSS1(2:NXY-2*NX-1,1);
CSS2=CSS2(2:NXY-2*NX-1,1);
CWW1=CWW1(2:NXY-2*NX-1,1);
CWW2=CWW2(2:NXY-2*NX-1,1);
CEE1=CEE1(2:NXY-2*NX-1,1);
CEE2=CEE2(2:NXY-2*NX-1,1);
CCC1=CCC1(2:NXY-2*NX-1,1);
CCC2=CCC2(2:NXY-2*NX-1,1);
%****************************************************************************

ss=sparse(NXY+NX+2:2*NXY-NX-1,2+2*NX:NXY-1,SS,DimS,DimS); 
nn=sparse(NXY+NX+2:2*NXY-NX-1,2:NXY-2*NX-1,NN,DimS,DimS);

ww1=sparse(NX+2:NXY-NX-1,NXY+NX+1:2*NXY-NX-2,WW1,DimS,DimS);
ee1=sparse(NX+2:NXY-NX-1,NXY+NX+3:2*NXY-NX,EE1,DimS,DimS);

css1=sparse(NX+2:NXY-NX-1,2*NX+2:NXY-1,CSS1,DimS,DimS);
css2=sparse(NXY+NX+2:2*NXY-NX-1,NXY+2*NX+2:2*NXY-1,CSS2,DimS,DimS);

cww1=sparse(NX+2:NXY-NX-1,NX+1:NXY-NX-2,CWW1,DimS,DimS);
cww2=sparse(NXY+NX+2:2*NXY-NX-1,NXY+NX+1:2*NXY-NX-2,CWW2,DimS,DimS);

cee1=sparse(NX+2:NXY-NX-1,NX+3:NXY-NX,CEE1,DimS,DimS);
cee2=sparse(NXY+NX+2:2*NXY-NX-1,NXY+NX+3:2*NXY-NX,CEE2,DimS,DimS);

cnn1=sparse(NX+2:NXY-NX-1,2:NXY-2*NX-1,CNN1,DimS,DimS);
cnn2=sparse(NXY+NX+2:2*NXY-NX-1,NXY+2:2*NXY-2*NX-1,CNN2,DimS,DimS);

ccc1=sparse(NX+2:NXY-NX-1,NX+2:NXY-NX-1,CCC1,DimS,DimS);
ccc2=sparse(NXY+NX+2:2*NXY-NX-1,NXY+NX+2:2*NXY-NX-1,CCC2,DimS,DimS);

CoeMat=css1+cww1+ccc1+cee1+cnn1+ee1+ww1+nn+ss+css2+cww2+cee2+cnn2+ccc2;

clear ss  ww1 ee1 nn css2 cww2 ccc1 ccc2  cee2 cnn2 css1 cww1 cee1 cnn1;
clear  WW1 EE1 NN CSS1 CWW1   CNN1 CSS2 CWW2   CNN2 CEE1 CEE2 SS CCC1 CCC2;
clear in_epsx in_epsy in_epsz DimS;

% Eigenvalues solutions
options.tol = 1e-7;
options.disp = 0;
[hfield_1d,D]=eigs(CoeMat,speye(size(CoeMat)),ModeNum,kfind^2*dx*dy,options);
clear AAA;
%display field distribution and propagation constant
neff=sqrt(D/k0^2/dx/dy);

Hx=cell(ModeNum);Hy=cell(ModeNum);

for l=1:ModeNum
Hx{l}=transpose(reshape(hfield_1d(1:NXY,l),NX,NY));
Hy{l}=transpose(reshape(hfield_1d(NXY+1:2*NXY,l),NX,NY));
end
mode.Hx=Hx;
mode.Hy=Hy;
mode.neff=neff;

end