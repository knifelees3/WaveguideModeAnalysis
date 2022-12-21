n_clad=2.0;%isotropic
nn1_clad=1.73;
nn2_clad=1;
ns_clad=1.454;
nwe_clad=1;
nx_guide=2.01;   ny_guide=2.01;  nz_guide=2.01;

t=2.0e-7;
delta=00e-8;
w=t+delta;
ww=w+9*t;
wx1=(ww-w)/2;
wx2=wx1+w;
wy1=4*t;
wy2=wy1+0.5*t;
wy3=wy2+t;

core_number=100;
dx=t/core_number;
dy=dx;

NX=floor((ww)/dx);
NY=floor(10*t/dx);

nx_square=zeros(NY,NX);hx
nz_square=zeros(NY,NX);
ny_square=zeros(NY,NX);



for r=1:NY
    for s=1:NX
       if dx*s>=wx1&&dx*s<=wx2&&dy*r>=wy2&&dy*r<=wy3
            nx_square(r,s)=nx_guide^2;
            ny_square(r,s)=ny_guide^2;
            nz_square(r,s)=nz_guide^2;
       else
            if dx*s<wx1&&dy*r>=wy2&&dy*r<=wy3
            nx_square(r,s)=nwe_clad^2;
            ny_square(r,s)=nwe_clad^2;
            nz_square(r,s)=nwe_clad^2;
            end
            if dx*s>wx2&&dy*r>=wy2&&dy*r<=wy3
            nx_square(r,s)=nwe_clad^2;
            ny_square(r,s)=nwe_clad^2;
            nz_square(r,s)=nwe_clad^2;
            end
            
       end
             if dy*r<=wy1
            nx_square(r,s)=ns_clad^2;
            ny_square(r,s)=ns_clad^2;
            nz_square(r,s)=ns_clad^2;
             end
              if dy*r>wy1&&dy*r<wy2
            nx_square(r,s)=nn1_clad^2;
            ny_square(r,s)=nn1_clad^2;
            nz_square(r,s)=nn1_clad^2;
              end
              if dy*r>wy3
            nx_square(r,s)=nn2_clad^2;
            ny_square(r,s)=nn2_clad^2;
            nz_square(r,s)=nn2_clad^2;
              end
                  
    end
end
GuideCross.nx_square=nx_square;
GuideCross.ny_square=ny_square;
GuideCross.nz_square=nz_square;
GuideCross.NX=NX;
GuideCross.NY=NY;
NXY=NX*NY;
GuideCross.dx=dx;
GuideCross.dy=dy;