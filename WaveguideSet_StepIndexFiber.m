n_clad=1.45;
nx_guide=1.46;   ny_guide=1.46;  nz_guide=1.46;


radius=500e-9; % radius of fiber
w=8*radius; % Simulating range

% Geometry and grid
core_number=50; % Grid numbers
dx=radius/core_number;
dy=dx;
NX=floor((w)/dx);
NY=NX;

nx_square=zeros(NY,NX);
nz_square=zeros(NY,NX);
ny_square=zeros(NY,NX);

for r=1:NY
    for s=1:NX
       if (dx*s-4*radius)^2+(dy*r-4*radius)^2<=radius^2
            nx_square(r,s)=nx_guide^2;
            ny_square(r,s)=ny_guide^2;
            nz_square(r,s)=nz_guide^2;
       else
            nx_square(r,s)=n_clad^2;
            ny_square(r,s)=n_clad^2;
            nz_square(r,s)=n_clad^2;
            
            
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