% The main program for simulating the eigen mode analysis

% Set the material and geometry
% % the step index fiber
% WaveguideSet_StepIndexFiber;
% Rectangle waveguide
WaveguideSet_RectGuide;
% Wavelength
lamda=650e-9; % Wavelength
k0=2*pi/lamda;

kfind=k0*sqrt(max(GuideCross.nx_square(:)));
ModeNum=2;
%% Get the eigen field
mode=FDFDMode(GuideCross,k0,kfind,ModeNum);
%%  plot the results
x_list=dx:dx:dx*NX;
y_list=dy:dy:dy*NY;
[xgrid,ygrid]=meshgrid(x_list/1e-9,y_list/1e-9);

l=1;
figure
sgtitle(['Mode NO.',num2str(l),' neff=',num2str(mode.neff(l,l),3)]);
subplot(131)
pcolor(xgrid,ygrid,nx_square);
colormap jet; shading interp;
hold on
colorbar();title('refactive index');
xlabel('X (nm)');ylabel('Y (nm)');
subplot(132)
pcolor(xgrid,ygrid,real(mode.Hx{l}));
colormap jet; shading interp;
hold on
colorbar();title('Hx');
xlabel('X (nm)');ylabel('Y (nm)');
subplot(133)
pcolor(xgrid,ygrid,real(mode.Hy{l}));
colormap jet; shading interp;
colorbar();title('Hy');
xlabel('X (nm)');ylabel('Y (nm)');
set(gcf,'units','normalized', 'position',[0.2 0.2 0.8 0.4])