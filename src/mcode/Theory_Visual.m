%% Initial parameters
% Constants %
uo=4*pi*10^-7; % Diamagnetic vac const.H/mm
I=8; % Max current A
% Magnet field %
tic
% Estimated values %
a=99.5*10^-3; % Solenoid heigth 
b=41.8*10^-3; % Solenoid width 
ri=30*10^-3; % Solenoid inner radius
Ne=1000; % Number of windings
rm=ri+a/2;
ce=sqrt(-(b^2-a^2)/12);
re=rm*(1+(a^2)/(24*rm^2));
em=9.1*10^-31; % Electron mass
eq=1.6*10^-19; % Electron mass
vc=3*10^8;
ro=0.0272;
co=0;
No=541;
pz=sqrt(2*(3.5*10^6*eq)*em); % Given electron impulse
%
Bz = @(z, r,c, N) ((uo*N*I*((((r+i*c).^2)./(((z.^2)+(r+i*c).^2).^1.5))+(((r-i*c).^2)./(((z.^2)+(r-i*c).^2).^1.5)))/4));
Bzs =@(z, l, N, r) ((uo*N*I*((((z+l./2))./(sqrt((r.^2)+(z+l./2).^2)))-(((z-l./2))./(sqrt((r.^2)+(z-l./2).^2))))/2));
d1Bz=@(z,r,c,N) N.*(z.*1.0./((c.*1i-r).^2+z.^2).^(5.0./2.0).*(c.*1i-r).^2.*3.0+z.*1.0./(z.^2+(c.*1i+r).^2).^(5.0./2.0).*(c.*1i+r).^2.*3.0).*(-2.513274122871835e-6);
d2Bz=@(z,r,c,N)N.*(1.0./((c.*1i-r).^2+z.^2).^(5.0./2.0).*(c.*1i-r).^2.*3.0+1.0./(z.^2+(c.*1i+r).^2).^(5.0./2.0).*(c.*1i+r).^2.*3.0-z.^2.*1.0./((c.*1i-r).^2+z.^2).^(7.0./2.0).*(c.*1i-r).^2.*1.5e+1-z.^2.*1.0./(z.^2+(c.*1i+r).^2).^(7.0./2.0).*(c.*1i+r).^2.*1.5e+1).*(-2.513274122871835e-6);
F3= @(r,c, N) -integral(@(z) Bz3(z,r,c, N), -inf, inf)./2;
F4= @(r,c, N) integral(@(z) Bz4(z,r,c, N), -inf, inf);
F2= @(r,c, N) 2*integral(@(z) Bz2(z,r,c, N), 0, inf);
% 
xsol = linspace(-0.45,0.45,1024);
xsol = xsol(:);
figure(1)
plot( xsol, Bz(xsol, re, ce, Ne), xsol, Bzs(xsol, b, Ne, ri));
xlabel('z (m)');
ylabel('B_z on axis (T)');
legend('Multilayer','Single layer');
%
zsol= -0.4:0.0001:0.4;
xrsq= linspace(0000001,0.008,10^3);
xrsq=xrsq(:);
yc=linspace(0.000001,0.04,10^3);
yc=yc(:);
F3plot=@(r, c) abs(F3( r, c, No));
F3plot2=@(r, c) abs(F3(r,c, 1000));
F3plot05=@(r, c) abs(F3(r,c, 250));
figure( 3)
fmesh(F3plot, [0.02 0.04 0.000001 0.05])
hold on
fmesh(F3plot2, [0.02 0.04 0.000001 0.05])
fmesh(F3plot05, [0.02 0.04  0.000001 0.05])
xlabel('Rsq');
ylabel('c');
zlabel('F3'); 
hold off
figure(4)
plot(zsol, Bz(zsol, ro, co, No) )
figure( 5)
plot(zsol, d1Bz(zsol,  ro, co, No),zsol,d2Bz(zsol,  ro, co, No) )
legend('First der.','Second der.');