%% Magnetic Field Simulation %%
% Notes %
% Algorithm for optimization of Solenoid dimensions for focusing of
% electron beam. 
% Given restrictions/initial parameters:
% -focus length
% -effective field width
% -electron energy
% -approximate solenoid dimensions/manufacturing restrictions
% -Max. operating value of current
% Units of all retrieved values are in CI units, unless stated otherwise
%% Initial parameters
% Constants %
uo=4*pi*10^-7; % Diamagnetic vac const.H/mm
I=8; % Max current A
% Magnet field %
tic
% Estimated values  %
    a=99.5*10^-3; % Solenoid heigth 
    b=41.8*10^-3; % Solenoid width
    ri=30*10^-3; % Solenoid inner radius
    Ne=400; % Number of windings
    rm=ri+a/2;
    ce=sqrt(-(b^2-a^2)/12);
    re=rm*(1+(a^2)/(24*rm^2));
    em=9.1*10^-31; % Electron mass
    eq=1.6*10^-19; % Electron mass
    vc=3*10^8; % light speed
    pz=sqrt(2*(3.5*10^6*eq)*em); % Given electron impulse
    rmc=ri+b/2;
    rec=rmc*(1+(b^2)/(24*rmc^2));
    rb= 0.001; % Estimated beam radiuse m
    riad=0.015; % Adjusted Solenoid inner radius m
%
Bz = @(z, r,c, N) ((uo*N*I*((((r+i*c).^2)./(((z.^2)+(r+i*c).^2).^1.5))+(((r-i*c).^2)./(((z.^2)+(r-i*c).^2).^1.5)))/4)); % Magnet field z komponent
Bz2= @(z,r,c, N) Bz(z,r,c,N).^2;
%
syms z
d1Bz=@(z,r,c,N) N.*(z.*1.0./((c.*1i-r).^2+z.^2).^(5.0./2.0).*(c.*1i-r).^2.*3.0+z.*1.0./(z.^2+(c.*1i+r).^2).^(5.0./2.0).*(c.*1i+r).^2.*3.0)...
    .*(-2.513274122871835e-6); % First z derivative of Bz
d2Bz=@(z,r,c,N)N.*(1.0./((c.*1i-r).^2+z.^2).^(5.0./2.0).*(c.*1i-r).^2.*3.0+1.0./(z.^2+(c.*1i+r).^2).^(5.0./2.0).*(c.*1i+r).^2.*3.0-z.^2.*1.0./((c.*1i-r).^2+z.^2).^(7.0./2.0).*(c.*1i-r).^2.*1.5e+1-z.^2.*1.0./(z.^2+(c.*1i+r).^2).^(7.0./2.0)...
    .*(c.*1i+r).^2.*1.5e+1).*(-2.513274122871835e-6); % Second z derivative of Bz
Bz3= @(z,r,c, N) Bz(z,r,c, N).*d2Bz(z,r,c, N);
Bz4= @(z,r,c, N) Bz(z,r,c, N).^4;
%
% 
F3= @(r,c, N) -integral(@(z) Bz3(z,r,c, N), -inf, inf)./2;
F4= @(r,c, N) integral(@(z) Bz4(z,r,c, N), -inf, inf);
F2= @(r,c, N) 2*integral(@(z) Bz2(z,r,c, N), 0, inf);
%% Optimization %%
    F3un= @(x) F3(x(1),x(2),x(3));
    F4un= @(x) F4(x(1),x(2),x(3));
    FS= @(r,c) F3(r,c,N)+0.01*F4(r,c,N);
    FSun= @(x) FS(x(1),x(2),x(3));
    fun = @(x) F3(x(1),x(2),x(3))+eq^2/(3*pz^2)*F4(x(1),x(2),x(3)); % Normalized  spherical aberration
% Options %
options = optimoptions(@fmincon,'Algorithm','interior-point',...
    'ConstraintTolerance',1.0e-15,...
    'MaxIterations',10000, ...
    'MaxFunctionEvaluations',10000,...
    'OptimalityTolerance',1.0e-15,  ...
    'FunctionTolerance', 1.0e-15, ...
    'StepTolerance',1.0e-17);
%
nonlcon = @const; % Constraints function
x0 = [re ce Ne];
A = []; % No other linear constraints
bb = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x = fmincon(fun,x0,A,bb,Aeq,beq,lb,ub,nonlcon,options) % Vector of optimized parameters with F3 and F4 taken into consideration
fun1 = @(y) F3(y(1),y(2),y(3));
y = fmincon(fun1,x0,A,bb,Aeq,beq,lb,ub,nonlcon,options) % Vector of optimized parameters with F3 taken into consideration
ys = fmincon(fun1,x0,A,bb,Aeq,beq,lb,ub,[],options)
%% Double check of calculated features %%
syms z ax
rey=@(ax) (riad+ax/2)*(1+(ax^2)/(24*(riad+ax/2)^2));
% Output the parameters and solenoid properties that follow from them
    MagnetLengthA=vpasolve(rey(ax)==x(1),ax)
    MagnetLengthB=sqrt(-12*x(2).^2+MagnetLengthA(2)^2)
    HalbBreite=vpasolve(Bz(z,x(1), x(2), x(3))==0.5*Bz(0,x(1), x(2), x(3)),z)
    maxBz=Bz(0, x(1), x(2), x(3))
    F2=  @(x1,x2, x3) 2*integral(@(z) Bz2(z,x1,x2, x3), 0, inf);
    f=   @(x1,x2, x3) 1/(F2(x1,x2, x3).*(eq/(2*pz))^2);
    fokusf=f( x(1), x(2), x(3))
    syms z
    HalbBreiteF3=vpasolve(Bz(z,y(1), y(2),y(3))==0.5*Bz(0,y(1), y(2),y(3)),z)
    HalbBreiteS=vpasolve(Bz(z,ys(1), ys(2), ys(3))==0.5*Bz(0,ys(1), ys(2), ys(3)),z)
    maxBzF3=Bz(0, y(1), y(2),y(3))
    F2=  @(x1,x2, x3) 2*integral(@(z) Bz2(z,x1,x2, x3), 0, inf);
    f=   @(x1,x2, x3) 1/(F2(x1,x2, x3).*(eq/(2*pz))^2);
    SphericalAbbMinF3sumF4 = eq^2*rb^4*F3(x(1),x(2),x(3))/(4*pz^2)+(eq^2/(3*pz^2))*F4(x(1),x(2),x(3))*(eq^2*rb^4)/(4*pz^2) % Searched spher. ab. 
    SphericalAbbMinF3 = eq^2*rb^4*F3(y(1),y(2),y(3))/(4*pz^2) % Alternative spher. ab. with considering only F3
    SphericalAbbREGAE = eq^2*rb^4*F3(rec,ce,Ne)/(4*pz^2)+eq^4*rb^4*F4(rec,ce,Ne)/(12*pz^4) % Spher. ab. for REGAE
    emittance= eq^2*(rb/2)^4*F3(x(1),x(2),x(3))/(3*sqrt(2)*vc*em*pz)
    emittanceREGAE=eq^2*(rb/2)^4*F3(rec, ce, 2.5*Ne)/(3*sqrt(2)*vc*em*pz)
    % focal Lengths
    fokusfF3=f( y(1), y(2),y(3)) 
    fokusfS=f( ys(1), ys(2), ys(3))
    fokusregae= f( rec, ce, 2.5*Ne)
%% Visualization %% 
% Here ctr is put before all parameters to note the method they were
% estimated with
ctra=26.823*10^-3;
ctrb=47.544*10^-3;
ctrrmc=28.438*10^-3;
ctrNe=4000.781/8;
ctrce=sqrt(-(ctrb^2-ctra^2)/12);
ctrrec=ctrrmc*(1+(ctra^2)/(24*ctrrmc^2));
figure(2);
plot( xsol, Bz(xsol, x(1), x(2), x(3)), xsol, Bz(xsol, rec, ce, 2.5*Ne),xsol, Bz(xsol, ctrrec, ctrce, ctrNe));
line([-HalbBreite;HalbBreite],[0.055;0.055]); % Plot effective width of the field
line([-0.0176/2;0.0176/2],[0.05;0.05],'color','red'); % Plot magnet width(length)
legend('Retreived field for 8 A', 'Field for magnet with REGAE parameters','Constrained trust region','Effective field length- 38mm','Solenoid width- 18mm');
xlabel('z (m)');
ylabel('B_z on axis (T)');
title('Solenoid field at 8 A');
%% Part for testing the derivatives of Bz, and graphic display of integrals
zsol= -0.4:0.0001:0.4;
zsold= -0.2:0.0001:0.2;
xrsq= linspace(0000001,0.008,10^3);
xrsq=xrsq(:);
yc=linspace(0.000001,0.04,10^3);
yc=yc(:);
F3plot=@(r, c) abs(F3(r,c, x(3)));
F3plot2=@(r, c) abs(F3(r,c, 1000));
F3plot05=@(r, c) abs(F3(r,c, 250));
figure( 3)
fmesh(F3plot, [0.02 0.09 -0.02 0.05])
hold on
fmesh(F3plot2, [0.02 0.09 -0.02 0.05])
fmesh(F3plot05, [0.02 0.09  -0.02 0.05])
line([re;re],[ce;ce],[0;1.5],'color','red'); % Plot REGAE param.
line([x(1);x(1)],[0.001;0.001],[0;1.5],'color','green'); % Plot param. estimated with Interior point method
line([ctrrec;ctrrec],[ctrce;ctrce],[0;1.5],'color','blue'); % Plot the param. estimated with alt. method
xlabel('Rsq');
ylabel('c');
zlabel('F3'); 
hold off
figure( 5)
plot(zsold, d1Bz(zsold, x(1), x(2), x(3)),zsold,d2Bz(zsold, x(1), x(2), x(3)) )
legend('First der.','Second der.');
xlabel('z (m)');
ylabel('dB_z/dz on axis (T/m)');
toc
%% Nonlinear costraints %%
function [c,ceq] = const(x)
% Same constants as above:
uo=4*pi*10^-7;
I=8;
a=99.5*10^-3;
b=41.8*10^-3;
ri=30*10^-3;
Ne=400;
rm=ri+a/2;
ce=sqrt(-(b^2-a^2)/12);
re=rm*(1+(a^2)/(24*rm^2));
em=9.1*10^-31;
eq=1.6*10^-19;
pz=sqrt(2*(3.5*10^6*eq)*em);
Bz = @(z,x1,x2,x3) ((uo.*x3.*I.*((((x1+i*x2).^2)./(((z.^2)+(x1+i*x2).^2).^1.5))+(((x1-i*x2).^2)./(((z.^2)+(x1-i*x2).^2).^1.5)))./4));
Bz2= @(z,x1,x2,x3) Bz(z,x1,x2,x3).^2;
F2=  @(x1,x2, x3) 2*integral(@(z) Bz2(z,x1,x2, x3), 0, inf);
f=   @(x1,x2, x3) 1/(F2(x1,x2, x3).*(eq/(2*pz))^2);
% The restrictions themselves:
c(2) = -f(x(1),x(2), x(3))+0.5; % Minimal focal length
c(3) = -x(2)+0.0001; % Positivity the of parameter c=sqrt(...)
c(4) = -x(1); % Positivity of the parameter a
c(1)= abs(-Bz(0.025,x(1),x(2), x(3))+0.5*Bz(0,x(1),x(2), x(3)))-0.01; % Required Bz half width
ceq(1) = -Bz(0,x(1),x(2), x(3))+0.11; % Minimal Bz(z=0)

end