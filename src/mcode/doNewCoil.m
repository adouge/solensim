% copyright Dr. Thorsten Kamps, 2020
% All rights reserved.

% with beam pipe diameter 108 mm
% radius beam pipe = 108/2 mm = 54 mm
% coil can have min radius of 60 mm
% length of beam pipe = length of comp coil + 40 mm for vessel
% now with mu0 and ampereturns
xsol = linspace(-0.45,0.45,1024);
xsol = xsol(:);
Radius = 0.060; %0.060
Length = 0.120; % 0.090
Plus = Length./2 + xsol;
Minus = Length./2 - xsol;
Term1 = Plus./sqrt(Plus.^2 + Radius.^2);
Term2 = Minus./sqrt(Minus.^2 + Radius.^2);
LongCoil = Term1 + Term2;

mu0 = 4*pi.*10^(-7);
Current = 10;
NTurns = 5000;
NTurnsPerLength = NTurns./Length;
FieldCoil = 0.5.*mu0.*NTurnsPerLength.*Current.*LongCoil;


RadiusComp = 0.080; %0.100
% LengthComp = 0.240; %0.200
LengthComp = 0.240; %0.200
PlusComp = LengthComp./2 + xsol;
MinusComp = LengthComp./2 - xsol;
Term1Comp = PlusComp./sqrt(PlusComp.^2 + RadiusComp.^2);
Term2Comp = MinusComp./sqrt(MinusComp.^2 + RadiusComp.^2);
CompCoil = Term1Comp + Term2Comp;
% CompCoil = 0.24.*CompCoil./max(CompCoil);


% mu0 = 4*pi.*10^(-7);
CurrentComp = 10;
NTurnsComp = 0.200.*(LengthComp/Length).*NTurns;
NTurnsPerLengthComp = NTurnsComp./LengthComp;
FieldComp = 0.5.*mu0.*NTurnsPerLengthComp.*CurrentComp.*CompCoil;

NewCoil = FieldCoil-FieldComp;
% figure;
% subplot(2,1,1);
% plot(xsol,FieldCoil,xsol,FieldComp);
% subplot(2,1,2);
% plot(xsol,NewCoil);

figure;
plot(xsol,FieldCoil,xsol,FieldComp,xsol,NewCoil);
legend('Main coil','|Compensation coil|','Solenoid field');
xlabel('z (m)');
ylabel('B_z on axis (T)');
title('Solenoid field at 10 A');
axis([-0.4 0.4 -0.2 0.4])
line([-0.06;0.06],[-0.06;-0.06]);
line([-0.12;0.12],[-0.08;-0.08]);
line([-0.15;-0.15;+0.15;+0.15;-0.15],[-0.055;-0.15;-0.15;-0.055;-0.055]);
