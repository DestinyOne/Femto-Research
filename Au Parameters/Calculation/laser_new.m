clear;
clc;
%% CONSTANTS
h_bar=1.054d-34;         % J*s Plank's constant
Pi=3.141593;		      %
e_charge=1.602d-19;
m_e=9.11d-31;	          % kg electron effective mass
k_b=8.61733d-5; % eV/K
Na = 6.02214076e23; % mol^-1
epsilon0 = 8.854187817620d-12; % SI
v_c = 2.99792458e+8; % m/s
v_s = 3240; %m/s speed of sound
lambdal = 12.4e-7;%8.e-7; % m (laser wavelength)
f_laser = 2.*Pi*v_c/lambdal/1e15; % 1/fs

%% Material constants
newmat='Au';
Ti = 300; %K
T_me = 1/k_b;
T_melt = 1337; %account for over heating
A = 196.97; % g/mol
density = 19.25; % g/cm^3
n_ion = 1.e6*density/A*Na; %5.90e+28; %m-3
m_i = A*1.66054e-27; % kg
r0 = (3/4/pi/n_ion)^(1/3);%2.14e-10; % m, atomic radius
FE = h_bar^2/2/m_e*(3*pi^2*n_ion)^(2/3)/e_charge; % eV 
a0 = 4.065e-10; %lattice parameter
Tme = (2*pi*h_bar)^2/3/m_e/a0^2/e_charge/k_b;
% Permittivity
% epsilon1 = -22.104;
% epsilon2 = 1.78;
epsilon1 = -27.276;
epsilon2 = 1.0889;
%nu_eph0 = f_laser*epsilon2/(1-epsilon1);
nu_eph0 = 2.6e14;
%nu_eph0 = 1.23e11*300;

%% Names of files76
fZall = [newmat '_Zall.dat'];
fveff = [newmat '_veff.dat'];
lines = 1000;
%% Experimental data
Ted = [0.026, 0.63, 0.76, 0.89, 1.2, 1.4, 1.9, 2.4, 2.7, 3.4, 4.5, 4.9]';
Zeffd = [1, 1.06, 1.10, 1.14, 1.24, 1.33, 1.53, 1.72, 1.80, 2.09, 2.43, 2.58]';
veffd = [0.129, 0.5, 0.6, 0.8, 1.0, 1.2, 1.6, 2.1, 2.3, 2.9, 3.8, 4.2]';
datas = size(Ted);
datass = datas(1);
%% Temperature array
Te = zeros(lines,1);
Zeff = zeros(lines,1);
veff = zeros(lines,1);
Absorb = zeros(lines,1);
Reflect = zeros(lines,1);
length = zeros(lines,1);
phase = zeros(lines,1);
freflect = zeros(lines,1);
for i=1:lines
  Te(i) = 300 + (i-1)*50;
  Tek = Te(i)*k_b;
  k = 0;
  for j=1:datass-1
    if (Tek>Ted(j) && Tek<=Ted(j+1))
      k = j;
      break;
    end
  end
  if k == 0
    Zeff(i) = Zeffd(1);
    veff(i) = veffd(1);
  else
    Zeff(i) = Zeffd(k) + (Zeffd(k+1)-Zeffd(k))/(Ted(k+1)-Ted(k)) * (Tek-Ted(k));
    veff(i) = veffd(k) + (veffd(k+1)-veffd(k))/(Ted(k+1)-Ted(k)) * (Tek-Ted(k));
  end 
  
  
  f_plasma = sqrt(Zeff(i)*n_ion*e_charge^2/epsilon0/m_e)/1e15;
  Re_epsilon = 7.6 - f_plasma^2/(veff(i)^2 + f_laser^2);
  Im_epsilon = f_plasma^2*veff(i)/f_laser/(veff(i)^2 + f_laser^2);
  myk = sqrt(0.5*(sqrt(Re_epsilon*Re_epsilon + Im_epsilon*Im_epsilon) - Re_epsilon));
  myn = sqrt(0.5*(sqrt(Re_epsilon*Re_epsilon + Im_epsilon*Im_epsilon) + Re_epsilon));
  
  length(i) = v_c/(f_laser*1e15)/myk/2.0;
  Absorb(i) = 4*myn/( (myn+1)^2 + myk^2 );
  Reflect(i) = 1 - Absorb(i);
  phase(i) = atan(2*myk/(myn^2+myk^2-1));
  freflect(i) = sqrt(Reflect(i)/cos(phase(i)));
  
end


Nout = [newmat '_Reflec.dat'];
outpt = fopen(Nout,'w');
fprintf(outpt,'%d\n',lines);
for ii=1:lines
  fprintf(outpt,'%f %e\n',Te(ii),Reflect(ii));
end
fclose(outpt);

Nout = [newmat '_Pen.dat'];
outpt = fopen(Nout,'w');
fprintf(outpt,'%d\n',lines);
for ii=1:lines
  fprintf(outpt,'%f %e\n',Te(ii),length(ii)*1.e10);
end
fclose(outpt);

Nout = [newmat '_Z.dat'];
outpt = fopen(Nout,'w');
fprintf(outpt,'%d\n',lines);
for ii=1:lines
  fprintf(outpt,'%f %e\n',Te(ii),Zeff(ii));
end
fclose(outpt);

hold on;
grid();
set(gca,'FontSize',20);
title('Optical Parameters');
yyaxis left;
ylabel('R','FontWeight','bold','Fontsize',22);
plot(Te*k_b,Reflect,'LineWidth',3);
yyaxis right;
ylabel('l_p(nm)','FontWeight','bold','Fontsize',22);
plot(Te*k_b,length*1.e9,'--','LineWidth',3);
xlabel('Elctron Temperature(eV)','FontWeight','bold','Fontsize',22);
