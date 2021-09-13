clear all
close all
clc

%Given:
E = 200*10^9; %Pa
A = 0.001; %m^2
Grbm = [-1 0 1 0]; %same for all truss


%Element a
%Node I = Node A
%Node J = Node B
xi = -8; yi = 4;
xj = 0; yj = 4;
ida = [3 ;4 ;5; 6];
La = sqrt((xj-xi)^2 + (yj-yi)^2);
c = (xj-xi)/La;
s = (yj-yi)/La;
Grota = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
ka = ((E*A)/La)*[c*c s*c -c*c -s*c; s*c s*s -s*c -s*s; -c*c -s*c c*c s*c; -s*c -s*s s*c s*s]
%Element b
%Node I = Node B
%Node J = Node C
xi = 0; yi = 4;
xj = 8; yj = 0;
idb = [7 8 5 6]';
Lb = sqrt((xj-xi)^2 + (yj-yi)^2);
c = (xj-xi)/Lb;
s = (yj-yi)/Lb;
Grotb = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
kb = ((E*A)/Lb)*[c*c s*c -c*c -s*c; s*c s*s -s*c -s*s; -c*c -s*c c*c s*c; -s*c -s*s s*c s*s]
%Element c
%Node I = Node C
%Node J = Node D
xi = 8; yi = 0;
xj = 0; yj = 0;
idc = [7 8 1 2]';
Lc = sqrt((xj-xi)^2 + (yj-yi)^2);
c = (xj-xi)/Lc;
s = (yj-yi)/Lc;
Grotc = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
kc = ((E*A)/Lc)*[c*c s*c -c*c -s*c; s*c s*s -s*c -s*s; -c*c -s*c c*c s*c; -s*c -s*s s*c s*s]
%Element d
%Node I = Node
%Node J = Node C
xi = 0; yi = 4;
xj = 0; yj = 0;
idd = [5 6 2 1]';

Ld = sqrt((xj-xi)^2 + (yj-yi)^2);
c = (xj-xi)/Ld;
s = (yj-yi)/Ld;

Grotd = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
kd = ((E*A)/Ld)*[c*c s*c -c*c -s*c; s*c s*s -s*c -s*s; -c*c -s*c c*c s*c; -s*c -s*s s*c s*s]
%Element e
%Node I = Node A
%Node J = Node D
xi = -8; yi = 4;
xj = 0; yj = 0;
ide = [3 4 1 2]';

Le = sqrt((xj-xi)^2 + (yj-yi)^2);
c = (xj-xi)/Le;
s = (yj-yi)/Le;
Grote = [c s 0 0; -s c 0 0; 0 0 c s; 0 0 -s c];
ke = ((E*A)/Le)*[c*c s*c -c*c -s*c; s*c s*s -s*c -s*s; -c*c -s*c c*c s*c; -s*c -s*s s*c s*s]

% partition to obtain Kff and show your output

K = zeros(8,8);
K(ida,ida) = K(ida,ida) + ka;
K(idb,idb) = K(idb,idb) + kb;
K(idc,idc) = K(idc,idc) + kc;
K(idd,idd) = K(idd,idd) + kd;
K(ide,ide) = K(ide,ide) + ke;

% Partition
Kff = K([5 6 2 1], [5 6 2 1])
Kdf = K([8 7 4 3], [5 6 2 1]);
Kfd = K([5 6 2 1], [8 7 4 3]);
Kdd = K([8 7 4 3], [8 7 4 3]);
% displacements of the free DOFs and the reactions at supports.

% Load Vectors:
Pf = [0;0;-48;0]

%Displacements Outputs, U1, U2, U5, U6
Uf = inv(Kff).*(Pf)
Ud = zeros(4,1);
U = [Uf, Ud];

%Reactions numbered according to restrained DOFs 8, 7, 4, 5
Pd = Kdf*Uf

% Compute member forces

%Member Forces Element a

deltaa = U(ida) %Element global Disp.
Fa = ka*deltaa %Element global forces
deltaBPa = Grbm*Grota*deltaa %Basic deformation Delta Bar Prime (BP)
FBPa = ((E*A)/La)*deltaBPa %Basic Force Element a

%Member Forces Element b

deltab = U(idb) %Element global Disp.
Fb = kb*deltab %Element global forces
deltaBPb = Grbm*Grotb*deltab %Basic deformation Delta Bar Prime (BP)
FBPb = ((E*A)/Lb)*deltaBPb %Basic Force Element b

%Memeber Forces Element c

deltac = U(idc) %Element global Disp.
Fc = kc*deltac %Element global forces
deltaBPc = Grbm*Grotc*deltac %Basic deformation Delta Bar Prime (BP)
FBPc = ((E*A)/Lc)*deltaBPc %Basic Force Element c

%Member Forces Element d

deltad = U(idd) %Element global Disp.
Fd = kd*deltad %Element global forces
deltaBPd = Grbm*Grotd*deltad %Basic deformation Delta Bar Prime (BP)
FBPd = ((E*A)/Ld)*deltaBPd %Basic Force Element d

%Member Forces Element e

deltae = U(ide) %Element global Disp.
Fe = ke*deltae %Element global forces
deltaBPe = Grbm*Grote*deltae %Basic deformation Delta Bar Prime (BP)
FBPe = ((E*A)/Le)*deltaBPe %Basic Force Element