% Units are [kN] and [mm]
% Output
% Reaction forces,nodal displacements (+ up/right, ? down/left)
% Member basic forces (+ tension, ? compression)
clear all
close all
clc
format short G

% Step 1: Identify global DOFs, set inputs, set nodal coordinates
nEle = 9 ; % number of elements
fDOF = 7 ; % number of free DOFs
E = 200 * ones(nEle,1) ; % elastic modulus of each element
A = 3500 * ones(nEle,1) ; % cross?sectional area of each element
alpha = 6e-06 * ones(nEle,1) ; % coefficient of thermal expansion of each element
% Element 1 2 3 4 5 6 7 8 9
x = 1000 * [0 2 ; 2 3 ; 3 2 ; 2 0 ; 0 -1 ; -1 0 ; 0 2 ; 2 0 ; -1 3] ; % xi, xj coordinate
y = 1000 * [0 0 ; 0 2 ; 2 4 ; 4 4 ; 4 2 ; 2 0 ; 0 4 ; 0 4 ; 2 2] ; % yi yj coordinate

% Step 2: Build element stiffness matrices & assemble structure stiffness matrix
L = zeros(nEle,1) ; % initialize length of each element
c = zeros(nEle,1) ; % initialzie cosine of each element
s = zeros(nEle,1) ; % initialize sine of each element

Grbm = zeros(1,4,nEle) ; % initialize rigid body modes matrix of each element
Grot = zeros(4,4,nEle) ; % initialize rotation matrix of each element

KEleBasic = zeros(1,1,nEle) ; % initialize element stiffness matrix (basic)
KEleLocal = zeros(4,4,nEle) ; % initialize element stiffness matrix (local)
KEleGlobal = zeros(4,4,nEle) ; % initialize element stiffness matrix (global)

for n = 1 : nEle
L(n) = sqrt( (x(n,2) - x(n,1))^2 + (y(n,2) - y(n,1))^2 ) ; % length of each element

c(n) = (x(n,2) - x(n,1) ) / L(n) ; % cosine of each element
s(n) = (y(n,2) - y(n,1) ) / L(n) ; % sine of each element

Grbm(:,:,n) = [-1 0 1 0] ;% rigid body modes matrix of each element
Grot(:,:,n) = [ c(n) s(n) 0 0 ; -s(n) c(n) 0 0 ; 0 0 c(n) s(n) ; 0 0 -s(n) c(n)] ; % rotation matrix of each element

KEleBasic(:,:,n) = E(n)*A(n)/L(n) ; % element stiffness matrix of each element (basic)
KEleLocal(:,:,n) = Grbm(:,:,n)' * KEleBasic(:,:,n) * Grbm(:,:,n) ; % element stiffness matrix of each element (local)
KEleGlobal(:,:,n) = Grot(:,:,n)' * KEleLocal(:,:,n) * Grot(:,:,n) ; % element stiffness matrix of each element (globals)

end

% ID array
ID = [ 8 9 1 10 ; % Element 1
1 10 11 12 ; % Element 2
11 12 2 3 ; % Element 3
2 3 4 5 ; % Element 4
4 5 6 7 ; % Element 5
6 7 8 9 ; % Element 6
8 9 2 3 ; % Element 7
1 10 4 5 ; % Element 8
6 7 11 12 ; % Element 8
] ;

structSize = max(max(ID)) ; % total number of DOFs (free + fixed)
Kstr = zeros(structSize,structSize) ;% initialize structure stiffness matrix

for n = 1 : nEle
Kstr(ID(n,:),ID(n,:)) = Kstr(ID(n,:),ID(n,:)) + KEleGlobal(:,:,n) ; % build total structure stiffness matrix
end

% partitioned structure stiffness matrices
Kff = Kstr(1:fDOF ,1:fDOF) ;% free DOFs
Kdf = Kstr(fDOF+1:end ,1:fDOF) ;% mixed terms
Kfd = Kstr(1:fDOF ,fDOF+1:end) ;% mixed terms
Kdd = Kstr(fDOF+1:end ,fDOF+1:end) ;% fixed DOFs

% Step 3: Build load vector
Pf = [0 ; 0 ; 0 ; 0 ; 0 ; 0 ; -4] ;% specified nodal loads
P0 = [zeros(length(Kstr),1)] ; % initialize vector of structure fixed?end forces (free and fixed DOFs)

FFglobal = zeros(4,nEle) ;% initialize element fixed?end forces (global)
FFlocal = zeros(4,nEle) ;% initialize element fixed?end forces (local)
FFbasic = zeros(1,nEle) ;% initialize element fixed?end forces (basic)

dT = [0 0 100 50 0 0 0 0 0] ;% temperature change on each element
w = [0 0 0 0 0.5 0 0 0 0]/1000 ; % distributed axial load on each element

for n = 1:nEle
FFbasic(:,n) = KEleBasic(:,:,n) * -alpha(n)*dT(n)*L(n) ;% element fixed?end force (basic - temperature only)
FFglobal(:,n) = Grot(:,:,n)' * [1 0 1 0]' * w(n)*L(n)/2 + Grot(:,:,n)' * Grbm(:,:,n)' * FFbasic(:,n) ; % element fixed?end forces(global)

FFlocal(:,n) = Grot(:,:,n) * FFglobal(:,n) ;% element fixed?end forces (local)
P0(ID(n,:)) = P0(ID(n,:)) + FFglobal(:,n) ;% update P0 to include-fixed?end forces of element n
end

Ud = [0 ; -3 ; 0 ; 0 ; 0] ;% support settlements
Pf0 = P0(1:fDOF) + Kfd*Ud ; % FFs acting on free DOFs
Pd0 = P0(fDOF+1:end) ;% FFs acting on fixed DOFs

% Step 4: Solve for structure nodal displacements
Uf = Kff\(Pf - Pf0) ; % solve for disp. at free DOFs
U = [Uf ; Ud] ;% total disp. vector

% Step 5: Solve for element forces
UEleGlobal = zeros(4,nEle) ; % initialize element nodal displacements (global)
UEleLocal = zeros(4,nEle) ; % initialize element nodal displacements (local)
UEleBasic = zeros(1,nEle) ; % initialize element nodal displacements (basic)

FEleGlobal = zeros(4,nEle) ; % initialize element end forces (global)
FEleLocal = zeros(4,nEle) ; % initialize element end forces (local)
FEleBasic = zeros(1,nEle) ; % initialize element end forces (basic)

for n = 1 : nEle
for j = 1 : 4
    
UEleGlobal(j,n) = U(ID(n,j)) ; % element nodal displacements, including support settlement (global)

end
FEleGlobal(:,n) = KEleGlobal(:,:,n) * UEleGlobal(:,n) + FFglobal(:,n) ; % element end forces(global)
UEleLocal(:,n) = Grot(:,:,n) * UEleGlobal(:,n) ; % element nodal displacements (local)

FEleLocal(:,n) = KEleLocal(:,:,n) * UEleLocal(:,n) + FFlocal(:,n) ;% element end forces (local)
UEleBasic(:,n) = Grbm(:,:,n) * UEleLocal(:,n) ;% element nodal displacements (basic)

FEleBasic(:,n) = KEleBasic(:,:,n) * UEleBasic(:,n) + FFbasic(:,n) ;% element end forces (basic)
end

% Step 6: Solve for reactions
Pd = Kdf*Uf + Kdd*Ud + Pd0 ;
% display results
% display element Kff
disp('Element stiffness matrices (global coordinates):')
disp(' ')
for n = 1 : nEle
kEleMsg = ['K_Ele',num2str(n),' = '] ;
disp(kEleMsg)
disp(' ')
disp(KEleGlobal(:,:,n))
end
% display Kff
disp('Structure stiffness matrix at free DOFs:')
Kff
disp(' ')
% display laod vector
disp('Load vectors at free DOFs:')
Pf
disp(' ')
Pf0
disp(' ')
% display nodal displacements
disp('Nodal displacements:')
for j = 1 : fDOF
dispMsg = [' ','DOF ',num2str(j),': ',sprintf('%.4f',Uf(j)),' mm'] ;
disp(dispMsg)
end
disp(' ')
% display reaction forces
disp('Reaction forces:')
disp([' ','R1X = ',sprintf('%.1f',Pd(1)),' kN'])
disp([' ','R1Y = ',sprintf('%.1f',Pd(2)),' kN'])
disp([' ','R2Y = ',sprintf('%.1f',Pd(3)),' kN'])
disp([' ','R3X = ',sprintf('%.1f',Pd(4)),' kN'])
disp([' ','R3Y = ',sprintf('%.1f',Pd(5)),' kN'])
disp(' ')
% display element forces
disp('Bar forces:')
for n = 1 : nEle
forceMsg = [' ','Element ',num2str(n),': ',sprintf('%.1f',FEleBasic(n)),'kN'] ;
disp(forceMsg)
end
disp(' ')

