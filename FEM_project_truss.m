close all
clear all 
format long
clc
ele_nod=[2 1;1 3];
%number of elements
num_ele=size(ele_nod,1);
%number of nodes
num_nod=3;
% nodes coordinates
nod_coor=[0 0;-0.5 0;0.4 -0.3];
% elements degree of freedom (DOF) 
ele_dof=[3 4 1 2;1 2 5 6];

rho = 7800;
A = 750e-06;
E = 200e09;
L = 0.5;
% initial zero matrix for all matrices
displacement=zeros(2*num_nod,1);
force=zeros(2*num_nod,1);
stiffness=zeros(2*num_nod);
stress = zeros(num_nod - 1, 1);
mass = zeros(2*num_nod);
%applied loads at DOFs
force(2) = -20e03;
displacement(3,1) = 0.0;
displacement(4,1) = 0.0;
displacement(5,1) = 0.0;
displacement(6,1) = 0.0;

% computation of the system stiffness matrix
for e=1:num_ele
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L;
 C
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L;
 S
 k=(A*E/L)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;...
   -C*C -C*S C*C C*S; -C*S -S*S C*S S*S];
 k
 m = (rho*A*L/6)*[2 0 1 0;0 2 0 1;1 0 2 0;0 1 0 2];
 m
% extract the rows of ele_dof (for each element e)
 ele_dof_vec=ele_dof(e,:);
 ele_dof_vec
    for i=1:4
        for j=1:4
            stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))= stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);
            mass(ele_dof_vec(1,i),ele_dof_vec(1,j))= mass(ele_dof_vec(1,i),ele_dof_vec(1,j))+ m(i,j);
        end
    end
end

% known force array
known_f_a=[1;2];
for i=1:size(known_f_a,1)
   dis_new(i,1)=displacement(known_f_a(i,1),1);
   force_new(i,1)=force(known_f_a(i,1),1);
end
for i=1:size(known_f_a,1)
    for j=1:size(known_f_a,1)
    stiff_new(i,j)=stiffness(known_f_a(i,1),known_f_a(j,1));
end
end
% solving the partitioned matrix 
dis_new=linsolve(stiff_new,force_new);
for i=1:size(known_f_a,1)
  displacement(known_f_a(i,1),1)=dis_new(i,1);
end

%Plot undeformed graph
for e=1:num_ele
    x=[nod_coor(ele_nod(e,1),1) nod_coor(ele_nod(e,2),1)];
    y=[nod_coor(ele_nod(e,1),2) nod_coor(ele_nod(e,2),2)];
plot(x,y,'b')
hold on
end

% known dicplacement array
known_dis_a=[3;4;5;6];
for i=1:size(known_dis_a,1)
    force(known_dis_a(i,1),1)=stiffness(known_dis_a(i,1),:)...
    *displacement;
end
% stress in elements
for e=1:num_ele
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L;
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L;
 stress(e)=(E/L)*[-C -S C S]*displacement((ele_dof(e,:))');
end
disp('stiffness')
stiffness
disp('displacement')
displacement
disp('force')
force
disp('stress')
stress'
k=0;
for i=1:3
    for j=1:2
    k=k+1;
    nod_coor_def(i,j)=nod_coor(i,j)+ displacement(k,1);
    end
end
% deformed truss plot
for e=1:num_ele
    x=[nod_coor_def(ele_nod(e,1),1) nod_coor_def(ele_nod(e,2),1)];
    y=[nod_coor_def(ele_nod(e,1),2) nod_coor_def(ele_nod(e,2),2)];
plot(x,y,'r')
hold on
end