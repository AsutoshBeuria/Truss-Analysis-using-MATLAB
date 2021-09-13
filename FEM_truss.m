% 2D truss program                                                            
% this program solves the 2D truss problems using Finite Element Method (FEM)                                        
% Copyright Farzad Mohebbi
% created by Farzad Mohebbi April 2018  
% https://www.researchgate.net/profile/Farzad_Mohebbi
% see attached PDF file for more details
% This is a free program. 
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY.
% This program can be modified based on user requirements.
clc
clear all % clear memory
close all
%tic       % starts a stopwatch timer
% elements nodes
ele_nod=[1 3;3 2;3 4;2 4];
%number of elements
num_ele=size(ele_nod,1);
%number of nodes
num_nod =4;
% nodes coordinates
nod_coor=[0 0;100 0;50 50;200 100];
% elements degree of freedom (DOF) 
ele_dof=[1 2 5 6;5 6 3 4;5 6 7 8;3 4 7 8];
% A, E, L are cross sectional area, Young's modulus, length of elements,respectively.
A = [2; 2; 1; 1];
E = 2e06;
% initial zero matrix for all matrices
displacement=zeros(2*num_nod,1);
force=zeros(2*num_nod,1);
stiffness=zeros(2*num_nod);
%applied loads at DOFs
force(8) = -1e03;
force(7) = -1e05;
%Boundary conditions
displacement (1,1)=0.0;
displacement (2,1)=0.0;
displacement (3,1)=0.0;
displacement(4,1) = 0.0;
% computation of the system stiffness matrix
for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 k=(A(e)*E/L(e)*[C*C C*S -C*C -C*S;C*S S*S -C*S -S*S;...
   -C*C -C*S C*C C*S; -C*S -S*S C*S S*S]);
   
% extract the rows of ele_dof (for each element e)
ele_dof_vec=ele_dof(e,:);
    for i=1:4
        for j=1:4
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))=...
  stiffness(ele_dof_vec(1,i),ele_dof_vec(1,j))+k(i,j);
        end
    end
end
% known force array
known_f_a=[5;6;7;8];
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
dis_new= stiff_new \ force_new;
for i=1:size(known_f_a,1)
  displacement(known_f_a(i,1),1)=dis_new(i,1);
end
% undeformed truss plot
for e=1:num_ele
    x=[nod_coor(ele_nod(e,1),1) nod_coor(ele_nod(e,2),1)];
    y=[nod_coor(ele_nod(e,1),2) nod_coor(ele_nod(e,2),2)];
plot(x,y,'b')
hold on
end
% known dicplacement array
known_dis_a=[1;2;3;4];
for i=1:size(known_dis_a,1)
    force(known_dis_a(i,1),1)=stiffness(known_dis_a(i,1),:)...
    *displacement;
end
% stress in elements
for e=1:num_ele
 L(e)=sqrt((nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))^2+...
      (nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))^2);
 C=(nod_coor(ele_nod(e,2),1)-nod_coor(ele_nod(e,1),1))/L(e);
 S=(nod_coor(ele_nod(e,2),2)-nod_coor(ele_nod(e,1),2))/L(e);
 stress(e)=(E/L(e))*[-C -S C S]*displacement((ele_dof(e,:))');
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
for i=1:4
    for j=1:2
    k=k+1;
    nod_coor_def(i,j)=nod_coor(i,j)+displacement(k,1);
    end
end
% deformed truss plot
for e=1:num_ele
    x=[nod_coor_def(ele_nod(e,1),1) nod_coor_def(ele_nod(e,2),1)];
    y=[nod_coor_def(ele_nod(e,1),2) nod_coor_def(ele_nod(e,2),2)];
plot(x,y,'r')
hold on
end
%toc