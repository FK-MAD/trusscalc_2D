%% 2D Truss Calculator
% Solves a plannar (2D) truss consisting on n nodes and N members.
%
% INPUT: Coordinates of nodes, connectivity, Young's modulus, cross
% sectional area for every member, loads, supports, scaling constant for
% plots
% OUTPUT: Displacements and reaction forces for every truss node. Force, strain and stress for every truss member
%
% Author: Filippos Katsimalis, fkatsimalis@uth.gr, fkatsimalis@gmail.com

close all; clear; clc;
%% INPUT: Geometry data

% INPUT: coordinates of truss nodes in meters [m]. "coord" is a nx2 matrix
% where cells (i,1) and (i,2) hold the x and y coordinates, respectively,
% of node i in meters [m].
coord=...
    [0,0
    .5,.5
    1,0
    1.5,.5
    2,0
    2.5,.5
    3,0
    3.5,.5
    4,0];

% INPUT: connectivity of truss members. "conn" is a Nx2 matrix where cells
% (i,1) and (i,2) hold the nodes that connect to truss member i.
conn=...
    [1,2
    1,3
    2,3
    2,4
    3,4
    3,5
    4,5
    4,6
    5,6
    5,7
    6,7
    6,8
    7,8
    7,9
    8,9];


conn=sort(conn,2);

n=size(coord,1); % number of truss nodes
N=size(conn,1); % number of truss members

figure
plot_truss(coord,conn,n,N,'-',1) % plot of unloaded truss with numbered nodes and members

%% INPUT: Young's modulus and area
E=200*10^9; % INPUT: Young's modulus [Pa] of all truss members
A=10^-2*ones(N,1); % INPUT: cross sectional area for each truss member [m^2]. "A" is a vector of length N where element i
% holds the cross sectional area of truss member i

%% INPUT: BCs
force=zeros(n,2); % INPUT: x,y applied load at every truss node [N] (all zero by default). "force" is a nx2 matrix where cells
% (i,1) and (i,2) hold the x and y force, respectively, applied on node i
% in Newton [N].
force(3,:)=[0,-50000]; % Fx_3=0, Fy_3=-500000 [N]
force(5,:)=[0,-50000]; % Fx_5=0, Fy_5=-500000 [N]
force(7,:)=[0,-50000]; % Fx_7=0, Fy_7=-500000 [N]

disp=inf(n,2); % INPUT: x,y specified dispacement of every truss node [m] (inf where there is no specified diplacement - all inf by default)
% "disp" is a nx2 matrix where cells (i,1) and (i,2) hold the x and y displacement, respectively, 
% applied on node i in meters [m]. By default, all nodes are free to move (inf displacement).
disp(1,:)=0; % node 1 -> pinned (x=0,y=0)
disp(9,2)=0; % node 9 -> x roller (x=inf,y=0)

%% Calculations
[K,L,theta]=K_calc(coord,conn,n,N,E,A); % calculation of global matrix K, length and angle for every truss member

[D,F]=DF_calc(disp,force,K,n); % calculation of displacements and forces for each truss node

[member_force,strain,stress]=fss_calc(conn,D,E,A,N,L,theta); % calculation of force, strain and stress for each truss member

%% Results
write_results(n,N,D,F,E,A,member_force,strain,stress); % write results in a text file

%% Plots & INPUT: Scaling constant for plots
figure
plot_truss(coord,conn,n,N,':',0); % plot of the unloaded truss
scale=500; % INPUT: scaling constant to display displacements more clearly
plot_truss(coord+D*scale,conn,n,N,'-',0); % loaded truss with scaled displacements

%% Functions
function [K,L,theta]=K_calc(coord,conn,n,N,E,A)
K=zeros(2*n);
L=zeros(N,1);
theta=zeros(N,1);
m=zeros(N,1);
for k=1:n
    [r,~]=find(conn==k);
    k_adj=conn(r,:);
    k_adj(k_adj==k)=[];
    for adj=k_adj
        x=coord(k,1);
        y=coord(k,2);
        x_adj=coord(adj,1);
        y_adj=coord(adj,2);
        
        member=sum(sort([k,adj])==conn,2)==2;
        
        if m(member)==0
            m(member)=1;
            L(member)=sqrt((x_adj-x)^2+(y_adj-y)^2);
            inpr_x=(x_adj-x)/L(member);
            inpr_y=(y_adj-y)/L(member);
            if inpr_y/inpr_x==-inf
                theta(member)=atan(inf);
            else
                theta(member)=atan(inpr_y/inpr_x);
            end
        else
            inpr_x=(x_adj-x)/L(member);
            inpr_y=(y_adj-y)/L(member);
        end
        
        alpha=zeros(1,2*n);
        alpha(2*k)=-inpr_y;
        alpha(2*k-1)=-inpr_x;
        alpha(2*adj)=inpr_y;
        alpha(2*adj-1)=inpr_x;
        alpha=-alpha;
        
        K(2*k-1,:)=K(2*k-1,:)+(E*A(member)*alpha*inpr_x/L(member));
        K(2*k,:)=K(2*k,:)+(E*A(member)*alpha*inpr_y/L(member));
    end
end
end

function [D,F]=DF_calc(disp,force,K,n)
disp_lin=reshape(disp',[],1)';
force_lin=reshape(force',[],1)';

f_disp=find(disp_lin==inf);
n_f=length(f_disp);

s_disp=setdiff(1:2*n,f_disp);

swap=[f_disp,s_disp];

K_swaped=K(swap,swap);
K_ff=K_swaped(1:n_f,1:n_f);
K_fs=K_swaped(1:n_f,n_f+1:end);
K_ss=K_swaped(n_f+1:end,n_f+1:end);

F=force_lin(swap)';
F_f=F(1:n_f);

D_s=disp_lin(s_disp)';

D_f=K_ff\(F_f-K_fs*D_s);
F_s=K_fs'*D_f+K_ss*D_s;

D=zeros(1,2*n);
D(s_disp)=disp_lin(s_disp);
D(f_disp)=D_f;
D=reshape(D,2,n)';

F=zeros(1,2*n);
F(s_disp)=F_s;
F(f_disp)=force_lin(f_disp);
F=reshape(F,2,n)';

end

function write_results(n,N,D,F,E,A,member_force,strain,stress)
fid=fopen('results.txt','w');

fprintf(fid,'%10s %20s %20s %15s %15s\n',...
    'Node','x displacement [m]','y displacement [m]','x Force [N]','y Force [N]');
fprintf(fid,'%10d %20.4g %20.4g %15.4g %15.4g\n',[1:n;D';F']);

fprintf(fid,'\n%s\n\n',repmat('-',1,85));

fprintf(fid,'%10s %12s %12s %15s %15s %15s\n',...
    'Member','Area [m^2]','E [Pa]','Force [N]','Strain','Stress [Pa]');
fprintf(fid,'%10d %12.2e %12.2e %15.4g %15.4g %15.4g\n',...
    [1:N;A';repmat(E,1,N);member_force';strain';stress']);

fclose(fid);
end

function [member_force,strain,stress]=fss_calc(conn,D,E,A,N,L,theta)
stress=zeros(N,1);
strain=zeros(N,1);
for k=1:N
    th=theta(k);
    n_1=conn(k,1);
    n_2=conn(k,2);
    u_1=D(n_1,1)*cos(th)+D(n_1,2)*sin(th);
    u_2=D(n_2,1)*cos(th)+D(n_2,2)*sin(th);
    strain(k)=(u_2-u_1)/L(k);
    stress(k)=E*strain(k);
end
stress=E*strain;
member_force=stress.*A;
end

function plot_truss(coord,conn,n,N,style,txt)
x_offset=0.01*max(coord(:,1));
y_offset=0.01*max(coord(:,2));
hold on
for k=1:N
    x=[coord(conn(k,1),1),coord(conn(k,2),1)];
    y=[coord(conn(k,1),2),coord(conn(k,2),2)];
    line(x,y,'linestyle',style,'color','k','linewidth',1);
    if txt==1
        text(x_offset+x(1)+(x(2)-x(1))/2,y_offset+y(1)+(y(2)-y(1))/2,num2str(k),...
            'color','k','fontsize',15);
    end
end

for k=1:n
    x=coord(k,1);
    y=coord(k,2);
    plot(x,y,'or','markerfacecolor','r','linewidth',0.5);
    if txt==1
        text(x_offset+x,y_offset+y,num2str(k),'color','r','fontsize',15);
    end
end
xlabel('\it x [m]','fontsize',20)
ylabel('\it y [m]','fontsize',20)
axis equal

end