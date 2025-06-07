clc
clear

%% Model and Mesh
gap = 0.02;
r = 0.1;
theta = linspace(2*pi,0,100)';
xr = r*cos(theta);
yr = r*sin(theta);
efflc = ~(abs(yr)<gap/2 & xr>0);
xyr = [xr(efflc) + 0.5,yr(efflc) + 0.5];
V = [0 0; 1 0; 1 0.5-gap/2; 0.5+r 0.5-gap/2; xyr ; 0.5+r 0.5+gap/2; 1 0.5+gap/2; 1 1; 0 1];
R1 = [3,size(V,1),V(:,1)',V(:,2)']';
model = createpde();
geometryFromEdges(model, decsg(R1)); 
mesh = generateMesh(model,Hmax=0.01,Hgrad=2.0,GeometricOrder='linear'); 
nodes = mesh.Nodes';
elements = mesh.Elements';

%% Boundary Information
boundary = findNodes(model.Mesh, 'region', 'Edge', 1:model.Geometry.NumEdges);
bc1_nodes = find(nodes(:,2)==0);
bc2_nodes = find(nodes(:,2)==1);
free_nodes = setdiff(boundary,[bc1_nodes;bc2_nodes]);
interior_nodes = setdiff(1:size(nodes,1),boundary);
Edges = unique(sort([elements(:,[1,2]);elements(:,[2,3]);elements(:,[1,3])],2),'rows');
bc1_edges = Edges(ismember(Edges(:,1),bc1_nodes) & ismember(Edges(:,2),bc1_nodes),:);
bc2_edges = Edges(ismember(Edges(:,1),bc2_nodes) & ismember(Edges(:,2),bc2_nodes),:);
free_edges = Edges(ismember(Edges(:,1),free_nodes) & ismember(Edges(:,2),free_nodes),:);
boundary_edges = Edges(ismember(Edges(:,1),boundary) & ismember(Edges(:,2),boundary),:);

%% Stiffness Matrix
K = sparse(size(nodes,1),size(nodes,1)); 
for e = 1:size(elements,1)
  nc = elements(e,1:3);
  Ke = LKe(nodes(nc,:)); 
  K(nc,nc)= K(nc,nc)+Ke;
end  

%% General eigenfunctions
numberofbases = 8;
[e,~] = eigs(sparse(K),numberofbases, 'smallestabs');
e = e*10; 

%% Specific eigenfunctions
bc = [bc1_nodes;bc2_nodes];

Ku = K;
Ku(bc,:) = 0;
Ku(bc,bc) = eye(length(bc));

[eu,v] = eigs(sparse(Ku),numberofbases, 'smallestabs');
eu = eu*10;

%% Hybrid eigenfunctions
e = e + eu;

%% Information Collection
Nn = length(nodes);
phi = zeros(Nn,numberofbases);
dphi_dx = zeros(Nn,numberofbases);
dphi_dy = zeros(Nn,numberofbases);
d2phi_dx2 = zeros(Nn,numberofbases);
d2phi_dy2 = zeros(Nn,numberofbases);

coords(:,4) = 0; 
for i = 1:numberofbases
    p = e(:,i);
    [dp_dx,dp_dy] = derivs(nodes,elements,p);
    [d2p_dx2,d2p_dxdy] = derivs(nodes,elements,dp_dx);
    [d2p_dydx,d2p_dy2] = derivs(nodes,elements,dp_dy);
    phi(:,i) = p;
    dphi_dx(:,i) = dp_dx;
    dphi_dy(:,i) = dp_dy;
    d2phi_dx2(:,i) = d2p_dx2;
    d2phi_dy2(:,i) = d2p_dy2;
end

%% FEM Results
E = 1; mu = 0.3;
K = assembleGlobalStiffness(model, E, mu);
bc = [2*bc1_nodes-1;2*bc1_nodes];
mbc = 2*bc2_nodes;
K(bc,:) = 0;
K(:,bc) = 0;
K(bc,bc) = eye(length(bc));
F = sparse(2*size(nodes,1),1);
F(mbc) = 1.0/1000;
F([mbc(1),mbc(2)]) = F([mbc(1),mbc(2)])/2;
U = full(K\F);
U1 = U(1:2:end);U2 = U(2:2:end);
U = [U1,U2];
L_label = [5534 389 5080 7385 8579 10972 5413 2781 7872 9832];
U_label = U(L_label,:);

%% Saving
save('Package_Open_Notch_Strong_Form', ...
    'nodes','bc1_nodes','bc2_nodes','free_nodes','interior_nodes', ...
    'phi','dphi_dx','dphi_dy','d2phi_dx2','d2phi_dy2',...
    'L_label','U_label','U');

%% Plot
figure(1);
pdeplot(model, 'NodeLabels', 'off');
axis equal; axis tight; axis off;
title('Mesh');
figure(2);
tiledlayout(2,4)
for i = 1:numberofbases
    disp(i)
    nexttile(i)
    patch('Faces',elements(:,1:3),'Vertices',nodes,...
          'FaceVertexCData', phi(:,i), 'FaceColor', 'interp', 'edgecolor','None');
    colormap jet; title(['$\mathbf{\phi_{',num2str(i),'}}$'],'Interpreter','latex')
    axis equal; axis tight; axis off; hold on
end


