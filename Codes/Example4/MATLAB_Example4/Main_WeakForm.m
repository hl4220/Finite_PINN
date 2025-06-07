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

%% Stiffness matrix
K = sparse(size(nodes,1),size(nodes,1)); 
for e = 1:size(elements,1)
  nc = elements(e,1:3);
  Ke = LKe(nodes(nc,:)); 
  K(nc,nc)= K(nc,nc)+Ke;
end  

%% General eigenfunctions
numberofbases = 4;
[e,~] = eigs(sparse(K),numberofbases, 'smallestabs');
e = e*10; 

%% Information Collection
Nn = length(nodes);
phi = zeros(Nn,numberofbases);
dphi_dx = zeros(Nn,numberofbases);
dphi_dy = zeros(Nn,numberofbases);

coords(:,4) = 0; 
for i = 1:numberofbases
    p = e(:,i);
    [dp_dx,dp_dy] = derivs(nodes,elements,p);
    phi(:,i) = p;
    dphi_dx(:,i) = dp_dx;
    dphi_dy(:,i) = dp_dy;
end

%% Saving
save('Package_Open_Notch_Weak_Form','nodes','phi','dphi_dx','dphi_dy');

%% Plot
figure(1);
pdeplot(model, 'NodeLabels', 'off');
axis equal; axis tight; axis off;
title('Mesh');
figure(2);
tiledlayout(1,4)
for i = 1:numberofbases
    disp(i)
    nexttile(i)
    patch('Faces',elements(:,1:3),'Vertices',nodes,...
          'FaceVertexCData', phi(:,i), 'FaceColor', 'interp', 'edgecolor','None');
    colormap jet; title(['$\mathbf{\phi_{',num2str(i),'}}$'],'Interpreter','latex')
    axis equal; axis tight; axis off; hold on
end


