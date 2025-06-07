clc
clear

%% Model and Mesh
V = [0 0; 1 0; 1 1; 0 1];
R1 = [3,size(V,1),V(:,1)',V(:,2)']';
C1 = [1,0.35,0.25,0.15]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
C2 = [1,0.75,0.5,0.10]';
C2 = [C2;zeros(length(R1) - length(C2),1)];
C3 = [1,0.3,0.7,0.20]';
C3 = [C3;zeros(length(R1) - length(C3),1)];
gm = [R1,C1,C2,C3];
sf = 'R1-C1-C2-C3';
ns = char('R1','C1','C2','C3');
ns = ns';
g = decsg(gm,sf,ns);
model = createpde();
geometryFromEdges(model, g); 
mesh = generateMesh(model,Hmax=0.0075,Hgrad=1.0,GeometricOrder='linear'); 
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
save('Package_Thermal_Expansion_Weak_Form','nodes','phi','dphi_dx','dphi_dy');

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


