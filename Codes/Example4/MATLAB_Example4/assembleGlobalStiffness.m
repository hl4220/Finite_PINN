function K_global = assembleGlobalStiffness(model, E, nu)
    % Assemble the global stiffness matrix for 2D elasticity problem
    %
    % Inputs:
    % model - PDEModel object containing the mesh information
    % E - Young's modulus
    % nu - Poisson's ratio
    %
    % Output:
    % K_global - Global stiffness matrix

    % Extract mesh information
    nodes = model.Mesh.Nodes;  % Node coordinates (2 x numNodes)
    elements = model.Mesh.Elements;  % Element connectivity (3 x numElements)
    numNodes = size(nodes, 2);  % Total number of nodes
    numElements = size(elements, 2);  % Total number of elements

    % Plane stress material matrix D
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];

    label1 = zeros(1,36*numElements);
    label2 = zeros(1,36*numElements);
    values = zeros(1,36*numElements);

    % Loop over each element
    for iElem = 1:numElements
        % Get the node indices for the current element
        nodeIndices = elements(:, iElem);

        % Extract the coordinates of the nodes of the current element
        x1 = nodes(1, nodeIndices(1));
        y1 = nodes(2, nodeIndices(1));
        x2 = nodes(1, nodeIndices(2));
        y2 = nodes(2, nodeIndices(2));
        x3 = nodes(1, nodeIndices(3));
        y3 = nodes(2, nodeIndices(3));

        % Compute the area of the triangle element
        Area = 0.5 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]);

        % Compute the B matrix for the current element
        b1 = y2 - y3;
        b2 = y3 - y1;
        b3 = y1 - y2;
        c1 = x3 - x2;
        c2 = x1 - x3;
        c3 = x2 - x1;

        B = (1 / (2 * Area)) * [b1, 0, b2, 0, b3, 0;
                                0, c1, 0, c2, 0, c3;
                                c1, b1, c2, b2, c3, b3];

        % Compute the elemental stiffness matrix
        K_elem = Area * (B' * D * B);

        % Assemble the elemental stiffness matrix into the global stiffness matrix
        dof = [2*nodeIndices(1)-1, 2*nodeIndices(1), ...
               2*nodeIndices(2)-1, 2*nodeIndices(2), ...
               2*nodeIndices(3)-1, 2*nodeIndices(3)];

        [dof1,dof2] = meshgrid(dof,dof);

        label1(36*iElem-35:36*iElem) = reshape(dof1,[],1);
        label2(36*iElem-35:36*iElem) = reshape(dof2,[],1);
        values(36*iElem-35:36*iElem) = reshape(K_elem,[],1);
    end

    K_global = sparse(label1,label2,values,numNodes*2,numNodes*2);
end
