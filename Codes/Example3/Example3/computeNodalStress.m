function stressAtNodes = computeNodalStress(model, displacement, E, nu)
    % Compute the stress at the nodes based on the given displacement
    %
    % Inputs:
    % model - PDEModel object containing the mesh information
    % displacement - Displacement vector (2*numNodes x 1)
    % E - Young's modulus
    % nu - Poisson's ratio
    %
    % Output:
    % stressAtNodes - Stress at nodes (3 x numNodes), where each column
    %                 contains [sigma_xx; sigma_yy; sigma_xy]

    % Extract mesh information
    nodes = model.Mesh.Nodes;  % Node coordinates (2 x numNodes)
    elements = model.Mesh.Elements;  % Element connectivity (3 x numElements)
    numNodes = size(nodes, 2);  % Total number of nodes
    numElements = size(elements, 2);  % Total number of elements

    % Initialize stress arrays
    stressAtNodes = zeros(3, numNodes);  % Stress at each node
    nodeContributionCount = zeros(1, numNodes);  % Number of elements contributing to each node

    % Plane stress material matrix D
    D = (E / (1 - nu^2)) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];

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

        % Extract the displacement vector for the current element
        u_elem = displacement([2*nodeIndices(1)-1, 2*nodeIndices(1), ...
                               2*nodeIndices(2)-1, 2*nodeIndices(2), ...
                               2*nodeIndices(3)-1, 2*nodeIndices(3)]);

        % Compute the strain in the element
        strain = B * u_elem;

        % Compute the stress in the element using Hooke's law
        stress = D * strain;

        % Distribute the stress to the nodes of the element
        for i = 1:3
            stressAtNodes(:, nodeIndices(i)) = stressAtNodes(:, nodeIndices(i)) + stress;
            nodeContributionCount(nodeIndices(i)) = nodeContributionCount(nodeIndices(i)) + 1;
        end
    end

    % Average the stress at each node
    for i = 1:numNodes
        if nodeContributionCount(i) > 0
            stressAtNodes(:, i) = stressAtNodes(:, i) / nodeContributionCount(i);
        end
    end
end
