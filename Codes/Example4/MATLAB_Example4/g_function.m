function Z = transfinite_interpolation_one_side(X, Y, Z_boundary)
    % X, Y: Meshgrid points in the domain
    % Z_boundary: Values along the left boundary (size ny)
    %
    % Z: Interpolated values in the domain

    [ny, nx] = size(X); % Get the size of the grid

    % Initialize the interpolated grid
    Z = zeros(ny, nx);

    % Perform transfinite interpolation
    for i = 1:ny
        for j = 1:nx
            xi = (j-1) / (nx-1);   % Parameter along the horizontal axis

            % Transfinite interpolation using only the left boundary
            Z(i, j) = (1 - xi) * Z_boundary(i);
        end
    end
end

% Example usage
nx = 50; % Number of points along x-axis
ny = 50; % Number of points along y-axis

x = linspace(0, 1, nx); % x-coordinates
y = linspace(0, 1, ny); % y-coordinates

[X, Y] = meshgrid(x, y); % Generate the grid

% Boundary values (just an example, can be any given data)
Z_boundary = linspace(0, 1, ny); % Values at the left boundary (e.g., linearly increasing)

% Perform transfinite interpolation
Z = transfinite_interpolation_one_side(X, Y, Z_boundary);

% Plot the interpolated surface
figure;
mesh(X, Y, Z);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('2D Transfinite Interpolation from Left Boundary');
