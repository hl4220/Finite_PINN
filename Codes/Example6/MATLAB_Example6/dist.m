function distance = dist(x1, y1, x2, y2)
    % Calculate the Euclidean distance between points (x1, y1) and (x2, y2)
    % x1, y1, x2, y2 can be arrays, and the output will be an array
    distance = sqrt((x2 - x1).^2 + (y2 - y1).^2);
end