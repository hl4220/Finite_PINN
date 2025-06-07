function phi_val = linseg(x, y, x1, y1, x2, y2)
    % Calculate the phi value for each pair of points (x, y) and the line segment defined by (x1, y1) and (x2, y2)
    % x, y, x1, y1, x2, y2 can be arrays, and the output will be an array
    L = dist(x1, y1, x2, y2);
    xc = (x1 + x2) / 2;
    yc = (y1 + y2) / 2;
    f = (1 ./ L) .* ((x - x1) .* (y2 - y1) - (y - y1) .* (x2 - x1));
    t = (1 ./ L) .* ((L / 2).^2 - dist(x, y, xc, yc).^2);
    varphi = sqrt(t.^2 + f.^4);
    phi_val = sqrt(f.^2 + (1/4) .* (varphi - t).^2);
end
