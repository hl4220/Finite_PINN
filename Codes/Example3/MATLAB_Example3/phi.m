function R = phi(x, y, segments)
    % Calculate the R value for each pair in (x, y) arrays for an array of line segments
    % x, y can be arrays, segments is an N x 4 matrix
    m = 1;
    R = zeros(size(x));
    num_segments = size(segments, 1);
    
    for i = 1:num_segments
        phi_val = linseg(x, y, segments(i, 1), segments(i, 2), segments(i, 3), segments(i, 4));
        R = R + 1 ./ phi_val.^m;
    end
    
    R = 1 ./ R.^(1/m);
end
