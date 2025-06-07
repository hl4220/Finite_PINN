function Ke = LKe(pk)
  Pe=[ones(3,1),pk]; % 3 by 3 matrix with rows=[1 xcorner ycorner] 
  Area=abs(det(Pe))/2; % area of triangle e = half of parallelogram area
  C=inv(Pe); % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
  grad=C(2:3,:);
  Ke=Area*grad'*grad; % element matrix from slopes b,c in grad
end