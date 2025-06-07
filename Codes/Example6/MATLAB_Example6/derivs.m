function [dU_dx,dU_dy] = derivs(nodes,elements,U)
N = size(nodes,1);
T = size(elements,1);
dU_dx = zeros(N,1);
dU_dy = zeros(N,1);
num = zeros(N,1);
for e=1:T
    nc = elements(e,:);
    Pe=[ones(3,1),nodes(nc,:)];
    Area=abs(det(Pe))/2; 
    C=inv(Pe); 
    grad=C(2:3,:);
    u = U(nc);
    du = grad * u;
    dU_dx(nc) = dU_dx(nc) +  du(1);
    dU_dy(nc) = dU_dy(nc) +  du(2);
    num(nc) = num(nc)+1;
end
dU_dx = dU_dx./num;
dU_dy = dU_dy./num;
end