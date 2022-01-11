function [tnp1,ynp1] = bdf2ni(tn,yn,ynm1,h,f,dfdy)
  %aka trapezoidal with newton iteration
f1 = feval(f,tn,yn);
f2 = feval(f,tn-h,ynm1);
p = yn + h*(3*f1-f2)/2;
tnp1 = tn + h;
J = feval(dfdy,tn,yn);
[L,U] = lu(eye(size(J)) - 2/3*h*J);
for i = 1:10
    f3 = feval(f,tnp1,p);
    c = (4*yn - ynm1 + 2*h*f3)/3;
    residual = c - p;
    step = U \ ( L \residual);
    cn = p + step;
    if norm(step) < 1e-3*norm(cn)
        ynp1 = cn;
        return
    end
    p = cn;
end
error('Newton iteration failed to converge.')
