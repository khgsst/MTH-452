function [tnp1,ynp1] = AM2ni(tn,yn,h,f,dfdy)
  %aka trapezoidal with newton iteration
f1 = feval(f,tn,yn);
p = yn + h*f1;
tnp1 = tn + h;
J = feval(dfdy,tn,yn);
[L,U] = lu(eye(size(J)) - 0.5*h*J);
for i = 1:10
    f2 = feval(f,tnp1,p);
    c = yn + 0.5*h*(f1 + f2);
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
