function [tnp1,ynp1] = AM2si(tn,yn,h,f)
  %aka trapezoidal with simple fixed point iteration
f1 = feval(f,tn,yn);
p = yn + h*f1;
tnp1 = tn + h;
for i = 1:10
    f2 = feval(f,tnp1,p);
    c = yn + 0.5*h*(f1 + f2);
    if norm(c - p) < 1e-3*norm(c)
        ynp1 = c;
        return
    end
    p = c;
end
error('Simple iteration failed to converge.')
