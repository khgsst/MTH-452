function [tnp1,ynp1] = BDF2si(tn,yn,ynm1,h,f)
  %aka trapezoidal with simple fixed point iteration
f1 = feval(f,tn,yn);
f2 = feval(f,tn-h,ynm1);
p = yn + h*(3*f1-f2)/2;
tnp1 = tn + h;
for i = 1:10
    f3 = feval(f,tnp1,p);
    c = (4*yn - ynm1 + 2*h*f3)/3;
    if norm(c - p) < 1e-3*norm(c)
        ynp1 = c;
        return
    end
    p = c;
end

error('Simple iteration failed to converge.')
