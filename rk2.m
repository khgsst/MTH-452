function [tnp1,ynp1] = rk2(tn,yn,h,f)
%aka explicit midpoint
  f1 = feval(f,tn,yn);
yn1 = yn + 0.5*h*f1;
tnp1 = tn + 0.5*h;
f2 = feval(f,tnp1,yn1);
ynp1 = yn + h*f2;
