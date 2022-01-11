function [tnp1,ynp1] = Heun(tn,yn,h,f)
%aka explicit trapezoidal
  f1 = feval(f,tn,yn);
yn1 = yn + h*f1;
tnp1 = tn + h;
f2 = feval(f,tnp1,yn1);
ynp1 = yn + 0.5*h*(f1 + f2);
