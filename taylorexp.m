function r=taylorexp(n,z)
  r=1+z;
  fac=1;
  for i=2:n,
    fac=fac*i;
    r=r+z.^i/fac;
  end
  