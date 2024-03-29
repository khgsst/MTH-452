function mymakeplotS(method)
% makeplotS.m
%
% call plotS to plot the region of absolute stability S 
% (and the Order Star if plotOS=1 below) for a specific method.   
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter7  (2007)
%
% modified mar 2010 to add nystrom, m-simpson, bdf, taylorexp

if nargin < 1
  method = 'taylor';     % set to desired method, see choices below
end
plotOS = 0;           % set to 1 to also plot Order Star

switch method

  case 'taylor'
     % Taylor series method of order n
     n = 5;
     disp(['  order = ' num2str(n)])
     R = @(z) taylorexp(n,z);
     axisbox = [-5 5 -5 5];
     if n>7
       xy = 0.7*n;
       axisbox = [-xy xy -xy xy]; 
     end

  case 'trbdf2'
     R = @(z) (1 + 5/12 * z) ./ (1 - 7/12 * z + 1/12 * z.^2);
     axisbox = [-5 15 -10 10];

  case 'pade12'
     R = @(z) (1 + 1/3 * z) ./ (1 - 2/3 * z + 1/6 * z.^2);
     axisbox = [-2 8 -5 5];

  case 'pade22'
     R = @(z) (1 + 1/2 * z + 1/12 * z.^2) ./ (1 - 1/2 *z + 1/12 * z.^2);
     axisbox = [-5 5 -5 5];

  case 'pade13'
     R = @(z) (1 + 1/4 * z) ./ (1 - 3/4 * z + 1/4 * z.^2 - 1/24 * z.^3);
     axisbox = [-3 7 -5 5];

  case 'sdirk3'
     gamma = 2/3;
     R = @(z) (1 + z*(1-2*gamma) + z.^2 * (1/2 - 2*gamma + gamma^2)) ./ ...
              (1 - gamma*z).^2;
     axisbox = [-5 35 -20 20];

  case 'radau5'
     R = @(z) (1 + 2/5 * z + 1/20 * z.^2) ./ ...
              (1 - 3/5 * z + 3/20 * z.^2 - 1/60 * z.^3);
     axisbox = [-5 15 -10 10];

  otherwise
     disp(' *** method not recognized')
  end


% Now that R is set, plot S:
myplotS(R,plotOS,axisbox)

function r=taylorexp(n,z)
  r=1+z;
  fac=1;
  for i=2:n
    fac=fac*i;
    r=r+z.^i/fac;
  end
 