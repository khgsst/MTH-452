function fderror

  methods = {'Dp';'Dm'};
  plotstr = {'-b+';'--rx'};

  u = @(x) sin(x)
  x = 1

  powers = -1:-1:-8;
  h = 2.^powers;

  figure;
  for j=1:length(methods)
    method = methods{j};
    disp(['%%% ', method, ' %%%'])
    [approx,difference,ratio] = tableDf(method,u,x,h);

    fprintf('h  \tapprox       \terror     \tratio\n')
    for i=1:length(approx)
      fprintf('2^%d \t%0.12f \t%0.5e \t%g\n',powers(i),approx(i),difference(i),ratio(i))
    end
  
    loglog(h(2:end),difference(2:end),plotstr{j})
    hold on
  end
  axis tight
  xlabel('h')
  ylabel('error')
  legend(methods,'Location','Southeast');

  
function [approx,difference,ratio]=tableDf(method,u,x,h)
  
  n = length(h);
  approx = zeros(n,1);
  difference = zeros(n,1);
  ratio = zeros(n,1);
  
  for i=1:n
    hi = h(i);
    approx(i) = feval(method,u,x,hi);
  end
  difference(2:n) = abs(approx(2:n)-approx(1:n-1));
  ratio(3:n) = difference(2:n-1)./difference(3:n);

  
function du=Dp(u,x,h)
  du = (feval(u,x+h)-feval(u,x))/h;

  
function du=Dm(u,x,h)
  du = (feval(u,x)-feval(u,x-h))/h;

% third order accurate approximation given in (1.4)
%function du=D3(u,x,h)
%  du = (2*feval(u,x+h)+3*feval(u,x)-6*feval(u,x-h)+feval(u,x-2*h))/(6*h);
