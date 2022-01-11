function Rabs=myplotBL(alpha,beta,axisbox)
%
% plot the boundary locus for the region of absolute stability of a Linear
% Multistep Method with characteristic polynomials rho and sigma.
%
% axisbox = axis region 
%           (if omitted, default values are used that may be a poor choice)
%
% See also mymakeplotBL.m, from which this routine is called for a variety of
% specific LMMs.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter7  (2007)
%
% modified mar 2010 to add axis equal, shade internal region
  
if nargin<3
  axisbox=[];
end

plotSurf=1; % turn on surface shading

rho = @(z) polyval(alpha,z);
sigma = @(z) polyval(beta,z);

theta = linspace(0, 2*pi, 1001);
eitheta = exp(1i * theta);
z = rho(eitheta) ./ sigma(eitheta);
if (isempty(axisbox))
  xa=min(0,min(real(z)))-1; xb=max(0,max(real(z)))+1;
  ya=min(0,min(imag(z)))-1; yb=max(0,max(imag(z)))+1;
else
  xa = axisbox(1); xb = axisbox(2); ya = axisbox(3); yb = axisbox(4);
end

figure
if (plotSurf)
  %plot surf first because overlays white on outside of region
  nptsx = 51;
  nptsy = 51;
  x = linspace(xa,xb,nptsx);
  y = linspace(ya,yb,nptsy);
  [X,Y] = meshgrid(x,y);
  Z = X + 1i*Y;
  for i=1:nptsx
    for j=1:nptsy
      Rabs(i,j)=norm(roots(alpha-Z(i,j)*beta),inf);
    end
  end
  colormap([.7 .7 1;1 1 1])
  %colormap([.7 .7 .7;1 1 1])   % grayscale
  caxis([.9 1.1])
  contourf(x,y,Rabs,[0 1],'k')
end

% plot axes:
hold on
plot([xa xb],[0 0],'k')
plot([0 0],[ya yb],'k')
box on

axis equal 
axis([xa xb ya yb]) 
set(gca,'FontSize',15)

% highlight boundary locus
plot(z,'LineWidth', 4)

