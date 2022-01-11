h=0.01;
T=pi;
t=[0:h:pi+h];
N=length(t);
ybdf=0*t;
lambda=-100;
f=@(t,y) lambda*(y-sin(t))+cos(t);
df=@(t,y) lambda;
fexact=@(t) exp(lambda*t)+sin(t);
ynm1=1;
ynm1ab2=ynm1;
ybdf(1)=ynm1;
yn=fexact(h);
ynab2=yn;
ybdf(2)=yn;
yab2(2)=yn;
tnp1=t(2);
for j=3:N

  ynp1ab2=ynab2+h*(3*f(tnp1,ynab2)-f(tnp1-h,ynm1ab2))/2;
  yab2(j)=ynp1ab2;
  ynm1ab2=ynab2;
  ynab2=ynp1ab2;

  %[tnp1,ynp1]=bdf2si(tnp1,yn,ynm1,h,f);
  [tnp1,ynp1]=bdf2ni(tnp1,yn,ynm1,h,f,df);
  ybdf(j)=ynp1;
  ynm1=yn;
  yn=ynp1;

end
figure
plot(t,ybdf,'-o')
hold on
plot(t,fexact(t),'--*r')
plot(t,yab2,'-sg')
title(['\lambda =',num2str(lambda)])
legend('bdf2si','exact','ab2')