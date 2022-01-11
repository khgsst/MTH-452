h=0.01;
T=pi;
t=[0:h:pi+h];
N=length(t);
ybdf=0*t;
lambda=-10;
f=@(t,y) lambda*(y-sin(t))+cos(t);
df=@(t,y) lambda;
fexact=@(t) exp(lambda*t)+sin(t);
ynm1=1;
ybdf(1)=ynm1;
yn=fexact(h);
ybdf(2)=yn;
tnp1=t(2);
for j=3:N,
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
