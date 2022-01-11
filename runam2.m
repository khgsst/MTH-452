h=0.01;
T=pi;
t=[0:h:pi+h];
N=length(t);
yam=0*t;
lambda=-500;
f=@(t,y) lambda*(y-sin(t))+cos(t);
df=@(t,y) lambda;
fexact=@(t) exp(lambda*t)+sin(t);
yn=1;
yam(1)=yn;
tnp1=t(2);
for j=2:N
  %[tnp1,ynp1]=am2si(tnp1,yn,h,f);
  [tnp1,ynp1]=am2ni(tnp1,yn,h,f,df);
  yam(j)=ynp1;
  yn=ynp1;
end
figure
plot(t,yam,'-o')
hold on
plot(t,fexact(t),'--*r')
