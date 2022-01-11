function mymakeplotBL(method,r)
% mymakeplotBL.m
%
% call myplotBL to plot the boundary locus of the stability region S
% for a specific Linear Multistep Method.
%
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter7  (2007)
%
% modified mar 2010 to add nystrom, bdf

if nargin < 1
  method = 'AM';     % set to desired method or class, see choices below
end

if nargin < 2
  r = 1;             % for r-step classes of methods
end

axisbox = []; % empty by default (myplotBL will create)
nick='';

switch method
  
 case 'AB'
  %rho = @(z) (z-1) .* z.^(r-1);
  alpha=zeros(1,r+1);
  alpha(1:2)=[1,-1];
  
  switch r
   case 1 
    %sigma = @(z) 1;
    beta=[0,1];
    nick='Forward Euler';
   case 2 
    %sigma = @(z) (3*z - 1)/2;
    beta=[0,3,-1]/2;
   case 3 
    %sigma = @(z) ((23*z - 16) .* z +5) / 12;
    beta=[0,23,-16,5]/12;
   case 4 
    %sigma = @(z) (((55*z - 59).*z +37).*z - 9) / 24;
    beta=[0,55,-59,37,-9]/24;
   case 5 
    %sigma = @(z) ((((1901*z - 2774).*z +2616).*z -1274).*z + 251) / 720;
    beta=[0,1901,-2774,2616,-1274,251]/720;
   case 6 
    %sigma = @(z) (((((4277*z - 7923).*z +9982).*z -7298).*z + 2877).*z - 475) / 1440;
    beta=[0,4277,-7923,9982,-7298,2877,-475]/1440;
   otherwise
    disp(' *** order not recognized')
    return
  end
  
 case 'AM'
  %rho = @(z) (z-1) .* z.^(r-1);
  alpha=zeros(1,r+1);
  alpha(1:2)=[1,-1];
  
  switch r
   case 1 
    axisbox = [-4 4 -4 4];
    %sigma = @(z) (z + 1)/2;
    beta=[1,1]/2;
    nick='Trapezoidal';
   case 2 
    %sigma = @(z) ((5*z + 8).*z -1) / 12;
    beta=[5,8,-1]/12;
   case 3 
    %sigma = @(z) (((9*z + 19).*z -5).*z +1) / 24;
    beta=[9,19,-5,1]/24;
   case 4 
    %sigma = @(z) ((((251*z + 646).*z - 264).*z + 106) -19).*z / 720;
    beta=[251,646,-264,106,-19]/720;
   case 5 
    %sigma = @(z) (((((475*z + 1427).*z - 798).*z + 482).*z ...
    %		   - 173).*z + 27)/ 1440;
    beta=[475,1427,-798,482,-173,27]/1440;
   otherwise
    disp(' *** order not recognized')
    return
  end
  
 case 'Nystrom'
  alpha=zeros(1,r+1);
  alpha(1:3)=[1,0,-1];
  
  switch r
   case 2 
    beta=[0,2,0];
   case 3 
    beta=[0,7,-2,1]/3;
   case 4 
    beta=[0,8,-5,4,-1]/3;
   otherwise
    disp(' *** order not recognized')
    return
  end

 case 'Milne-Simpson'
  alpha=zeros(1,r+1);
  alpha(1:3)=[1,0,-1];
  
  switch r
   case 2 
    beta=[1,4,1]/3;
   case 3 
    beta=[1,4,1,0]/3;
   otherwise
    disp(' *** order not recognized')
    return
  end

 case 'BDF'
  beta=zeros(1,r+1);
  switch r
   case 1
    beta(1)=1;
    alpha=[1,-1];
    nick='Backward Euler';    
   case 2
    beta(1)=2/3;
    alpha=[3,-4,1]/3;
   case 3
    beta(1)=6/11;
    alpha=[11,-18,9,-2]/11;
   case 4
    beta(1)=12/25;
    alpha=[25,-48,36,-16,3]/25;
   case 5
    beta(1)=60/137;
    alpha=[137,-300,300,-200,75,-12]/137;
   case 6
    beta(1)=60/147;
    alpha=[147,-360,450,-400,225,-72,10]/147;
   otherwise
    disp(' *** order not recognized')
    return
  end
  rho = @(z) polyval(alpha,z);
  sigma = @(z) polyval(beta,z);
  
 otherwise
  disp(' *** method not recognized')
  return
end


%z=plotBL(rho,sigma,axisbox);
myplotBL(alpha,beta,axisbox);

if (isempty(nick))
  title([num2str(r),'-step ',method,' method'])
else
  title([num2str(r),'-step ',method,' method (',nick,')'])
end
