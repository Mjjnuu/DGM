a=-1;
b=1;
n=1000;
h=(b-a)/n;
t=a:h:b;
y=zeros(5,n+1);

A=-1;
B=3;
ratio=(B-A)/(b-a);
T=t*ratio + A - a*ratio;

AA=0;
BB=0.5;
Ratio=(BB-AA)/(b-a);
TT=(t-a)*Ratio + AA;

Ksi=(2*TT-AA-BB)/(BB-AA);

plot (t,legendreP(4,t));
hold on

plot (T,legendreP(4,Ksi));
plot (TT,legendreP(5,Ksi));

% 
% for i=1:5
%     y(i,:)=legendreP(i,t);
% end    
% 
% plot (t,legendreP(1,t));
% hold on
% 
% for j=2:5
%     plot (t,y(j,:)) 
% end    

