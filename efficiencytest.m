%clear all;
t=-5:0.0001:5;

maxP=5;
maxN=2*4;

Results=zeros(maxP,maxN/2);
func = @(x) exp(sin(x)).*cos(2.*x)+sin(3.*x-1);
hold on


plot(t,func(t))

% for p=0:maxP
%     for n=2:2:maxN
%         h=10/n;
%         X=-5:h:5;
%         Summa=0;
%         for k=1:5
%             tic 
%             Multipliers=ApprLP(func,t,X,p);
%             Y=plotLP(t,X,Multipliers);
%             Summa = Summa + toc;
%         end 
%         
%         Results(p+1,n/2)=Summa/5; 
%     end 
%     
% end    
% 
% 
% for power=1:maxP+1
%     plot(Results(power,:))
% end
% 
% figure
% hold on
% for elements=1:maxN/2
%     plot(Results(:,elements))
% end    
% plot(1:maxP+1,0.0012*(1:maxP+1).^2.2 + 0.009)
% Results

X=-5:10/4:5;
p=4;

totaltime=0;
for i=1:100
    tic
    Multipliers=ApprLP(func,X,p);
    Y=plotLP(t,X,Multipliers);
    %plot(t,Y)
    totaltime = totaltime + toc;
end    
%Multipliers=ApprLP(func,t,X,p);
%Y=plotLP(t,X,Multipliers);
disp(totaltime)
plot(t,Y)


% U=1:1/743:5;
% length(U);
% 
% quadweights = [0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, ...
%     0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.0812743883615744];
% 
% quadpoints = [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0, ...
%     0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261];
% 
% 
% time=0;
% for i=1:1000000
%     tic;
%     sum(quadweights.*sin(quadpoints));
%     time = time + toc;
% end
% disp("Sum time: " + time);
% 
% qp=quadpoints';
% 
% Time=0;
% for i=1:1000000
%     tic;
%     quadweights*quadpoints';
%     Time = Time + toc;
% end
% disp("Matrix product time: " + Time);












