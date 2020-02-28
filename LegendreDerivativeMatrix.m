function LPder=LegendreDerivativeMatrix(qP)
   LPder=zeros(11,length(qP));    
   LPder(1,:) = 0;
   LPder(2,:) = 1;
   LPder(3,:) = 3*qP;
   LPder(4,:) = 1/2*(15*qP.^2-3);
   LPder(5,:) = 1/8*(140*qP.^3-60*qP);
   LPder(6,:) = 1/8*(315*qP.^4-210*qP.^2+15);
   LPder(7,:) = 21/8*(33*qP.^5-30*qP.^3+5*qP);
   LPder(8,:) = 7/16*(429*qP.^6-495*qP.^4+135*qP.^2-5);
   LPder(9,:) = 9/16*(715*qP.^7-1001*qP.^5+385*qP.^3-35*qP);
   LPder(10,:) = 45/128*(2431*qP.^8-4004*qP.^6+2002*qP.^4-308*qP.^2+7);
   LPder(11,:) = 55/128*(4199*qP.^9-7956*qP.^7+4914*qP.^5-1092*qP.^3+63);

end