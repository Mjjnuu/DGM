function results=hyperbolicSolver(t,x,X,P,F,dF,S,bc,method,flux)
    %Input:
    %t, the time domain vector. Must satisfy t(i)-t(i-1)=t(j)-t(j-1)
    %for all i and j.
    %x, the spatial domain vector
    %X, partition of spatial domain
    %P, maximum degree of polynomial space
    %F, force function
    %dF, derivative of F
    %S, source function
    %bc, vector containing boundary conditions u(x(1),t)
    %length(bc)=length(t)
    
    
    %Solves the hyperbolic partial differential equation du/dt+dF(u)/dx=S(u) 
    %using discontinuous Galerkin method

    
    a=x(1);
    b=x(end);
    
    N=length(X)-1;
    startT=t(1);
    endT=t(end);
    nT=length(t);
    deltaT=(endT-startT)/nT;
    
    
    iM=diag(kron(ones(1,N),(2*(0:P)+1)/2));
    iJ=diag(kron(2./diff(X),ones(1,P+1)));
    
    [currentUpper, nextLower, currentLower, previousUpper, coefficients] = centralFluxMatrices(N,P);
    
    boundaryCondition=zeros(N*(P+1),1);
    boundaryCondition(1:P+1)=(-1).^(0:P);
    
    results=zeros(N*(P+1),nT);

    [qW,qP]=GaussLegendreQuad(21);
    LP=LegendreMatrix(qP);
    LPder=LegendreDerivativeMatrix(qP);
    

    switch method
        case 'Euler'
            switch flux
                case 'upwind'
                    for i=1:nT
                        W=results(:,i);
                        results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,W,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(upwindFlux(N,W,F) - F(bc(i,deltaT)).*boundaryCondition));
                        if i/floor(nT/10)==floor(i/floor(nT/10))
                            fprintf('%s','*');
                        end
                    end  %upwind Euler  
                    fprintf('%s\n', '');

                case 'Godunov'   
                    for i=1:nT
                        W=results(:,i);
                        results(:,i+1)=W + deltaT.*(iM*generalIntegral(N,W,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,W,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,W,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)...
                            ));
                        
                        if i/floor(nT/10)==floor(i/floor(nT/10))
                            fprintf('%s','*');
                        end
                    end  %Euler Godunov    
                    fprintf('%s\n', '');
                    
                otherwise
                    error("Wrong type of flux")
            end        
        
        case 'RK4'
        
            switch flux
            
                case 'upwind'
                    for i=1:nT
                        
                        U=results(:,i);
                        
                        k1 = deltaT.*(iM*generalIntegral(N,U,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,U,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(upwindFlux(N,U,F) - F(bc(i,deltaT)).*boundaryCondition));
                        
                        k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,U+k1/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(upwindFlux(N,U+k1/2,F) - F(bc(i,deltaT)).*boundaryCondition));
                        
                        k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,S,qW,qP,LP,LP) ...
                            + iM*iJ* generalIntegral(N,U+k2/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(upwindFlux(N,U+k2/2,F) - F(bc(i,deltaT)).*boundaryCondition));
                        
                        k4 = deltaT.*(iM*generalIntegral(N,U+k3,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,U+k3,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(upwindFlux(N,U+k3,F) - F(bc(i,deltaT)).*boundaryCondition));
                        
                        results(:,i+1) = U +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
                        
                        
                        if i/floor(nT/10)==floor(i/floor(nT/10))
                            fprintf('%s','*');
                        end
                        
                    end  %Runge-Kutta upwind
                    fprintf('%s\n', '');
                    
                case 'Godunov'
                    
                    for i=1:nT
                        
                        U=results(:,i);
                        
                        k1 = deltaT.*(iM*generalIntegral(N,U,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,U,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,U,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k2 = deltaT.*(iM*generalIntegral(N,U+k1/2,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,U+k1/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,U+k1/2,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k3 = deltaT.*(iM*generalIntegral(N,U+k2/2,S,qW,qP,LP,LP) ...
                            + iM*iJ* generalIntegral(N,U+k2/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,U+k2/2,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k4 = deltaT.*(iM*generalIntegral(N,U+k3,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,U+k3,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,U+k3,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                    
                        results(:,i+1) = U +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
                        
                        if i/floor(nT/10)==floor(i/floor(nT/10))
                            fprintf('%s','*');
                        end
                    end   %Runge-Kutta 4 Godunov
                    
                    fprintf('%s\n', '');
                    
                otherwise
                    error("Wrong type of flux")
                    
            end
            
            
        case 'AB4' %%%%%%%%%4th order Adams-Bashforth
            
            switch flux
                case 'upwind'
                    
                    
                case 'Godunov'
                    Memory=zeros(N*(P+1),4);
                    for i=1:3
                        
                        R=results(:,i);
                        k1 = deltaT.*(iM*generalIntegral(N,R,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,R,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,R,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k2 = deltaT.*(iM*generalIntegral(N,R+k1/2,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,R+k1/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,R+k1/2,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k3 = deltaT.*(iM*generalIntegral(N,R+k2/2,S,qW,qP,LP,LP) ...
                            + iM*iJ* generalIntegral(N,R+k2/2,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,R+k2/2,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        
                        k4 = deltaT.*(iM*generalIntegral(N,R+k3,S,qW,qP,LP,LP) ...
                            +  iM*iJ*generalIntegral(N,R+k3,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,R+k3,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        results(:,i+1) = R +(1/6)*(k1 + 2*k2 + 2*k3 + k4);
                    end
%                     Memory(:,2)=results(:,1);
%                     Memory(:,3)=results(:,2);
%                     Memory(:,4)=results(:,3);
                    
                    Memory(:,2)=deltaT.*(iM*generalIntegral(N,results(:,1),S,qW,qP,LP,LP) + ...
                        iM*iJ*generalIntegral(N,results(:,1),F,qW,qP,LP,LPder) ...
                        - iM*iJ*(GodunovFlux(N,results(:,1),F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                    Memory(:,3)=deltaT.*(iM*generalIntegral(N,results(:,2),S,qW,qP,LP,LP) + ...
                        iM*iJ*generalIntegral(N,results(:,2),F,qW,qP,LP,LPder) ...
                        - iM*iJ*(GodunovFlux(N,results(:,2),F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                    Memory(:,4)=deltaT.*(iM*generalIntegral(N,results(:,3),S,qW,qP,LP,LP) + ...
                        iM*iJ*generalIntegral(N,results(:,3),F,qW,qP,LP,LPder) ...
                        - iM*iJ*(GodunovFlux(N,results(:,3),F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                    for i=4:nT
                        R=results(:,i);
                        Memory(:,1)=Memory(:,2);
                        Memory(:,2)=Memory(:,3);
                        Memory(:,3)=Memory(:,4);
                        Memory(:,4)=deltaT.*(iM*generalIntegral(N,R,S,qW,qP,LP,LP) + ...
                            iM*iJ*generalIntegral(N,R,F,qW,qP,LP,LPder) ...
                            - iM*iJ*(GodunovFlux(N,R,F,bc(i,deltaT), currentUpper, nextLower, currentLower, previousUpper)));
                        results(:,i+1)=R+55/24*Memory(:,4)-59/24*Memory(:,3)+37/24*Memory(:,2)-9/24*Memory(:,1);
                        if i/floor(nT/10)==floor(i/floor(nT/10))
                            fprintf('%s','*');
                        end
                    end
                    fprintf('%s\n', '');
                    %disp("Adams-Bashfort: " + toc);
                    
                otherwise
                    error("Wrong type of flux")
                
            end %flux Adams-Bashforth    
         
        otherwise    
            error("Wrong type of method")
            
    end %switch       
            
end %function
                    
                    
                    
                    
                    
                    
                    
                    
                    