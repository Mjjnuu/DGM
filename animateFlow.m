function animateFlow(t,x,X,results,fps,multiplier,lower,upper)
    N=length(X)-1;
    P=size(results, 1)/N-1;
    startT=t(1);
    endT=t(end);
    nT=length(t);
    deltaT=(endT-startT)/(nT-1);
    frameT=max(floor(nT/(fps/multiplier*(endT-startT))),1);
    
    for i=1:frameT:nT
        tic
        U=reshape(results(:,i),P+1,N);
        Y=plotLP(x,X,U);
        plot(x,Y,'r');
        axis([x(1) x(end) lower upper]);
        text(1.02,x(end)+0.1,'Time');
        text(1.02,1,num2str((i-1)*(endT-startT)*1/nT));
        pause(max( 1/fps-toc,0 ));
    end
end    