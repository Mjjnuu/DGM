function [currentUpper, nextLower, currentLower, previousUpper, coefficients] = centralFluxMatrices(N,P)
    currentUpper=kron(diag(ones(1,N)), ones(1,P+1));
    currentLower=kron(diag(ones(1,N)),  (-1).^(0:P));
    nextLower=kron(toeplitz(zeros(1,N),[0 1 zeros(1,N-2)] ),  (-1).^(0:P));
    previousUpper=kron(toeplitz([0 1 zeros(1,N-2)], zeros(1,N)), ones(1,P+1));
    coefficients=kron(ones(N,1), (-1).^(0:P)');
end