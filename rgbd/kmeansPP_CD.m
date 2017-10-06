function [L,C] = kmeansPP_CD(X,k, maxNum)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 07/08/11
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

if(nargin<3)
    maxNum = 200;
end

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    D = zeros(1, size(X,2));
    for i = 2:k
        for ii = 1:i-1
            indx = find(L==ii);
            tdist = 1 - X'*C(:,ii).^2;
            D(indx) = tdist(indx);
        end
                
        D = cumsum(D);
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(C'*X);
    end
    
    numIter = 1;
    while any(L ~= L1)
        L1 = L;
        
        for i = 1:k
            l = L==i; 
            %             [~,~,C(:,i)] = svds(X(:,l)',1);
            C(:,i) = sum(X(:,l),2)/sum(l);
            C(:,i) = C(:,i) ./ norm(C(:,i));
        end
        
        [~,L] = max((C'*X).^2);
        
        numIter = numIter + 1;
        if(numIter>maxNum)
            break;
        end
    end
    
end
