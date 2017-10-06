function indx = getInvalidIndices(Params, mmtype)

if(nargin<2)
    mmtype = 'gmm';
end

indx = [];

if(strcmp(mmtype, 'gmm'))
    for i=1:length(Params.alpha)
        % check alpha
        if(Params.alpha(i)==0)
            indx = [indx i];
            continue;
        end
        
        % check mu
        mu = Params.mu(i,:);
        
        if(isnan(mu))
            indx = [indx i];
            continue;
        end
        
        % check sigma
        Sigma = Params.sigma(:,:,i);
        
        [U,p]= chol(Sigma);
        
        if p ~= 0
            indx = [indx i];
        end
    end
else
    for i=1:length(Params.alpha)
        % check alpha
        if(Params.alpha(i)==0)
            indx = [indx i];
            continue;
        end
        
        % check mu
        mu = Params.mu(i,:);
        
        if(isnan(mu))
            indx = [indx i];
            continue;
        end
    end
end