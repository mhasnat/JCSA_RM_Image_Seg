function rgnb = getRegionNeighbors(tmpLabImg)

[~, e] = imRAG(tmpLabImg);

regLabels = unique(e);

for i=1:length(unique(regLabels))
    segNum = regLabels(i);

    ind1 = find(e(:,1)==segNum);
    acind = e(ind1,2);
    
    ind2 = find(e(:,2)==segNum);
    acind = [acind; e(ind2,1)];
    
    rgnb{segNum} = acind;
end