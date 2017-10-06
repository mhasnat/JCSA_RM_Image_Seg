function newLabelsT = removeSmallObjectsFinally(tmpLabImg, parNormal, rgnb)

newLabelsT = tmpLabImg;
parNormalT = parNormal;
rgnbT = rgnb;

STATS = regionprops(newLabelsT, 'Area', 'PixelIdxList');
totalNumRegions = length(STATS);

totalNumPix = length(newLabelsT(:));

thNumPix = 100;
thPrb = thNumPix / totalNumPix;
segmentConsidered = zeros(1,length(STATS));

for t_segNum=1:totalNumRegions
    
    if(t_segNum==98)
        ee = 1;
    end
    
    if(STATS(t_segNum).Area==0) continue; end;
        
    % prb = parNormalT.alpha(t_segNum);
    prb = STATS(t_segNum).Area / totalNumPix;
    
    % Check the eligibility of the segment to survive
    if(prb < thPrb)
        % Remove the segment and assign it to the nearest one based on
        % minimum number of edge
        t_nnseg = rgnbT.nbsegs{t_segNum};
        
        t_nnseg(find(t_nnseg>totalNumRegions)) = [];
        nnsegstat = find(segmentConsidered(t_nnseg)==-1); % find the small segments which do not exist anymore
        t_nnseg(nnsegstat) = [];
        
        numEdgeNS = zeros(1, length(t_nnseg));
        
        for ti=1:length(t_nnseg)
            numEdgeNS(ti) = size(rgnbT.adjEdges{t_segNum,t_nnseg(ti)},1);
            % numEdgeNS(ti) = compute_Div_N(parNormalT, t_segNum, t_nnseg(ti));
        end
        
        % Decide which segment to be merged with - now based on the segment
        % with which it shares most number of edges
        [~, taa] = min(numEdgeNS);
        t_segToMerge = t_nnseg(taa);
        if(t_segToMerge>totalNumRegions) continue; end
        
        %% Update all necessary information
        tedgesIndices = rgnbT.adjEdges{t_segToMerge,t_segNum}; % add edge pixels
        if(~isempty(tedgesIndices))
            tedgesIndices = sub2ind(size(newLabelsT), tedgesIndices(:,1), tedgesIndices(:,2));
        end
        
        % Correct the labels in the image
        assocPixelList = STATS(t_segNum).PixelIdxList;
        assocPixelList = [assocPixelList; tedgesIndices]; % update pixel list
        
        rgnbT.adjEdges{t_segNum, t_segToMerge} = []; % remove the edge among them
        rgnbT.adjEdges{t_segToMerge, t_segNum} = []; % remove the edge among them
                
        newLabelsT(assocPixelList) = t_segToMerge;
        STATS(t_segToMerge).PixelIdxList = [STATS(t_segToMerge).PixelIdxList; assocPixelList];
        
        % Update parameters and divergences
        parNormalT = updateMergedRegionParameters(parNormalT, parNormalT, [t_segNum t_segToMerge], t_segToMerge);
        
        % Correct the neighbors and boundary information
        adjsToAdd = rgnbT.nbsegs{t_segNum};
        adjsToAdd(find(adjsToAdd==t_segToMerge)) = []; % remove the segment itself from list
        
        thisSegAdj = rgnbT.nbsegs{t_segToMerge};
        
        for caj = 1:length(adjsToAdd)
            tcadj = adjsToAdd(caj);
            
            adjExist = find(thisSegAdj == tcadj);
            if(isempty(adjExist))
                % Add node -: Means this is a new neighbor
                rgnbT.nbsegs{t_segToMerge}(end+1) = tcadj;
                rgnbT.nbsegs{tcadj}(end+1) = t_segToMerge;
                
                % Add edges of this new neighbor
                rgnbT.adjEdges{t_segToMerge, tcadj} = rgnbT.adjEdges{t_segNum,tcadj};
                rgnbT.adjEdges{tcadj,t_segToMerge} = rgnbT.adjEdges{tcadj,t_segNum};
                
                % Erase edgelists of the merged segment
                rgnbT.adjEdges{t_segNum, tcadj} = [];
                rgnbT.adjEdges{tcadj, t_segNum} = [];
                
                % Remove the merged segment from the neighbors of others
                ind2remove = find(rgnbT.nbsegs{tcadj}==t_segNum);
                rgnbT.nbsegs{tcadj}(ind2remove) = [];
                
                % % %                 % Compute BD_normal among the neighbors
                % % %                 rgnbT.Div_N(t_segToMerge, tcadj) = compute_Div_N(parNormalT, t_segToMerge, tcadj);
                % % %                 rgnbT.Div_N(tcadj, t_segToMerge) = rgnbT.Div_N(t_segToMerge, tcadj);
                % % %                 rgnbT.KappaMerged(t_segToMerge, tcadj) = mergeKappaWMM(parNormalT, [t_segToMerge, tcadj]);
                % % %                 rgnbT.KappaMerged(tcadj, t_segToMerge) = rgnbT.KappaMerged(t_segToMerge, tcadj);
            else
                % Update node -: Means this neighbor already exists in the list
                % Therefore it is necessary to update the edges
                tmedges = [rgnbT.adjEdges{t_segToMerge, tcadj}; rgnbT.adjEdges{t_segNum, tcadj}];
                rgnbT.adjEdges{t_segToMerge, tcadj} = tmedges;
                rgnbT.adjEdges{tcadj,t_segToMerge} = tmedges;
                
                % Erase edgelists of the merged segment
                rgnbT.adjEdges{t_segNum, tcadj} = [];
                rgnbT.adjEdges{tcadj, t_segNum} = [];
                
                % Remove the merged segment from the neighbors of others
                ind2remove = find(rgnbT.nbsegs{tcadj}==t_segNum);
                rgnbT.nbsegs{tcadj}(ind2remove) = [];
                
                % % %                 % Update BD_normal among the neighbors
                % % %                 rgnbT.Div_N(t_segToMerge, tcadj) = compute_Div_N(parNormalT, t_segToMerge, tcadj);
                % % %                 rgnbT.Div_N(tcadj, t_segToMerge) = rgnbT.Div_N(t_segToMerge, tcadj);
                % % %
                % % %                 rgnbT.KappaMerged(t_segToMerge, tcadj) = mergeKappaWMM(parNormalT, [t_segToMerge, tcadj]);
                % % %                 rgnbT.KappaMerged(tcadj, t_segToMerge) = rgnbT.KappaMerged(t_segToMerge, tcadj);
            end
        end
        
        segmentConsidered(t_segNum) = -1;
        
    end
end
