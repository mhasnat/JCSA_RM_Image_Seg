function eval = getRAGMergeSegments_Edge_Plane_gui(combFeat, tmpLabImg, thOptions, parNormal, parColor, rgnb, allInfoSI, scale)

rgnbT = rgnb; % Region information
parNormalT = parNormal; % Parameters for each region
parColorT = parColor; % Color Parameters for each region

tNS = length(rgnbT.nbsegs); % total number of segments
segmentConsidered = zeros(1,tNS);
newLabels = tmpLabImg;
[numR, numC] = size(newLabels);
STATS = regionprops(tmpLabImg, 'PixelIdxList');

% All thresholds
thDivNormalMax = thOptions.thDivNormalMax;
thDivNormalMin = thOptions.thDivNormalMin;
planarityTh = thOptions.planarityTh;
thKappa = thOptions.thKappa;
edgeStrengthTh = thOptions.edgeStrengthTh;

% Collect necessary information
Igf = getRGBDGradient(double(allInfoSI.rgbImg) / 255, allInfoSI.depImg / max(allInfoSI.depImg(:)));
magEdgeSC = Igf(1:scale:end, 1:scale:end);
magEdgeSC = magEdgeSC ./ max(magEdgeSC(:));

% magEdgeSC = allInfoSI.magGrdCD(1:2:end, 1:2:end);
tNS = length(rgnbT.nbsegs); % total number of segments

for segNum = 1:tNS
    if(segmentConsidered(segNum)) continue; end
    
    while(~segmentConsidered(segNum))
        
        % Check the eligibility of region merging based on Kappa value (>thKappa)
        if(parNormal.kappa(segNum)<thKappa)
            segmentConsidered(segNum) = 1;
            continue;
        end
        
        acind = rgnbT.nbsegs{segNum};
        Div_N = ones(1, length(acind))*500;
        
        for ti=1:length(acind)
            if(parNormal.kappa(acind(ti))>thKappa)
                % Div_N(ti) = (rgnbT.Div_N(segNum, acind(ti)) + rgnbT.Div_N(acind(ti), segNum))/2; % average
                Div_N(ti) = min(rgnbT.Div_N(segNum, acind(ti)), rgnbT.Div_N(acind(ti), segNum)); % minimum
                % Div_N(ti) = max(rgnbT.Div_N(segNum, acind(ti)), rgnbT.Div_N(acind(ti), segNum)); % maximum
            end
        end
        
        % Check -1 the validity to merge based on threshold BD_max
        thChk1 = Div_N<=thDivNormalMax;
        
        if(sum(thChk1)==0)
            segmentConsidered(segNum) = 1;
            continue;
        end
        
        for cn=1:length(thChk1)
            % Check - 3 : Edge
            if(thChk1(cn))
                % NEW * Edge aware connectivity
                % Check the strength of the edge among the two regions
                tedgesIndices = rgnbT.adjEdges{segNum,acind(cn)}; % add edge pixels
                if(~isempty(tedgesIndices))
                    tedgesIndices_main = sub2ind(size(newLabels), tedgesIndices(:,1), tedgesIndices(:,2));
                    
                    tedgesIndices(tedgesIndices(:,1)==1 | tedgesIndices(:,1)==numR, :) = [];
                    tedgesIndices(tedgesIndices(:,2)==1 | tedgesIndices(:,2)==numC, :) = [];
                    
                    tedgesIndices_clipped = sub2ind(size(newLabels), tedgesIndices(:,1), tedgesIndices(:,2));
                    
                    % Modification
                    tedgesIndices_allD = [];
                    tedgesIndices_allD(:,1) = tedgesIndices_clipped + numR; % columnwise forward
                    tedgesIndices_allD(:,2) = tedgesIndices_clipped - numR; % columnwise backward
                    tedgesIndices_allD(:,3) = tedgesIndices_clipped + 1; % rowwise forward
                    tedgesIndices_allD(:,4) = tedgesIndices_clipped - 1; % rowwise backward
                    tedgesIndices_allD(:,5) = tedgesIndices_allD(:,1) + 1; % lower diagonal forward
                    tedgesIndices_allD(:,6) = tedgesIndices_allD(:,1) - 1; % upper diagonal forward
                    tedgesIndices_allD(:,7) = tedgesIndices_allD(:,2) + 1; % lower diagonal backward
                    tedgesIndices_allD(:,8) = tedgesIndices_allD(:,2) - 1; % upper diagonal backward
                    
                    uniqueIndices = unique([tedgesIndices_main; tedgesIndices_allD(:)]);
                else
                    uniqueIndices = [];
                    tedgesIndices_main = [];
                end
                
                % Apply edge aware thresholding
                tedgeGbp = mean(magEdgeSC(uniqueIndices))*3;
                % tedgeGbp = mean(magEdgeSC(tedgesIndices_main));
                if(tedgeGbp > edgeStrengthTh)
                    thChk1(cn) = 0; continue;
                end
            end
        end
        
        thChk = find(thChk1==1);
        
        if(isempty(thChk))
            segmentConsidered(segNum) = 1;
            continue;
        end
        
        %% >> FOr speed up
        labelRep = [];
        [~, aa] = sort(Div_N(thChk));
        
        for tcnt = 1:length(aa)
            tactInd = thChk(aa(tcnt));
            % Check - 2 : Planarity
            if(Div_N(tactInd)>thDivNormalMin)
                % Check the planar condition - it is remains a plane after merging
                % these two regions
                outlierRatio = getPlanarityScore([segNum, acind(tactInd)], STATS, combFeat(:, 4:6));
                if(outlierRatio<planarityTh)
                    continue;
                else
                    labelRep = acind(thChk(aa(tcnt)));
                    break;
                end
            else
                labelRep = acind(thChk(aa(tcnt)));
                break;
            end
        end
        
        if(isempty(labelRep))
            segmentConsidered(segNum) = 1;
            continue;
        end
        
        % Correct the labels in the image
        assocPixelList = STATS(labelRep).PixelIdxList;
        assocPixelList = [assocPixelList; tedgesIndices_main]; % update pixel list
        
        rgnbT.adjEdges{segNum,labelRep} = []; % remove the edge among them
        
        newLabels(assocPixelList) = segNum;
        STATS(segNum).PixelIdxList = [STATS(segNum).PixelIdxList; assocPixelList];
        
        % Update parameters and divergences
        parNormalT = updateMergedRegionParameters(parNormalT, parNormalT, [segNum labelRep], segNum); % Normal
        parColorT = updateMergedRegionParametersColor(parColorT, parColorT, [segNum labelRep], segNum); % Color
        
        % Correct the neighbors and boundary information
        adjsToAdd = rgnbT.nbsegs{labelRep};
        adjsToAdd(find(adjsToAdd==segNum)) = []; % remove the segment itself from list
        
        thisSegAdj = rgnbT.nbsegs{segNum};
        
        for caj = 1:length(adjsToAdd)
            tcadj = adjsToAdd(caj);
            
            if(tcadj==31)
                kk = 1;
            end
            
            adjExist = find(thisSegAdj == tcadj);
            if(isempty(adjExist))
                % Add node -: Means this is a new neighbor
                rgnbT.nbsegs{segNum}(end+1) = tcadj;
                rgnbT.nbsegs{tcadj}(end+1) = segNum;
                
                % Add edges of this new neighbor
                rgnbT.adjEdges{segNum, tcadj} = rgnbT.adjEdges{labelRep,tcadj};
                rgnbT.adjEdges{tcadj,segNum} = rgnbT.adjEdges{tcadj,labelRep};
                
                % Erase edgelists of the merged segment
                rgnbT.adjEdges{labelRep, tcadj} = [];
                rgnbT.adjEdges{tcadj, labelRep} = [];
                
                % Remove the merged segment from the neighbors of others
                ind2remove = find(rgnbT.nbsegs{tcadj}==labelRep);
                rgnbT.nbsegs{tcadj}(ind2remove) = [];
                
                % Compute BD_normal among the neighbors
                rgnbT.Div_N(segNum, tcadj) = compute_Div_N(parNormalT, segNum, tcadj);
                rgnbT.Div_N(tcadj, segNum) = compute_Div_N(parNormalT, tcadj, segNum);
                rgnbT.KappaMerged(segNum, tcadj) = mergeKappaWMM(parNormalT, [segNum, tcadj]);
                rgnbT.KappaMerged(tcadj, segNum) = rgnbT.KappaMerged(segNum, tcadj);
                
                % Compute BD_color among the neighbors
                rgnbT.Div_C(segNum, tcadj) = compute_Div_C(parColorT, segNum, tcadj);
                rgnbT.Div_C(tcadj, segNum) = compute_Div_C(parColorT, tcadj, segNum);
            else
                % Update node -: Means this neighbor already exists in the list
                % Therefore it is necessary to update the edges
                tmedges = [rgnbT.adjEdges{segNum, tcadj}; rgnbT.adjEdges{labelRep, tcadj}];
                rgnbT.adjEdges{segNum, tcadj} = tmedges;
                rgnbT.adjEdges{tcadj,segNum} = tmedges;
                
                % Erase edgelists of the merged segment
                rgnbT.adjEdges{labelRep, tcadj} = [];
                rgnbT.adjEdges{tcadj, labelRep} = [];
                
                % Remove the merged segment from the neighbors of others
                ind2remove = find(rgnbT.nbsegs{tcadj}==labelRep);
                rgnbT.nbsegs{tcadj}(ind2remove) = [];
                
                % Update BD_normal among the neighbors
                rgnbT.Div_N(segNum, tcadj) = compute_Div_N(parNormalT, segNum, tcadj);
                rgnbT.Div_N(tcadj, segNum) = compute_Div_N(parNormalT, tcadj, segNum);
                rgnbT.KappaMerged(segNum, tcadj) = mergeKappaWMM(parNormalT, [segNum, tcadj]);
                rgnbT.KappaMerged(tcadj, segNum) = rgnbT.KappaMerged(segNum, tcadj);
                
                % Update BD_color among the neighbors
                rgnbT.Div_C(segNum, tcadj) = compute_Div_C(parColorT, segNum, tcadj);
                rgnbT.Div_C(tcadj, segNum) = compute_Div_C(parColorT, tcadj, segNum);
            end
        end
        segmentConsidered(labelRep) = -1;
        
        % Finally remove the segment from the list
        ind2remove = find(rgnbT.nbsegs{segNum}==labelRep);
        rgnbT.nbsegs{segNum}(ind2remove) = [];
    end
end

% eval.newLabels = removeSmallObjectsFinally(newLabels, parNormalT, rgnbT);
eval = removeSmallObjectsFinally(newLabels, parNormalT, rgnbT);