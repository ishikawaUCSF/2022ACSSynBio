function L = MeasureLine3D2(skelBW)

% Measure length of polyline in 3D. (08/01/2017)
% Input image should include a skeltonized polyline without branching using
% 'Skeleton3D'.
% This function calculate length from the end point step by step. 
% Need imaging tool box and functions 'Skel2Graph3D' and 'Neighbor26Matrix1'.
% Skel2Graph3D: https://github.com/phi-max/skel2graph3d-matlab


%%
pxSize = sum(skelBW(:));
[~,node,~] = Skel2Graph3D(skelBW,1);
numNode = length(node);

ep1 = cell2mat({node.ep});
endPointID = find(ep1);
branchID = find(ep1 == 0);
numEndPoint = length(endPointID);
numBranch = length(branchID);
sizeImg = size(skelBW);
comment = ['node: ', num2str(numNode), ', endopoint: ', num2str(numEndPoint), ' branch: ', num2str(numBranch)];
disp(comment); 
L = [0 0 0 0 0]; % x, y, z, idx, Euclidean distance


if numBranch == 0
    if numEndPoint == 2
        ep2 = endPointID(1);
        PxList1 = zeros(pxSize, 5);
        PxList1(1,1:4) = [round(node(1,ep2).comx), round(node(1,ep2).comy), round(node(1,ep2).comz), round(node(1,ep2).idx)];
        skelBW2 = skelBW;
        
        for n1 = 1:pxSize
            % crop neighbor image at 3x3
            skelBW2(PxList1(n1,1),PxList1(n1,2),PxList1(n1,3)) = 0;
            neighborIndex1 = find(skelBW2(PxList1(n1,1)-1:PxList1(n1,1)+1, PxList1(n1,2)-1:PxList1(n1,2)+1, PxList1(n1,3)-1:PxList1(n1,3)+1));
            [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

            ImgBW1 = false(sizeImg);
            ImgBW1(PxList1(n1,1)+x2-2, PxList1(n1,2)+y2-2, PxList1(n1,3)+z2-2) = 1;
            index1 = find(ImgBW1);
            sizInd1 = length(index1);
            if sizInd1 == 1
                [a2, b2, c2] = ind2sub(sizeImg, index1);
                PxList1(n1+1,1:4) = [a2, b2, c2, index1];
                PxList1(n1,5) = Neighbor26Matrix1(neighborIndex1);
            elseif sizInd1 > 1
                leng1 = zeros(sizInd1,1);
                for n6 = 1:sizInd1
                    leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                end
                [~, Ind3] = min(leng1);
                [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                PxList1(n1+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
            else
                break
            end
        end
        L = PxList1;
    elseif numEndPoint > 2
        for n3 = 1:numNode
            if length(node(1,n3).idx) == 1
                PxList1 = zeros(pxSize, 5);
                PxList1(1,1:4) = [round(node(1,n3).comx), round(node(1,n3).comy), round(node(1,n3).comz), round(node(1,n3).idx)];
                skelBW2 = skelBW;
                
                for n4 = 1:pxSize
                    if ~isempty(find(PxList1(n4,4) == L(:,4), 1))
                        break
                    else
                        % crop neighbor image at 3x3
                        skelBW2(PxList1(n4,1),PxList1(n4,2),PxList1(n4,3)) = 0;
                        neighborIndex1 = find(skelBW2(PxList1(n4,1)-1:PxList1(n4,1)+1, PxList1(n4,2)-1:PxList1(n4,2)+1, PxList1(n4,3)-1:PxList1(n4,3)+1));
                        [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

                        ImgBW1 = false(sizeImg);
                        if ~isempty(neighborIndex1)
                            for n5 = 1:length(neighborIndex1)
                                ImgBW1(PxList1(n4,1)+x2(n5)-2, PxList1(n4,2)+y2(n5)-2, PxList1(n4,3)+z2(n5)-2) = 1;
                            end
                        end
                        index1 = find(ImgBW1);
                        sizInd1 = length(index1);
                        if sizInd1 == 1
                            [a2, b2, c2] = ind2sub(sizeImg, index1);
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1);
                        elseif sizInd1 > 1
                            leng1 = zeros(sizInd1,1);
                            for n6 = 1:sizInd1
                                leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                            end
                            [~, Ind3] = min(leng1);
                            [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                            PxList1(n4+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
                        else
                            break
                        end
                       
                    end

                end
            L = vertcat(L,PxList1);  
            end
                
        end
        empty1 = L(:,1) == 0;
        L(empty1,:) = [];
    else
        display('please check the skeleton image');
        
    end
        
elseif numBranch == 1
    bp1 = find(branchID);
    bp_x = round(node(1,bp1).comx);
    bp_y = round(node(1,bp1).comy);
    bp_z = round(node(1,bp1).comz);
    skelBW1 = skelBW;
    skelBW1(bp_x-1:bp_x+1, bp_y-1:bp_y+1, bp_z-1:bp_z+1) = 0;
    bpInd = sub2ind(sizeImg, bp_x, bp_y, bp_z);
    CC_FLAGE1 = bwconncomp(skelBW1, 26);
    
    if CC_FLAGE1.NumObjects == 3
        segments = cellfun(@length, CC_FLAGE1.PixelIdxList);
        [~, IND] = min(segments(:));
        for n3 = 1:numNode
            if node(1,n3).ep == 1 && isempty(find(CC_FLAGE1.PixelIdxList{1,IND} == node(1,n3).idx, 1))
                PxList1 = zeros(pxSize, 5);
                PxList1(1,1:4) = [round(node(1,n3).comx), round(node(1,n3).comy), round(node(1,n3).comz), round(node(1,n3).idx)];
                skelBW2 = skelBW;
                
                for n4 = 1:pxSize
                    if PxList1(n4,4) == bpInd  || ~isempty(find(PxList1(n4,4) == L(:,4), 1))
                        break
                    else
                        % crop neighbor image at 3x3
                        skelBW2(PxList1(n4,1),PxList1(n4,2),PxList1(n4,3)) = 0;
                        neighborIndex1 = find(skelBW2(PxList1(n4,1)-1:PxList1(n4,1)+1, PxList1(n4,2)-1:PxList1(n4,2)+1, PxList1(n4,3)-1:PxList1(n4,3)+1));
                        [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

                        ImgBW1 = false(sizeImg);
                        if ~isempty(neighborIndex1)
                            for n5 = 1:length(neighborIndex1)
                                ImgBW1(PxList1(n4,1)+x2(n5)-2, PxList1(n4,2)+y2(n5)-2, PxList1(n4,3)+z2(n5)-2) = 1;
                            end
                        end
                        index1 = find(ImgBW1);
                        sizInd1 = length(index1);
                        bp2 = find(bpInd == index1);
                        if sizInd1 == 1
                            [a2, b2, c2] = ind2sub(sizeImg, index1);
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1);
                        elseif sizInd1 > 1 && ~isempty(bp2)
                            [a2, b2, c2] = ind2sub(sizeImg, index1(bp2,1));
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1(bp2,1)];
                        elseif sizInd1 > 1
                            leng1 = zeros(sizInd1,1);
                            for n6 = 1:sizInd1
                                leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                            end
                            [~, Ind3] = min(leng1);
                            [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                            PxList1(n4+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
                        else
                            break
                        end
                    end

                end
            L = vertcat(L,PxList1);  
            end
     

        end
        empty1 = L(:,1) == 0;
        L(empty1,:) = [];
    elseif CC_FLAGE1.NumObjects == 2
        
        for n3 = 1:numNode
            if node(1,n3).ep == 1
                PxList1 = zeros(pxSize, 5);
                PxList1(1,1:4) = [round(node(1,n3).comx), round(node(1,n3).comy), round(node(1,n3).comz), round(node(1,n3).idx)];
                skelBW2 = skelBW;
                
                for n4 = 1:pxSize
                    if PxList1(n4,4) == bpInd  || ~isempty(find(PxList1(n4,4) == L(:,4), 1))
                        break
                    else
                        % crop neighbor image at 3x3
                        skelBW2(PxList1(n4,1),PxList1(n4,2),PxList1(n4,3)) = 0;
                        neighborIndex1 = find(skelBW2(PxList1(n4,1)-1:PxList1(n4,1)+1, PxList1(n4,2)-1:PxList1(n4,2)+1, PxList1(n4,3)-1:PxList1(n4,3)+1));
                        [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

                        ImgBW1 = false(sizeImg);
                        if ~isempty(neighborIndex1)
                            for n5 = 1:length(neighborIndex1)
                                ImgBW1(PxList1(n4,1)+x2(n5)-2, PxList1(n4,2)+y2(n5)-2, PxList1(n4,3)+z2(n5)-2) = 1;
                            end
                        end
                        index1 = find(ImgBW1);
                        sizInd1 = length(index1);
                        bp2 = find(bpInd == index1);
                        if sizInd1 == 1
                            [a2, b2, c2] = ind2sub(sizeImg, index1);
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1);
                        elseif sizInd1 > 1 && ~isempty(bp2)
                            [a2, b2, c2] = ind2sub(sizeImg, index1(bp2,1));
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1(bp2,1)];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1(bp2));
                        elseif sizInd1 > 1
                            leng1 = zeros(sizInd1,1);
                            for n6 = 1:sizInd1
                                leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                            end
                            [~, Ind3] = min(leng1);
                            [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                            PxList1(n4+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
                        else
                            break
                        end
                    end

                end
            L = vertcat(L,PxList1);  
            end
     

        end
        empty1 = L(:,1) == 0;
        L(empty1,:) = [];
    end
elseif numBranch == 2
    skelBW1 = skelBW;
    bp1 = find(branchID);
    bpIndx = zeros(2,5);
    for n7 = 1:2
        bp_x = round(node(1,bp1(n7)).comx);
        bp_y = round(node(1,bp1(n7)).comy);
        bp_z = round(node(1,bp1(n7)).comz);
        skelBW1(bp_x-1:bp_x+1, bp_y-1:bp_y+1, bp_z-1:bp_z+1) = 0;
        bpInd = sub2ind(sizeImg, bp_x, bp_y, bp_z);
        bpIndx(n7,1:4) = [bp_x, bp_y, bp_z, bpInd];
    end
    a3 = abs(bpIndx(1,1)-bpIndx(2,1));
    b3 = abs(bpIndx(1,2)-bpIndx(2,2));
    c3 = abs(bpIndx(1,3)-bpIndx(2,3));
    leng2 = sqrt(a3^2 + b3^2 + c3^2);
    CC_FLAGE1 = bwconncomp(skelBW1, 26);

    if CC_FLAGE1.NumObjects == 3
        segments = cellfun(@length, CC_FLAGE1.PixelIdxList);
        [~, IND] = min(segments(:));
        for n3 = 1:numNode
            if node(1,n3).ep == 1 && ~isempty(find(CC_FLAGE1.PixelIdxList{1,IND} == node(1,n3).idx, 1))
                PxList1 = zeros(pxSize, 5);
                PxList1(1,1:4) = [round(node(1,n3).comx), round(node(1,n3).comy), round(node(1,n3).comz), round(node(1,n3).idx)];
                skelBW2 = skelBW;

                for n4 = 1:pxSize
                    if ~isempty(find(PxList1(n4,4) == bpIndx(:,4), 1))  || ~isempty(find(PxList1(n4,4) == L(:,4), 1))
                        break
                    else
                        % crop neighbor image at 3x3
                        skelBW2(PxList1(n4,1),PxList1(n4,2),PxList1(n4,3)) = 0;
                        neighborIndex1 = find(skelBW2(PxList1(n4,1)-1:PxList1(n4,1)+1, PxList1(n4,2)-1:PxList1(n4,2)+1, PxList1(n4,3)-1:PxList1(n4,3)+1));
                        [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

                        ImgBW1 = false(sizeImg);
                        if ~isempty(neighborIndex1)
                            for n5 = 1:length(neighborIndex1)
                                ImgBW1(PxList1(n4,1)+x2(n5)-2, PxList1(n4,2)+y2(n5)-2, PxList1(n4,3)+z2(n5)-2) = 1;
                            end
                        end
                        index1 = find(ImgBW1);
                        sizInd1 = length(index1);
                        bp2 = find(bpIndx(:,4) == index1);
                        if sizInd1 == 1
                            [a2, b2, c2] = ind2sub(sizeImg, index1);
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1);
                        elseif sizInd1 > 1 && isempty(bp2)
                            [a2, b2, c2] = ind2sub(sizeImg, index1(bp2,1));
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1(bp2,1)];
                                                elseif sizInd1 > 1
                            leng1 = zeros(sizInd1,1);
                            for n6 = 1:sizInd1
                                leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                            end
                            [~, Ind3] = min(leng1);
                            [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                            PxList1(n4+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
                        else
                            break
                        end
                    end

                end
            L = vertcat(L,PxList1);  
            end


        end
        empty1 = L(:,1) == 0;
        L(empty1,:) = [];
        L(end+1,5) = leng2;
    elseif CC_FLAGE1.NumObjects == 2

        for n3 = 1:numNode
            if node(1,n3).ep == 1
                PxList1 = zeros(pxSize, 5);
                PxList1(1,1:4) = [round(node(1,n3).comx), round(node(1,n3).comy), round(node(1,n3).comz), round(node(1,n3).idx)];
                skelBW2 = skelBW;

                for n4 = 1:pxSize
                    if ~isempty(find(PxList1(n4,4) == bpIndx(:,4), 1))  || ~isempty(find(PxList1(n4,4) == L(:,4), 1))
                        break
                    else
                        % crop neighbor image at 3x3
                        skelBW2(PxList1(n4,1),PxList1(n4,2),PxList1(n4,3)) = 0;
                        neighborIndex1 = find(skelBW2(PxList1(n4,1)-1:PxList1(n4,1)+1, PxList1(n4,2)-1:PxList1(n4,2)+1, PxList1(n4,3)-1:PxList1(n4,3)+1));
                        [x2, y2, z2] = ind2sub([3,3,3], neighborIndex1);

                        ImgBW1 = false(sizeImg);
                        if ~isempty(neighborIndex1)
                            for n5 = 1:length(neighborIndex1)
                                ImgBW1(PxList1(n4,1)+x2(n5)-2, PxList1(n4,2)+y2(n5)-2, PxList1(n4,3)+z2(n5)-2) = 1;
                            end
                        end
                        index1 = find(ImgBW1);
                        sizInd1 = length(index1);
                        bp2 = find(bpIndx(:,4) == index1);
                        if sizInd1 == 1
                            [a2, b2, c2] = ind2sub(sizeImg, index1);
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1);
                        elseif sizInd1 > 1 && ~isempty(bp2)
                            [a2, b2, c2] = ind2sub(sizeImg, index1(bp2,1));
                            PxList1(n4+1,1:4) = [a2, b2, c2, index1(bp2,1)];
                            PxList1(n4,5) = Neighbor26Matrix1(neighborIndex1(bp2));
                        elseif sizInd1 > 1
                            leng1 = zeros(sizInd1,1);
                            for n6 = 1:sizInd1
                                leng1(n6,1) = Neighbor26Matrix1(neighborIndex1(n6));
                            end
                            [~, Ind3] = min(leng1);
                            [a2, b2, c2] = ind2sub(sizeImg, index1(Ind3));
                            PxList1(n4+1,1:5) = [a2, b2, c2, index1(Ind3), leng1(Ind3)];
                        else
                            break
                        end
                    end

                end
            L = vertcat(L,PxList1);  
            end


        end
        empty1 = L(:,1) == 0;
        L(empty1,:) = [];
        L(end+1,5) = leng2;
    end
    display('please check the skeleton image just in case');
else
    display('please check the skeleton image');
    
end



