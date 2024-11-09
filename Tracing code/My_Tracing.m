RecVol = importdata('Volume-ICNP1.mat');

BoxSize1=5; 
[boxCoordinates, BoxCenter, BoxRadius, sphere] = create_box(BoxSize1); 
M = find_maxima3D(RecVol,sphere);
M(RecVol<10)=0;
M=M.*(My_paddzero(ones(size(RecVol)-10),size(RecVol)));
PossiblePosNum=sum(M(:));
fprintf('There are possiblely %d atoms in this particle.\n',PossiblePosNum);
Peaks = RecVol.*M; 
PeaksBox3=zeros(size(RecVol));

ind_find=find(M);
[yy,xx,zz] = ind2sub(size(M),ind_find);
tempPeakList=zeros(1,PossiblePosNum);
for kkk=1:PossiblePosNum
    IntensityBox3 = sum(sum(sum(RecVol(yy(kkk)-1:yy(kkk)+1,xx(kkk)-1:xx(kkk)+1,zz(kkk)-1:zz(kkk)+1))));
    PeaksBox3(yy(kkk),xx(kkk),zz(kkk))=IntensityBox3;
    tempPeakList(kkk)=IntensityBox3;
end
[PeakList, ind] = sort(tempPeakList,'descend');
PeakPosList=[yy(ind),xx(ind),zz(ind)]';


NumOfBins=100;
[Count,BinEdges] = histcounts(PeakList,NumOfBins);
BinLocation=(BinEdges(1:end-1)+BinEdges(2:end))/2;
TF = islocalmin(smooth(Count,5)); 
K=find(TF);
figure;plot(BinLocation,Count,BinLocation(TF),Count(TF),'r*');title(['histogram of Local Maximum in Reconstruction, with ',num2str(length(K)),' local minimums exist'])
flg=input('\n choose which local minimum as threshold? \n(input interger 1,2,3,etc, default as 1; if not satisfy with this result, input 0)');
if isempty(flg)
    fprintf('Use the first maximum By Default\n');
    flg=1; 
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
elseif flg==0
    figure;histogram(PeakList,NumOfBins)
    threshold=round(input('Manual choose a threshold?\n'));
% % elseif flg<1 || flg>length(K)
    fprintf('Warning! Use the first minimum in default\n');
    flg=1;
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
else
    flg=round(flg);
    tempLocation=K(flg);
    threshold=BinLocation(tempLocation);
end
Th = threshold/2;
figure;histogram(PeakList,NumOfBins);hold on;plot([threshold,threshold],[1,max(Count)],"LineWidth",1);plot([Th,Th],[1,max(Count)],"LineWidth",1);hold off
figure;plot(PeakList);hold on;plot([1,PossiblePosNum],[threshold,threshold],"LineWidth",1);plot([1,PossiblePosNum],[Th,Th],"LineWidth",1);hold off

NumofRoughPos=sum(PeakList>threshold);
fprintf('There are approximately %d atoms in this particle.\n',NumofRoughPos);
RoughAtomPosList=PeakPosList(:,1:NumofRoughPos);
RoughAtomintensityList=PeakList(1:NumofRoughPos);


Recvol_medfilt = medfilt3(RecVol,[5 5 5]);
V=rescale(Recvol_medfilt);
[Count,~] = histcounts(V);
otsu_th=otsuthresh(Count);
BW = imbinarize(V,otsu_th); 
BW=bwareaopen(BW,round(length(RecVol(:))/100));
BW=imclose(BW,strel('sphere',5));
BW=imcomplement(bwareaopen(imcomplement(BW),round(length(RecVol(:))/100)));

while length(table2array(regionprops3(bwconncomp(BW),"Volume")))>1 % 还是有多个连通区域，那就可能是表面一些独立的区域，需要手动仔细检查
    BW=imopen(BW,strel('sphere',10));
    BW=bwareaopen(BW,round(length(RecVol(:))/100));
end

RemoveList=[];
for i=1:NumofRoughPos
    if BW(RoughAtomPosList(1,i),RoughAtomPosList(2,i),RoughAtomPosList(3,i))==0
        RemoveList=[RemoveList,i];
    end
end
fprintf('%d atom postions are deleted because they are far away from the bulk particle\n',length(RemoveList))
AtomtypeForDeleteShow=ones(1,NumofRoughPos);AtomtypeForDeleteShow(RemoveList)=3;
RoughAtomintensityList1=RoughAtomintensityList;
RoughAtomintensityList1(RemoveList)=[];
RoughAtomPosList1=RoughAtomPosList;
RoughAtomPosList1(:,RemoveList)=[];
NumofRoughPos1=length(RoughAtomintensityList1);

Res = 0.3434; % 5.5 Mx 512 pixels for Themis G2
RoughAtomPosListForTrace=RoughAtomPosList1';
MaxIter = 14;
CritIter = 7;
minDist = 2/Res;
SearchRad = 3;
Q = 0.5;
Alpha = 1;
cropHalfSize = SearchRad;
% get polynomial power array
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];                %#ok<SAGROW>
                else                
                    fitCoeff(end+1,:) = [i j k 0];                 %#ok<SAGROW>
                end
            end
        end
    end
end

[X,Y,Z]      = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd    = find(X.^2+Y.^2+Z.^2 <=(SearchRad+0.5)^2);
XYZdata.X    = X(SphereInd); XYZdata.Y = Y(SphereInd); XYZdata.Z = Z(SphereInd);

Orders       = fitCoeff(:,1:3);
PosArr       = zeros(size(RoughAtomPosListForTrace));
TotPosArr    = zeros(size(RoughAtomPosListForTrace));

exitFlagArr  = zeros(1, size(RoughAtomPosListForTrace,1));
CoeffArr     = repmat(fitCoeff(:,4),[1 size(RoughAtomPosListForTrace,1)]);
tempStr='';
tempReverseStr='';


% perform the main tracing loop
for i=1:size(RoughAtomPosListForTrace,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4; 
          endFlag = 1;
        end
        cropXind = RoughAtomPosListForTrace(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = RoughAtomPosListForTrace(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = RoughAtomPosListForTrace(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = RecVol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX == -100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1; 
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2; 
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+RoughAtomPosListForTrace(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                        exitFlagArr(i) = -3; % too close
                    end
                        TotPosArr(i,:) = PosArr(i,:) + RoughAtomPosListForTrace(i,:);
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    
    if exitFlagArr(i)~=0
        fprintf('peak %d, flag %d \n',i,exitFlagArr(i));
    end
    
    if mod(i,100)==0
        fprintf(tempReverseStr);
        tempStr=sprintf('%d atoms have been traced',i);
        fprintf(tempStr);
        tempReverseStr=repmat('\b',1,length(tempStr)); 
    end
%     if mod(i,100)==0
%         fprintf('\n--%d atoms have been traced--\n',i);
%     end
end

TracedPos=TotPosArr((exitFlagArr==0),:)';
AtomTypeListForAutoTracingShow=(exitFlagArr==0);
AtomTypeListForAutoTracingShow=AtomTypeListForAutoTracingShow(AtomTypeListForAutoTracingShow>0);

save('.\Trced_positions.mat','RecVol','TracedPos','AtomTypeListForAutoTracingShow');