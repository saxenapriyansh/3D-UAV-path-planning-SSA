function CurrentNeighbours = Neighbours(CurrentNode, sizeX, sizeY, sizeZ,Boundaryinitial, Boundaryfinal, AllPoints)
% sizeX = size(DummyX,1);
% sizeY = size(DummyY,1);
% sizeZ = size(DummyZ,1);
% CurrentNode = 2;
% CurrentNode = 25275;
if sizeZ==0
    A = sizeY*sizeX;
else
    A = sizeZ*sizeY*sizeX;
end
CurrentNeighbours = [ CurrentNode+1; CurrentNode-1; CurrentNode + sizeY; ...
                        CurrentNode - sizeY; CurrentNode + sizeY*sizeX;...
                        CurrentNode - sizeY*sizeX];

if AllPoints(CurrentNode,1)==Boundaryinitial(1)
    CurrentNeighbours(4) = NaN;
end
if AllPoints(CurrentNode,1)==Boundaryfinal(1)
    CurrentNeighbours(3) = NaN;
end
if AllPoints(CurrentNode,2)==Boundaryinitial(2)
    CurrentNeighbours(2) = NaN;
end
if AllPoints(CurrentNode,2)==Boundaryfinal(2)
    CurrentNeighbours(1) = NaN;
end
if AllPoints(CurrentNode,3)==Boundaryinitial(3)
    CurrentNeighbours(6) = NaN;
end
if AllPoints(CurrentNode,3)==Boundaryfinal(3)
    CurrentNeighbours(5) = NaN;
end
CurrentNeighbours(isnan(CurrentNeighbours(:,1)),:)=[];
CurrentNeighbours = (CurrentNeighbours(CurrentNeighbours>=0 & CurrentNeighbours<=A));

end
% for i=1:size(CurrentNeighbours,1)
% %     a = sum(abs(AllNodes(1:j,1:3)-ones(j,1)*CurrentNeighbours(i,1:3)),2)<1e-5;
% %     a = ismember(AllNodes(:,1:3),CurrentNeighbours(i,1:3),'rows');
%     if (CurrentNeighbours(i,1)<map(1,4)||CurrentNeighbours(i,1)>map(2,4)||...
%             CurrentNeighbours(i,2)<map(1,5)||CurrentNeighbours(i,2)>map(2,5)||...
%             CurrentNeighbours(i,3)<map(1,6)||CurrentNeighbours(i,3)>map(2,6))
%         CurrentNeighbours(i,:) = NaN;
%     elseif any(a==1) 
%         CurrentNeighbours(i,:) = NaN;
%     else
%         check = trycollide(map,CurrentNeighbours(i,:));
%         if check == 1
%             CurrentNeighbours(i,:) = NaN;
%         end
%     end
% end
% CurrentNeighbours(isnan(CurrentNeighbours(:,1)),:)=[];
% CurrentNeighbours(:,4) = 0;
% % j=find(AllNodes(:,4), 1, 'last' );
% % CurrentNeighbours(:,4) = j+1:j+size(CurrentNeighbours,1);
% % AllNodes(j+1:j+size(CurrentNeighbours,1),:) = CurrentNeighbours;
% % AllNodes = unique(AllNodes,'rows');
% end
