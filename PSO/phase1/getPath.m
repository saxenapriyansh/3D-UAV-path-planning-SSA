function [ paths,CurrentNode,NVNodes ] = getPath(CollisionTest, StartNode,GoalNode,AllPoints,nodes,Boundaryinitial,Boundaryfinal,s1,s2,s3 )
%GETPATH Summary of this function goes here
%   Detailed explanation goes here

    MaxNoofNodes = size(AllPoints,1);
    NodesEvaluated = zeros(MaxNoofNodes,1);
    NodesEvaluated(StartNode) = 1;
    CurrentNode = StartNode;
    paths = StartNode;
%     curNodeIndex = 1;
    NVNodes = 1;
    prevNode = GoalNode;
    while lt(NVNodes,nodes) && ~isequal(CurrentNode,GoalNode)
        
        mov_x = AllPoints(CurrentNode,1)>AllPoints(GoalNode,1);
        mov_y = AllPoints(CurrentNode,2)>AllPoints(GoalNode,2);
        mov_z = AllPoints(CurrentNode,3)>AllPoints(GoalNode,3);
        
        if(isequal(AllPoints(CurrentNode,1),AllPoints(GoalNode,1)))
            mov_x = 2;
        end
        if(isequal(AllPoints(CurrentNode,2),AllPoints(GoalNode,2)))
            mov_y = 2;
        end
        if(isequal(AllPoints(CurrentNode,3),AllPoints(GoalNode,3)))
            mov_z = 2;
        end
        directions = [mov_x;mov_y;mov_z];
        
        [CurrentNeighbours,ViableNeighbours] = Neighbours(CurrentNode,s1,s2,s3, Boundaryinitial, Boundaryfinal,AllPoints,directions);
        
        %Collision test on CurrentNeighbours and ViableNeighbours
        CurrentNeighbours = CurrentNeighbours(~CollisionTest(CurrentNeighbours));
        ViableNeighbours = ViableNeighbours(~CollisionTest(ViableNeighbours));
        
        % Delete already visited nodes from CurrentNeighbours and ViableNeighbours
        CurrentNeighboursABS = CurrentNeighbours(NodesEvaluated(CurrentNeighbours) == 0 ); % not visited current neighbours
        ViableNeighbours = ViableNeighbours(ViableNeighbours ~= prevNode);%delete previous node
        ViableNeighboursABS = ViableNeighbours(NodesEvaluated(ViableNeighbours) == 0);% not visited viable neighbours
        
%         flag = 0;
        if ~isempty(ViableNeighboursABS)
            pos = ceil(rand*size(ViableNeighboursABS,1));
            pos = ViableNeighboursABS(pos);
        elseif isempty(ViableNeighboursABS) && ~isempty(CurrentNeighboursABS)
            pos = ceil(rand*size(CurrentNeighboursABS,1));
            pos = CurrentNeighboursABS(pos);
%             flag = 1;
        elseif isempty(CurrentNeighboursABS) && ~isempty(ViableNeighbours)
            pos = ceil(rand*size(ViableNeighbours,1));
            pos = ViableNeighbours(pos);
        elseif isempty(ViableNeighbours)
            pos = ceil(rand*size(CurrentNeighbours,1));
            pos = CurrentNeighbours(pos);
%             flag = 1;
        else 
            break;
        end
%             if flag == 1 
                NVNodes = NVNodes + 1; 
                %fprintf('%d %d %d\n',AllPoints(pos,:));
%             end
            paths = [paths;pos];
            NodesEvaluated(pos) = 1;
            prevNode = CurrentNode;
            CurrentNode = pos;  
%             curNodeIndex = curNodeIndex + 1;
        %fprintf('%d %d %d\n',curNodeIndex,nodes,~isequal(CurrentNode,GoalNode));
    end


end

