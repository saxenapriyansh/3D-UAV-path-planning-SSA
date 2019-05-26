function [path,visited_nodes] = iba(map, start, goal,visited_nodes)
% GSO Find the shortest path from start to goal.
%   PATH = GSO(map, start, goal) returns an M-by-3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path.  The first
%   row is start and the last row is goal.  If no path is found, PATH is a
%   0-by-3 matrix.  Consecutive points in PATH should not be farther apart than
%   neighboring cells in the map (e.g.., if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = GSO(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = GSO(...) returns the path as well as
%   the number of points that were visited while performing the search.

persistent qd;

if qd 
    qd = qd+1;
else 
    qd =1;
end
xy_res = map(3,4);
z_res = map(4,4);

startcollisioncheck = collide(map,start);
goalcollisioncheck = collide(map,goal);
if startcollisioncheck == 1 || goalcollisioncheck == 1
    fprintf('Start or Goal within obstacle\n');
    path = [];
    return;
end

Boundaryinitial = map(1,4:6);
Boundaryfinal = map(2,4:6);
if (Boundaryfinal(3)-Boundaryinitial(3))<z_res
    z_res = Boundaryfinal(3)-Boundaryinitial(3);
end
AllPoints = mygrid(Boundaryinitial,Boundaryfinal,xy_res,z_res);
DummyX = (Boundaryinitial(:,1):xy_res:Boundaryfinal(:,1))';
DummyY = (Boundaryinitial(:,2):xy_res:Boundaryfinal(:,2))';
DummyZ = (Boundaryinitial(:,3):z_res:Boundaryfinal(:,3))';
CollisionTest = collide(map,AllPoints);

if mod(goal(:,1),xy_res)
    goal(:,1) = goal(:,1) - mod(goal(:,1),xy_res);
end
if mod(goal(:,2),xy_res)
%     goal(:,2) = goal(:,2) + (xy_res-mod((goal(:,2)-start(:,2)),xy_res));
    goal(:,2) = goal(:,2) - mod(goal(:,2),xy_res);
end
if mod(goal(:,3),z_res)
    goal(:,3) = goal(:,3) - mod(goal(:,3),z_res);
%     goal(:,3) = goal(:,3) - mod((goal(:,3)-start(:,3)),z_res);
end
if mod(start(:,1),xy_res)
    start(:,1) = start(:,1) - mod(start(:,1),xy_res);
end
if mod(start(:,2),xy_res)
    start(:,2) = start(:,2) - mod(start(:,2),xy_res);
end
if mod(start(:,3),z_res)
    start(:,3) = start(:,3) - mod(start(:,3),z_res);
%     start(:,3) = start(:,3) - mod((start(:,3)-start(:,3)),z_res);
end
% DistanceFromStart = sqrt(sum(bsxfun(@minus, AllPoints, start).^2, 2));
MaxNoofNodes = size(AllPoints,1);
GenerateNodes = (1:MaxNoofNodes)';
% StartNode = GenerateNodes(ismember(AllPoints,start,'rows'));
StartNode = GenerateNodes(sum(abs(AllPoints(:,:)-ones(size(AllPoints,1),1)*start),2)<eps);
% GoalNode = GenerateNodes(ismember(AllPoints,goal,'rows'));
GoalNode = GenerateNodes(sum(abs(AllPoints(:,:)-ones(size(AllPoints,1),1)*goal),2)<eps);


%First Task it to generate population of solutions

nPop = 20;
MaxIt = 25;

Q = 0.5;
r = 0.5;
A = 0.5;

% Create Empty bat Structure
empty_bat.Position=[];
empty_bat.velocity=[];
empty_bat.NV = [];

% Initialize Global Best
GlobalBest.Cost=inf;
GlobalBest.Position = [];

% Create bat Matrix
bat=repmat(empty_bat,nPop,1);
nodes = 200;
i = 1;

while i<=nPop
    [paths,CurrentNode,NVNodes] = getPath(CollisionTest,StartNode,GoalNode,AllPoints(:,:),nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1),visited_nodes);
    %fprintf('%d %d\n',curNodeIndex, nodes);
       if CurrentNode == GoalNode 
           bat(i).Position = AllPoints(paths,:);
       elseif CurrentNode ~= GoalNode && NVNodes == (nodes + 1)
               bat(i).Position = AllPoints(paths,:);
       elseif NVNodes < nodes
           continue;
       else 
           bat(i).Position = AllPoints(paths,:);
       end
       bat(i).Position(:,4) = paths(:,1);
       bat(i).velocity = 0;
       bat(i).NV = NVNodes;
       bat(i).Cost = size(paths,1) + 100*pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:));
       %fprintf('%d %d %d %d %d %d %d \n',AllPoints(CurrentNode,:),AllPoints(GoalNode,:),pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)));
       if bat(i).Cost < GlobalBest.Cost
           GlobalBest.Cost = bat(i).Cost; 
           GlobalBest.Position = bat(i).Position;
       end
       i = i + 1;
end
% Algorithm initialization

for it=1:MaxIt
    bat = nestedSortStruct(bat,'Cost');
    for i=1:nPop
        numbers = randi([1 nPop],1,4);
%         for j = 1:size(numbers,1)
%             if numbers(j) == i
%                 numbers(j) = 
%             end
%         end
        r1  = numbers(1);
        r2 = numbers(2);
        r3 = numbers(3);
        r4 = numbers(4);
        
        bat(i).velocity = bat(i).velocity + (bat(i).Cost - GlobalBest.Cost)*Q;
        
        [newPath,Cost] = changePath(bat(i).Position(:,4),GlobalBest.Position(:,4),AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1),visited_nodes);
         newPath = removeLoops(newPath);
        bat(i).Position = AllPoints(newPath(:,1),:);
        bat(i).Position(:,4) = newPath(:,1);
        bat(i).Cost = Cost;
        batU = empty_bat; %Offspring Bat solution 
        if( rand > r )
            [newPath,Cost] = changePath(GlobalBest.Position(:,4),[],AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1),visited_nodes);                 
        else
            [newPath,~] = changePath(bat(r2).Position(:,4),bat(r3).Position(:,4),AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1),visited_nodes);    
            [newPath,Cost] = changePath(bat(r1).Position(:,4),newPath,AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1),visited_nodes);    
        end
        newPath = removeLoops(newPath);
        batU.Position = AllPoints(newPath(:,1),:);
        batU.Position(:,4) = newPath(:,1);
        batU.Cost = Cost;
        
        
        if bat(r4).Cost <= bat(i).Cost && bat(r4).Cost <= batU.Cost
            batK  = bat(r4);
        elseif bat(i).Cost <= bat(r4).Cost && bat(i).Cost <= batU.Cost
            batK = bat(i);
        elseif batU.Cost <= bat(r4).Cost && batU.Cost <= bat(i).Cost
            batK = batU;
        end
        
        if rand < A
            bat(r4) = batK;
        end
         if bat(i).Cost < GlobalBest.Cost
            GlobalBest.Cost = bat(i).Cost; 
            GlobalBest.Position = bat(i).Position;
         end
         if bat(r4).Cost < GlobalBest.Cost
            GlobalBest.Cost = bat(r4).Cost; 
            GlobalBest.Position = bat(r4).Position;
        end
    end 
end
path = GlobalBest.Position(:,1:3);

if path(end,:) ~= AllPoints(GoalNode,:)
    fprintf('Path for vehicle %d not reached goal\n',qd);
end
visited_nodes = [visited_nodes;GlobalBest.Position(:,4)];
visited_nodes = unique(visited_nodes);
%fprintf('%d %d %d',size(path));
plot_path(map,path);
% plot3(path(:,1), path(:,2), path(:,3),'b');
% hold on;
% grid on;    
end

