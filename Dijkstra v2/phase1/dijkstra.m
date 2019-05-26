function [path] = dijkstra(map, start, goal)
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
startcopy = start;
goalcopy = goal;
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
NodesEvaluated = zeros(MaxNoofNodes,1);
NodesTobeEvaluated = zeros(MaxNoofNodes,1);
TotalScore(1:MaxNoofNodes,1) = (Inf);
GoalScore(1:MaxNoofNodes,1) = (Inf);
GoalScore(StartNode) = 0;
NodesTobeEvaluated(StartNode) = 1;
PreviousNodes = zeros(MaxNoofNodes,1);
NodeCounter = 1;
TotalScore(StartNode,1) = GoalScore(StartNode,1) + Heuristic(StartNode,GoalNode,AllPoints,true);

while ~isempty(NodesTobeEvaluated)
%     A = (GenerateNodes(NodesTobeEvaluated>0));
%     B = TotalScore(GenerateNodes(NodesTobeEvaluated>0));
    [~,CurrentNode] = min(TotalScore);
%     [~,CurrentNode] = min(B);
%     CurrentNode = A(CurrentNode);
%     AllPoints(CurrentNode,:)
%     TotalScore(CurrentNode,:)
    if NodesEvaluated(CurrentNode) == 1
        continue;
    end
    if CurrentNode==GoalNode      
        fprintf('Goal Reached\n');
        path = ReconstructPath(PreviousNodes,CurrentNode);
        break;        
    end   
    NodesEvaluated(CurrentNode) = 1;
    NodesTobeEvaluated(CurrentNode) = 0;
    NodeCounter = NodeCounter + 1;
    CurrentNeighbours = Neighbours(CurrentNode,size(DummyX,1),size(DummyY,1),size(DummyZ,1), Boundaryinitial, Boundaryfinal,AllPoints); 
    CurrentNeighbours = CurrentNeighbours(~CollisionTest(CurrentNeighbours));    
    if isempty(CurrentNeighbours)
        continue;
    end    
    for i=1:size(CurrentNeighbours)
        CurrentNeighbourNodeIndex = CurrentNeighbours(i);         
        if NodesEvaluated(CurrentNeighbourNodeIndex)==1
%             TotalScore(CurrentNeighbourNodeIndex,:) = NaN;
            continue;
        end         
%         TempGoalScore = pdist2(AllPoints(CurrentNode,:),AllPoints(StartNode,:)) + Heuristic(CurrentNode,CurrentNeighbourNodeIndex,AllPoints,astar);        
        TempGoalScore =GoalScore(CurrentNode) + sum(abs(AllPoints(CurrentNode,:) - AllPoints(CurrentNeighbourNodeIndex,:)));
        b = NodesTobeEvaluated(CurrentNeighbourNodeIndex);
        if ~b||TempGoalScore<GoalScore(CurrentNeighbourNodeIndex)
           PreviousNodes(CurrentNeighbourNodeIndex) = CurrentNode;
           GoalScore(CurrentNeighbourNodeIndex) = TempGoalScore;
           TotalScore(CurrentNeighbourNodeIndex)= GoalScore(CurrentNeighbourNodeIndex) + Heuristic(CurrentNeighbourNodeIndex,GoalNode,AllPoints,true);
%            GoalScore(CurrentNeighbourNodeIndex) + Heuristic(CurrentNeighbourNodeIndex,GoalNode,AllPoints,astar);
%            if TotalScore(CurrentNeighbourNodeIndex,:) == 0
%             TotalScore(CurrentNeighbourNodeIndex,:) = NaN;
%            end
           if ~b
               NodesTobeEvaluated(CurrentNeighbourNodeIndex) = 1;               
           end
        else
%             PreviousNodes(CurrentNeighbourNodeIndex,:) = CurrentNodeIndex;
%             TotalScore(CurrentNeighbourNodeIndex) = NaN;
        end
        
    end    
%     title([TotalScore(CurrentNode)]);
%     disp([TotalScore(CurrentNode)]);
%      plot3(AllPoints(CurrentNode,1), AllPoints(CurrentNode,2), AllPoints(CurrentNode,3),'g.');
%      pause(0.000000001);
    TotalScore(CurrentNode) = NaN;
%     [AllPoints(CurrentNode,1), AllPoints(CurrentNode,2), AllPoints(CurrentNode,3)];
   
%     hold on;
%     grid on;
    
%     x = x+1;
end
path = AllPoints(path,:);
path(1,:) = startcopy;
path(end,:) = goalcopy;
plot_path(map,path);
% plot3(path(:,1), path(:,2), path(:,3),'b');
% hold on;
% grid on;    
end

