function path = bbo(map, start, goal)
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

%BBO algorithm
%First Task it to generate population of solutions

nPop = 25;
MaxIt = 40;

%Constants for cost function
k1 = 1;
k2 = 80;
k3 = 0.3;

KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*nPop);     % Number of Kept Habitats

I = 1; % max immigration rate for each island
E = 1; % max emigration rate, for each island

% Migration Rates
mu=linspace(1,0,nPop);          % Emmigration Rates
Prob = linspace(1,0,nPop); 
ProbDot = linspace(1,0,nPop); 
pMutation(1:nPop) = 0.75;

for j = 1 : nPop
    Prob(j) = 1 / nPop; 
end

% Empty Habitat
empty_habitat.Position=[];
empty_habitat.Cost=[];
empty_habitat.NV = [];

% Initialize Global Best
GlobalBest.Cost=inf;
GlobalBest.Position = [];

% Create glowworms Matrix
habitat=repmat(empty_habitat,nPop,1);
nodes = 200;
i = 1;
while i<=nPop
    [paths,CurrentNode,NVNodes] = getPath(CollisionTest,StartNode,GoalNode,AllPoints(:,:),nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
    %fprintf('%d %d\n',curNodeIndex, nodes);
       if CurrentNode == GoalNode 
               habitat(i).Position = AllPoints(paths,:);
       elseif CurrentNode ~= GoalNode && NVNodes == (nodes + 1)
               habitat(i).Position = AllPoints(paths,:);
       elseif NVNodes < nodes
           continue;
       else 
           habitat(i).Position = AllPoints(paths,:);
       end
       habitat(i).Position(:,4) = paths(:,1);
       habitat(i).NV = NVNodes;
       turns  = calc_turns(habitat(i).Position(:,1:3));
       habitat(i).Cost = k1*size(paths,1) + k2*pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)) + k3*turns;
         
       %fprintf('%d %d %d %d %d %d %d \n',AllPoints(CurrentNode,:),AllPoints(GoalNode,:),pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)));
       if habitat(i).Cost < GlobalBest.Cost
           GlobalBest.Cost = habitat(i).Cost; 
           GlobalBest.Position = habitat(i).Position;
       end
       
       i = i + 1;
end
% Algorithm initialization

for it=1:MaxIt
    maxCost = 0;
    for i = 1:nPop
        if habitat(i).Cost > maxCost
            maxCost = habitat(i).Cost;
        end
    end
    
    for i = 1:nPop
        mu(i) = habitat(i).Cost/maxCost;
    end
    lambda = 1-mu;
    
    
    for i=1:nPop
        newPath = habitat(i).Position(:,4);
         Cost = habitat(i).Cost;
         if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;
                EP=EP/sum(EP);
                
                % Select Source Habitat
                j=RouletteWheelSelection(EP);
                mhabitat = habitat(i).Position;
                toward = habitat(j).Position;
                
                
                % Migration
                [newPath,Cost] = changePath(mhabitat(:,4),toward(:,4),AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
                 newPath = removeLoops(newPath);
        end
         
         habitat(i).Position = AllPoints(newPath,:);
            habitat(i).Position(:,4) = newPath;
            habitat(i).Cost = Cost;
            
        if habitat(i).Cost < GlobalBest.Cost
            GlobalBest.Cost = habitat(i).Cost; 
            GlobalBest.Position = habitat(i).Position;
        end                 
    end 
    %calculation of mutation probability
        maxCost = 0;
    for i = 1:nPop
        if habitat(i).Cost > maxCost
            maxCost = habitat(i).Cost;
        end
    end
    
    for i = 1:nPop
        mu(i) = habitat(i).Cost/maxCost;
    end
    lambda = 1-mu;
    if j < nPop
         ProbPlus = Prob(j+1);
          muPlus = mu(j+1);
    else
        ProbPlus = 0;
          muPlus = 0;
    end
    if i > 1
        ProbMinus = Prob(i-1);
        lambdaMinus = lambda(i-1);
    else
       ProbMinus = 0;
        lambdaMinus = 0;
    end
    ProbDot(i) = -(lambda(i) + mu(i)) * Prob(i) + lambdaMinus * ProbMinus + muPlus * ProbPlus;
    Prob = max(Prob, 0);
    Prob = Prob / sum(Prob); 
    
    P = size(habitat(i).Position,1);
    pMutation = pMutation .* (1 - Prob/P);
    for i = 1:nPop
        newPath = habitat(i).Position(:,4);
         Cost = habitat(i).Cost;
         % Mutation
        if rand<=pMutation(j)
             [newPath,Cost] = changePath(newPath,[],AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
             newPath = removeLoops(newPath);
        end
         habitat(i).Position = AllPoints(newPath,:);
              habitat(i).Position(:,4) = newPath;
              habitat(i).Cost = Cost;
        
        if habitat(i).Cost < GlobalBest.Cost
            GlobalBest.Cost = habitat(i).Cost; 
            GlobalBest.Position = habitat(i).Position;
        end      
    end
%         if mod(it,5) && it ~= 1
%         continue;
%     end
%     fprintf('Iteration Number %d : %d \n',it, GlobalBest.Cost);
end
path = GlobalBest.Position(:,1:3);
fprintf('%d \n', GlobalBest.Cost);
%fprintf('%d %d %d',size(path));
plot_path(map,path);
% plot3(path(:,1), path(:,2), path(:,3),'b');
% hold on;
% grid on;    
end

