function path = gso(map, start, goal)
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


%First Task it to generate population of solutions

nPop = 25;
MaxIt = 40;

%RANGE
range_init = 5.0;
range_boundary = 50.2;

%LUCIFERIN
luciferin_init = 25;
luciferin_decay = 0.4;
luciferin_enhancement = 0.6;

%Neighbors
k_neigh = 20;
beta = 0.5;

%Constants for cost function
k1 = 1;
k2 = 80;
k3 = 0.3;

% Create Empty Glowworm Structure
empty_glowworm.Position=[];
empty_glowworm.range=[];
empty_glowworm.luciferin=[];
empty_glowworm.NV = [];

% Initialize Global Best
GlobalBest.Cost=inf;
GlobalBest.Position = [];

% Create glowworms Matrix
glowworm=repmat(empty_glowworm,nPop,1);
nodes = 50;
i = 1;
while i<=nPop
    [paths,CurrentNode,NVNodes] = getPath(CollisionTest,StartNode,GoalNode,AllPoints(:,:),nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
    %fprintf('%d %d\n',curNodeIndex, nodes);
       if CurrentNode == GoalNode 
           glowworm(i).Position = AllPoints(paths,:);
       elseif CurrentNode ~= GoalNode && NVNodes == (nodes + 1)
               glowworm(i).Position = AllPoints(paths,:);
       elseif NVNodes < nodes
           continue;
       else 
           glowworm(i).Position = AllPoints(paths,:);
       end
       glowworm(i).Position(:,4) = paths(:,1);
       glowworm(i).range = range_init;
       glowworm(i).NV = NVNodes;
       glowworm(i).luciferin = luciferin_init;
       turns  = calc_turns(glowworm(i).Position(:,1:3));
       glowworm(i).Cost = k1*size(paths,1) + k2*pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)) + k3*turns;
       %fprintf('%d %d %d %d %d %d %d \n',AllPoints(CurrentNode,:),AllPoints(GoalNode,:),pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)));
       if glowworm(i).Cost < GlobalBest.Cost
           GlobalBest.Cost = glowworm(i).Cost; 
           GlobalBest.Position = glowworm(i).Position;
       end
       
       i = i + 1;
end
% Algorithm initialization

for it=1:MaxIt
    
    for i=1:nPop
        % Update luciferin
        %disp(glowworm(i).luciferin.x);
        glowworm(i).luciferin = (1-luciferin_decay).*glowworm(i).luciferin + luciferin_enhancement.*(glowworm(i).Cost/10);
        neighbors = [];
        for k =1:nPop
            %dist = pdist2(glowworm(i).Position(end,:),glowworm(k).Position(end,:));
            dist = abs(glowworm(i).Cost - glowworm(k).Cost);
			%if it is in it's range of sigth and it's brightness is higher
            if (dist ~= 0) && (dist <= glowworm(i).range) && (glowworm(i).luciferin <= glowworm(k).luciferin)
                neighbors = [neighbors ; k];
            end
        end
        
        if size(neighbors,1) > 0
         % find the node in the direction of which the glowworm should
         % follow
            li = glowworm(i).luciferin;
            sum_lk = sum(glowworm(i).luciferin);
            neighbors_index = size(neighbors,1);
            %calc probabilties for each neighbor been followed
            probs = zeros(1,neighbors_index);
            for j = 1:neighbors_index
                probs(j) = abs(sum(glowworm(j).luciferin) - sum(li));
            end
            probs = probs./abs( sum_lk - sum(size(probs,2)*li));
           
            %calc prob range
            acc = 0;
            wheel(1,size(probs,2)) = 0;
            for val = 1:size(probs,2)
                acc = acc + probs(val);
                wheel(val) = acc;
            end
            
            %wheel(-1) = 1 ;
            %randomly choice a value for wheel selection method
            rand_val = rand;
            following = i;
            for k = 1:size(wheel,1)
                if rand_val <= wheel(k)
                    following = k;
                end
            end

            toward_index = following;
            %Position update 
            glowworms = glowworm(i).Position;
            toward = glowworm(toward_index).Position;
            
            
            [newPath,Cost] = changePath(glowworms(:,4),toward(:,4),AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
            newPath = removeLoops(newPath);
            glowworm(i).Position = AllPoints(newPath(:,1),:);
            glowworm(i).Position(:,4) = newPath(:,1);
            glowworm(i).Cost = Cost;
            %new_position = glowworms + step_size.*(toward-glowworms)./normV;
           % glowworm(i).Position.x = new_position;
        elseif size(neighbors)== 0
          [newPath,Cost] = changePath(glowworm(i).Position(:,4),[],AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
          newPath = removeLoops(newPath);  
           glowworm(i).Position = AllPoints(newPath(:,1),:);
            glowworm(i).Position(:,4) = newPath(:,1);
            glowworm(i).Cost = Cost;
            
        end
        glowworm(i).range = min(range_boundary,max(0.1,glowworm(i).range + (beta*(k_neigh-size(neighbors,1)))));
        
        if glowworm(i).Cost < GlobalBest.Cost
            GlobalBest.Cost = glowworm(i).Cost; 
            GlobalBest.Position = glowworm(i).Position;
        end
    end
%     if mod(it,5) && it ~= 1
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

