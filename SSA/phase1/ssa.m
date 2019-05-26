function path = pso(map, start, goal)
% PSO Find the shortest path from start to goal.
%   PATH = PSO(map, start, goal) returns an M-by-3 matrix, where each row
%   consists of the (x, y, z) coordinates of a point on the path.  The first
%   row is start and the last row is goal.  If no path is found, PATH is a
%   0-by-3 matrix.  Consecutive points in PATH should not be farther apart than
%   neighboring cells in the map (e.g.., if 5 consecutive points in PATH are
%   co-linear, don't simplify PATH by removing the 3 intermediate points).
%
%   PATH = PSO(map, start, goal, astar) finds the path using euclidean
%   distance to goal as a heuristic if astar is true.
%
%   [PATH, NUM_EXPANDED] = PSO(...) returns the path as well as
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


nPop = 25;          % SSA - Number of search agents
MaxIt = 40;         % SSA - Maximum numbef of iterations



% w=1;                % Inertia Weight
% wdamp=0.98;         % Inertia Weight Damping Ratio
% c1=1.5;             % Personal Learning Coefficient
% c2=1.5;             % Global Learning Coefficient

% % Constriction Coefficient
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;               % Inertia Weight
% wdamp=1;             % Inertia Weight Damping Ratio
% c1=chi*phi1;         % Personal Learning Coefficient
% c2=chi*phi2;         % Global Learning Coefficient

% alpha=0.1;
% VelMax = 100;    % Maximum Velocity
% VelMin = 0;                  % Minimum Velocity
%Constants for cost function
k1 = 1;
k2 = 80;
k3 = 0.3;

% Create Empty Particle Structure
empty_particle.Position=[];
% empty_particle.Velocity=[];
empty_particle.Cost=[];
% empty_particle.Best.Position=[];
% empty_particle.Best.Cost=[];
empty_particle.NV = [];

% Initialize Global Best
GlobalBest.Cost=inf;            %SSA - Food Fitness
GlobalBest.Position = [];       %SSA - Food Position

div = 20;
% Create Particles Matrix
particle=repmat(empty_particle,nPop,1);     %SSA - Salps
nodes = 200;
i = 1;
while i<=nPop
    [paths,CurrentNode,NVNodes] = getPath(CollisionTest,StartNode,GoalNode,AllPoints(:,:),nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));
    %fprintf('%d %d\n',curNodeIndex, nodes);
       if CurrentNode == GoalNode 
           particle(i).Position = AllPoints(paths,:);
       elseif CurrentNode ~= GoalNode && NVNodes == (nodes + 1)
               particle(i).Position = AllPoints(paths,:);
       elseif NVNodes < nodes
           continue;
       else 
           particle(i).Position = AllPoints(paths,:);
       end
       particle(i).Position(:,4) = paths(:,1);
%        particle(i).Velocity = 0;
       particle(i).NV = NVNodes;
       turns  = calc_turns(particle(i).Position(:,1:3));
       particle(i).Cost = k1*size(paths,1) + k2*pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)) + k3*turns;
       %fprintf('%d %d %d %d %d %d %d \n',AllPoints(CurrentNode,:),AllPoints(GoalNode,:),pdist2(AllPoints(CurrentNode,:),AllPoints(GoalNode,:)));
       
         % Update Personal Best
        particle(i).Best.Position=particle(i).Position;
        particle(i).Best.Cost=particle(i).Cost;
       
       if particle(i).Cost < GlobalBest.Cost
           GlobalBest.Cost = particle(i).Cost; 
           GlobalBest.Position = particle(i).Position;
       end
       
       i = i + 1;
end


[sorted_salps_fitness,sortIdx] = sort([particle.Cost]);
particle = particle(sortIdx);

% Algorithm initialization
for it=1:MaxIt
    
    c1 = 2*exp(-(4*l/Max_iter)^2); %SSA - Eq. (3.2) in the paper
    
    for i=1:nPop
        %SalpPositions= SalpPositions';
        % -------------Neeed to look above line-----------------%
        
        
%     	particle(i).Velocity = abs(w*particle(i).Velocity ...
%             + c1*rand*((particle(i).Best.Cost-particle(i).Cost)/10) ...
%             + c2*rand*((GlobalBest.Cost-particle(i).Cost)/10));
%          % Update Velocity Bounds
%         particle(i).Velocity = max(particle(i).Velocity,VelMin);
%         particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
%         vel = particle(i).Velocity;
          newPath = particle(i).Position(:,4);
         
         % Update Position
         
         if it<=nPop/2
            for j=1:1:3
                c2=rand();
                c3=rand();
                %%%%%%%%%%%%% % Eq. (3.1) in the paper %%%%%%%%%%%%%%
                if c3<0.5 
                    SalpPositions(j,i)=GlobalBest.Position(j)+c1*((ub(j)-lb(j))*c2+lb(j));
                else
                    SalpPositions(j,i)=GlobalBest.Position(j)-c1*((ub(j)-lb(j))*c2+lb(j));
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
         
         
         
         [newPath,Cost] = changePath(newPath,vel,AllPoints(:,:),CollisionTest,StartNode,GoalNode,nodes,Boundaryinitial,Boundaryfinal,size(DummyX,1),size(DummyY,1),size(DummyZ,1));   
         newPath = removeLoops(newPath);
            particle(i).Position = AllPoints(newPath(:,1),:);
            particle(i).Position(:,4) = newPath(:,1);
            particle(i).Cost = Cost;
         % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
        end
        if particle(i).Cost < GlobalBest.Cost
           GlobalBest.Cost = particle(i).Cost; 
           GlobalBest.Position = particle(i).Position;
       end


    end
  % Inertia Weight Damping
    w=w*wdamp;
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

