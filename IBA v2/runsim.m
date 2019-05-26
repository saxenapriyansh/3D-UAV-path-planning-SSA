close all;
clear all;
clc;
addpath(genpath('./'));

%% Plan path
disp('Planning ...');
map = load_map('maps/map3.txt', 0.1, 0.5, 0.25);
 %start = { [2 10 2],[1 -4 1],[9.2 17 3],[9.2 10 3],[0.1 10 2]};
% stop  = {[1 -4 1],[0.1 17 3],[9 -4 1],[0.9 -4 5],[9 10 2]};
start = { [0 1 5],[0 2 5],[0 3 5],[19 4 5],[19 5 5]};
stop  = {[19 0 5],[19 5 5],[19 4 5],[0 3 5],[0 1 5]};

visited_nodes = [];
nquad = length(start);
for qn = 1:nquad
    v = cputime;
    [path{qn},visited_nodes] = iba(map, start{qn}, stop{qn},visited_nodes);
    c = cputime - v;
    fprintf('Algo Execution time = %d \n',c);
end
if nquad == 1
    plot_path(map, path{1});
else
    % you could modify your plot_path to handle cell input for multiple robots
    for qn = 1:nquad
         plot_path(map, path{qn});
    end
end

%% Additional init script
init_script;

%% Run trajectory
 v = cputime;
trajectory = test_trajectory(start, stop, map, path, true); % with visualization
 c = cputime - v;
 fprintf('Simulation Execution time = %d \n',c);