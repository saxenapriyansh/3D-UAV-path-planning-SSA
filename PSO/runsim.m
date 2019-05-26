close all;
clear all;
clc;
addpath(genpath('./'));

%% Plan path
disp('Planning ...');
map = load_map('maps/map2.txt', 0.1, 0.5, 0.25);
%  start = { [1 -4 1]};
%  stop  = {[0.1 17 3]};

start = { [0.5 7 1]};
stop  = {[8.5 3 3.5]}; % MAP 2

% start = {[0 1 5]};
% stop = {[19 1 5]};
nquad = length(start);
% for q=1:10
for qn = 1:nquad
    v = cputime;
    path{qn} = pso(map, start{qn}, stop{qn});
    c = cputime - v;
    fprintf('Algo Execution time = %d \n',c);
end
% end
if nquad == 1
    plot_path(map, path{1});
else
    % you could modify your plot_path to handle cell input for multiple robots
end

%% Additional init script
init_script;

%% Run trajectory
trajectory = test_trajectory(start, stop, map, path, true); % with visualization
