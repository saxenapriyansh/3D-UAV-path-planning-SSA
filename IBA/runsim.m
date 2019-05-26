close all;
clear all;
clc;
addpath(genpath('./'));

%% Plan path
disp('Planning ...');
map = load_map('maps/map4.txt', 0.1, 0.5, 0.25);
start = { [1 -4 1]};
stop  = {[0.1 17 3]};

%start = {[0 1 5]};
%stop = {[19 1 5]};
nquad = length(start);
for qn = 1:nquad
    v = cputime;
    path{qn} = iba(map, start{qn}, stop{qn});
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
trajectory = test_trajectory(start, stop, map, path, true); % with visualization