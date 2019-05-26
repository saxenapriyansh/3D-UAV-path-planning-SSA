% Add additional initialization if you want to.
% You can use this space to pre-compute the trajectory to avoid
% repeatedly computing the same trajectory in every call of the
% "trajectory_generator" function

% Generate trajectory
disp('Generating Trajectory ...');
waypoints = cell2mat(path(1));
% path = waypoints(1:25,:);

% clear all
% load waypoints.mat;
xyres = map(3,4);
zres = map(4,4);
modifiedpath = waypoints;
newpath = [];
err = 10e-5;
hisxinc = 0;
hisyinc = 0;
hiszinc = 0;

for i = 1:size(modifiedpath,1)-1  
    temp = waypoints(i,:);
    xinc = 0;
    yinc = 0;
    zinc = 0;
    if abs((waypoints(i+1,1)- temp(:,1))) < (xyres+err) && abs((waypoints(i+1,1)- temp(:,1))) > (xyres-err)
        xinc = 1;
        
    elseif abs((waypoints(i+1,2)- temp(:,2))) < (xyres+err) && abs((waypoints(i+1,2)- temp(:,2))) > (xyres-err)
        yinc = 1;
        
    elseif abs((waypoints(i+1,3)- temp(:,3))) < (zres+err) && abs((waypoints(i+1,3)- temp(:,3))) > (zres-err)
        zinc = 1;
        
    end
    if hisxinc && xinc    
        excluded = temp;
    elseif hisyinc && yinc
        excluded = temp;
    elseif hiszinc && zinc
        excluded = temp;
    else
        newpath(end+1,:) = temp;
    end
    hisxinc = xinc;
    hisyinc = yinc;
    hiszinc = zinc;
           
end
newpath(end+1,:) = waypoints(end,:);
iter1newpath = newpath;
% for i = 1:size(iter1newpath,1)-1
%     temp1 = iter1newpath(i,:);
%     temp2 = iter1newpath(i+1,:);
%     if pdist2(iter1newpath(i,:),iter1newpath(i+1,:))<(0.1+err)
%         iter1newpath(i+1,:) = NaN;
%     end
% end

iter2newpath = [];
iter2newpath(1,:) = iter1newpath(1,:);
i = 1;

while sum(iter2newpath(end,:) ~= iter1newpath(end,:))>0
    j = 0;
    C = 1;
    while C==1
        newpoint = iter1newpath(end-j,:);
        currentpoint = iter2newpath(i,:);
        checkpoints = mygrid(currentpoint,newpoint,xyres,zres);        
        C = collide(map,checkpoints);
        C = max(C);
        j = j + 1;
    end
    iter2newpath(end+1,:) = newpoint;
    i = i + 1;
end

for i = 1:size(iter2newpath,1)-1
    distbtwpts(i,1) =  pdist2(iter2newpath(i,:),iter2newpath(i+1,:));
end
totaldist = sum(distbtwpts); 
tmin = 0;
velmax = 1.4;
% totaltime = totaldist/velmax;
equations = [];
for i = 1:size(iter2newpath,1)-1
    ratio = distbtwpts(i)/totaldist;
    tmax = tmin + 1.3*sqrt(distbtwpts(i));
    currentstart = iter2newpath(i,:);
    currentstop = iter2newpath(i+1,:);
    currentequation = CalculateEquations(currentstart,currentstop,tmin,tmax);
    equations = [equations;currentequation];
    tmin = tmax;
end

save('allequations.mat','equations');
path = iter2newpath;
trajectory_generator([], [], map, path);