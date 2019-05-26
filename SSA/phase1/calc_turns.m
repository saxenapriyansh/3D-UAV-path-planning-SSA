function turns = calc_turns(path)
turns = 0;
lim = size(path,1);
for i = 2:lim-1
    if isequal(path(i,:)-path(i-1,:),path(i+1,:)-path(i,:))
        turns = turns + 1; 
    end
end

end