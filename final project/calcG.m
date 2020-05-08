function gab = calcG(twists,angles,init)
    gab = init;
    
    for i = size(twists,2):-1:1
        gab = twistExp(twists(:,i),angles(i))*gab;
    end
end