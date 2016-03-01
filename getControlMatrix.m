function controlMatrix=getControlMatrix(connectivity,state,deltaT,threshHold)
    width = sqrt(size(connectivity,1)); 
    controlMatrix = zeros(size(connectivity,1),size(connectivity,1));
    for i =1:size(controlMatrix,1)
       if state(i)< threshHold
        if connectivity(i,1)~=0
            controlMatrix(i,i+1)=1+offset(state(i+1),deltaT);
        elseif connectivity(i,2)~=0
            controlMatrix(i,i+width)=1+offset(state(i+width),deltaT);
        elseif connectivity(i,3)~=0
            controlMatrix(i,i-1)=1+offset(state(i-1),deltaT);
        elseif connectivity(i,4)~=0
            controlMatrix(i,i-width)=1+offset(state(i-width),deltaT);
        end
       end
    end
end

        
    
    