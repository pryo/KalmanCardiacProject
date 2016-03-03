function controlMatrix=getControlMatrix(vref,connectivity,state,deltaT,excitable_threshold,exciting_threshold)
    width = sqrt(size(connectivity,1)); 
    excitedV=44.6639;
    controlMatrix = zeros(size(connectivity,1),size(connectivity,1));
    for i =1:size(controlMatrix,1)
       if state(i)< excitable_threshold
        if connectivity(i,1)~=0&&state(i+1)<exciting_threshold(1)&&state(i+1)>exciting_threshold(2)%&&state(i+1)>exciting_threshold
            %controlMatrix(i,i+1)=1+offset(state(i+1),deltaT,vref);
            controlMatrix(i,i+1) = (excitedV-state(i))/state(i+1);
        elseif connectivity(i,2)~=0&&state(i+width)<exciting_threshold(1)&&state(i+width)>exciting_threshold(2)%%&&state(i+width)>exciting_threshold
            %controlMatrix(i,i+width)=1+offset(state(i+width),deltaT,vref);
            controlMatrix(i,i+width) = (excitedV-state(i))/state(i+width);
        elseif connectivity(i,3)~=0&&state(i-1)<exciting_threshold(1)&&state(i-1)>exciting_threshold(2)%&&state(i-1)>exciting_threshold
            %controlMatrix(i,i-1)=1+offset(state(i-1),deltaT,vref);
            controlMatrix(i,i-1) = (excitedV-state(i))/state(i-1);
        elseif connectivity(i,4)~=0&&state(i-width)<exciting_threshold(1)&&state(i-width)>exciting_threshold(2)%&&state(i-width)>exciting_threshold
            %controlMatrix(i,i-width)=1+offset(state(i-width),deltaT,vref);
            controlMatrix(i,i-width) = (excitedV-state(i))/state(i-width);
        end
       end
    end
end
 
    
    