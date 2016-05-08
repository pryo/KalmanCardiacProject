function dropRate = getDropRate(vref,state,deltaIndex)
dropRate = zeros(1,size(state,1));


for i = 1:length(state)
    
    index=findV(vref,state(i));
    if (index+deltaIndex)<5000
        nextV = vref(index+deltaIndex);
        dropRate(i) = nextV/vref(index);
    else
        dropRate(i) = 1;
        
    end
    
end
end
