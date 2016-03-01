function dropRate = getDropRate(vref,state,deltaIndex)
dropRate = zeros(1,size(state,1));
for i = 1:length(state)
    
    index=findV(vref,state(i));
    if (index+deltaIndex)<5000
        nextV = vref(index+deltaIndex);
    else
        nextV = vref(5000);
    end
    dropRate(i) = nextV/vref(index);
end
end
