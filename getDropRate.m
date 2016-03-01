function dropRate = getDropRate(vref,V,deltaT)
index=findV(vref,V);
nextV = Vref(index+deltaT);
dropRate = nextV/V;
end
