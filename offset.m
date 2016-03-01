function offSet=offset(V,deltaT,vref)
index =findV(vref,V)+deltaT;
if index<5000
    
    offSet=vref(index)-V;
else offSet = vref(5000)-V;
end
end