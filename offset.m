function offSet=offset(V,deltaT,vref)
offSet=vref(findV(vref,V)+deltaT)-V;