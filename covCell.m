mean=sum(Vref)/length(Vref);
diff =abs( Vref(:)-mean);
squre = diff*diff';
cov= squre/length(diff);