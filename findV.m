function index= findV(vref, voltage)
%this  function take Vref sequence and voltage then return the index of
%voltage in Vref

head = 1;
tail = length(vref);
seq = vref;
while length(seq)>2
    midIndex = round((head+tail)/2);
    if vref(midIndex)>voltage
        head = midIndex;
    elseif vref(midIndex)<voltage
        tail = midIndex;
    else
        head = midIndex;
        tail = midIndex;
    end
seq =vref(head:tail);
end
if abs(seq(1)-voltage)-abs(seq(2)-voltage)>0
    index =tail;
else 
    index = head;
end
end



        

