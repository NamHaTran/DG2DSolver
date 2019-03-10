function out=minmod(vector)
numberArgument=length(vector);
absVectorArguments=[];
isSameSign=true;
for i=1:numberArgument-1
    if (getSignOfDouble(vector(i)) ~= getSignOfDouble(vector(i+1)))
        isSameSign=false;
        break;
    end
    absVectorArguments=[absVectorArguments abs(vector(i))];
end
absVectorArguments=[absVectorArguments abs(vector(numberArgument))];
if isSameSign
    out=min(absVectorArguments)*getSignOfDouble(vector(1));
else
    out=0;
end