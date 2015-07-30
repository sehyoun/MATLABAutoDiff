function x = mpower(x,y)
% In Package myAD - Automatic Differentiation
% by Martin Fink, June 2006
% martinfink 'at' gmx.at

if (length(x)~=1 & length(y)~=1)
    error('Matrix multiplication not implemented, yet. Use ".*" if possible.');
end

if isa(y, 'myAD')
    if isa(x, 'myAD')
        temp1 = x.values.^(y.values);
        temp2 = temp1.*log(x.values);
        temp3 = y.values.*x.values.^(y.values-1);
        ssOnes = ones(size(x.derivatives,2),1);
        x.derivatives = temp3(:,ssOnes).*x.derivatives + temp2(:,ssOnes).*y.derivatives;
        x.values = temp1;
    else
        y.values = x.^y.values;
        temp = y.values.*log(x);
        y.derivatives = temp(:,ones(size(y.derivatives,1),1)).*y.derivatives;
        x = y;
    end
else
    temp = y.*x.values.^(y-1);
    x.derivatives = temp(:,ones(size(x.derivatives,1),1)).*x.derivatives;
    x.values = x.values.^y;
end
