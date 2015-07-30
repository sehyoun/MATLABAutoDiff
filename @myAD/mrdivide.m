function x = mrdivide(x,y)
% In Package myAD - Automatic Differentiation
% by Martin Fink, May 2007
% martinfink 'at' gmx.at

%%
if isa(y,'myAD')
    error('I will code this later');
else
    if numel(y)==1
        x.derivatives=x.derivatives/y;
        x.values=x.values/y;
    else
        error('Sorry, I have not implemented this yet');
    end
end