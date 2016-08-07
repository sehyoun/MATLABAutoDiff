function x=permute(x,order)
% by SeHyoun Ahn, July 2016

if ~issparse(x.values)
    aux=size(x.values);
    p=permute(reshape(1:prod(aux),aux),order);
    x.derivatives=x.derivatives(p(:),:);
    x.values=permute(x.values,order);
else
    if (size(order) == [1,2])
        if (order == [2,1])
            x = x';
        elseif (order == [1,2])
        else
            error('ORDER contains an invalid permutation; try using transpose instead.');
        end
    elseif (size(order) == [2,1])
        if (order == [2;1])
            x=x';
        elseif (order == [1;2])
        else
            error('ORDER contains an invalid permutation; try using transpose instead.');
        end
    else
        error('ORDER contains an invalid permutation; try using transpose instead.');
    end
end