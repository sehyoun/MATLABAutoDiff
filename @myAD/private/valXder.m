function der = valXder(val, der)
    der = val(:,ones(size(der,2),1)) .* der;
