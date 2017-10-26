function A = sparse(i,j,v,n,m);
% by SeHyoun Ahn, Oct 2017

% This function does not support all use cases of find.
    val = getvalues(v);
    der = getderivs(v);

    [id, jd, vd] = find(der);
    l = size(der, 2);

    retval = sparse(i, j, val, n, m);
    retder = sparse(n*(j(id)-1) + i(id), jd, vd, n*m, l);

    A = myAD(retval, retder);
end
