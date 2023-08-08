function P = getProjectionMatrix(DLevelSetFn,params)
    n = size(params,2);
    dc = DLevelSetFn(params);
    P = eye(n) - dc' * inv(dc * dc') * dc;
end