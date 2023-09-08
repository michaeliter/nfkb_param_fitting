function p = invert_constraint(params,a,b)
%INVERT_CONSTRAINT undo the contraint put by constrain.m
p = log((params-a)./(b-params));
end

