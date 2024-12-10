function [prop,flag] = generateProposal(param, ValidParam,step_size)
flag = false;
prop = param;
for i=1:length(param)
    v = normrnd(0,step_size(i));
    prop(i) = prop(i) + v;
end
if(~ValidParam(prop))
    flag = true;
end
end