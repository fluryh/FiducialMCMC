function [prop,flag] = generateProposal(param, ValidParam)
flag = false;
prop = param;
for i=1:length(param)
    v = normrnd(0,.1);
    prop(i) = prop(i) + v;
end
if(~ValidParam(prop))
    flag = true;
end
end