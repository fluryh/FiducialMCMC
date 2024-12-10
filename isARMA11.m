function isARMA = isARMA11(param)
isARMA = (abs(param(1)) < 1) * (param(3) > 0);
end