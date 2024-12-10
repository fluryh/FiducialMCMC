function isARIMA = isARIMAParam(param)
isARIMA = (param(3) > 0)*(abs(param(1))<1);
end