function [lossval] = lossfunc(x,trainlearning,phasedata)
value=floor(trainlearning/(x+0.0001));
lossval=sum((value-phasedata).^2);
end

