function [channelArea, channelAreaDerivative] = calcChannelArea(thisX)
%calcChannelArea Calculate channel area for given x value
%   See exercise for more info on how to calculate channel area

if (0 <= thisX) && (thisX <= 5)
    channelArea = 1.0 + 1.5 * (1.0 - thisX / 5.0)^2;
    channelAreaDerivative = -3.0 * (1.0 - thisX / 5.0) / 5.0;
elseif (5 < thisX) && (thisX <= 10)
    channelArea = 1.0 + 0.5 * (1.0 - thisX / 5.0)^2;
    channelAreaDerivative = -1.0 * (1.0 - thisX / 5.0) / 5.0;
end

end

