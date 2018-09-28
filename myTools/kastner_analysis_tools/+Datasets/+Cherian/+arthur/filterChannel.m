function [filtChannel] = filterChannel(channel,filter)
    filtChannel = channel .* filter;
end
