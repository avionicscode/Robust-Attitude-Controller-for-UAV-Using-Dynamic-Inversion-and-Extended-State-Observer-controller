function out = fal( e,alfa,delta )

if abs(e)> delta 
out=abs(e)^alfa*sign(e);
elseif abs(e)<=delta
    out=e/(delta^(1-alfa));
end

