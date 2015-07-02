function value=ramp(time,t1,t2)

if(time<=t1)
    value=0; return
elseif(time>=t2)
    value=1; return
else
    value=(time-t1)/(t2-t1);
end

end