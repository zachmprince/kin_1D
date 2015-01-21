function aux=create_material_prop(time_dep,value1,value2,t1,t2,space_dep,space_function)


if strcmp(time_dep,'constant_in_time')
    aux.ft = inline(0.*t+value1);
elseif strcmp(time_dep,'ramp_in_time')
    % ramp tat transitions for val1 to val2 when t is between t1 and t2
    if(abs(t2-t1)<1e-10)
        error('must use diffrent time steps in ramp definition');
    end
    if(t2<t1)
        error('we must have t2>t1 in ramp definition')
    end
    aux.ft = value1 + (value2-value1)*ramp(t,t1,t2)
else
    error('unknown time dependence');
end

if strcmp(space_dep,'constant_in_space')
    aux.fx = inline(0.*t+1);
elseif strcmp(space_dep,'spatial_dep')
    aux.fx = space_function;
else
    error('unknown spatial dependence');
end

return
end
