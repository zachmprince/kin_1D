function y=nu_fission_delayed(curr_time,x)

global dat

y= dat.beta_tot * nu_fission(curr_time,x);

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
