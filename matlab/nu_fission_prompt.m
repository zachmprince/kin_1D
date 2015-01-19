function y=nu_fission_prompt(curr_time,x)

global dat

y=(1-dat.beta_tot) * nu_fission(curr_time,x);

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
