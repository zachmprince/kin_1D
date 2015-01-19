function y=sigma_a(curr_time,x)

global dat

y0=1.1; a=0.005;

% set value
y=y0;

% short cut for rod movement times
t_beg_1 = dat.rod_mov_times.t_beg_1;
t_end_1 = dat.rod_mov_times.t_end_1;
t_beg_2 = dat.rod_mov_times.t_beg_2;
t_end_2 = dat.rod_mov_times.t_end_2;

% hard-coded locations
pos_beg_1=0.2*dat.width;
pos_end_1=0.3*dat.width;
pos_beg_2=0.25*dat.width;
pos_end_2=0.3*dat.width;
pos_beg_3=0.7*dat.width;
pos_end_3=0.8*dat.width;

if(curr_time >= t_beg_1)
    if(x>=pos_beg_1 && x<=pos_end_1)
        y=y - a*min(curr_time-t_beg_1,t_end_1-t_beg_1)/(t_end_1-t_beg_1);
    end
end

if(curr_time >= t_beg_2)
    if((x>=pos_beg_2 && x<=pos_end_2) || (x>=pos_beg_3 && x<=pos_end_3) )
        y=y + a*min(curr_time-t_beg_2,t_end_2-t_beg_2)/(t_end_2-t_beg_2);
    end
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
