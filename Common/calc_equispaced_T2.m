
function T2_axis = calc_equispaced_T2(T2l,T2h,T2base,nT2s,min_T2_distance,round_f,plot_f)
	
t  = T2base;									% [sec]
Al = exp(-t/T2l);
Ah = exp(-t/T2h);
A  = linspace(Al,Ah,nT2s);

T2_axis = round(1e3*(-t ./ log(A)))/1e3;		% [sec]

% Round the T2 values and spread them at a minimal resolution of 'min_T2_distance' [ms]
if (round_f)
	T2_axis = round(T2_axis);
end;

if (min_T2_distance >= 1)
	for idx = 2:nT2s
		while ((T2_axis(idx) - T2_axis(idx-1)) < min_T2_distance)
			T2_axis(idx) = T2_axis(idx) + 1;
		end;
	end;
end;

if (plot_f)
figure;
subplot(211);   plot(T2_axis*1e3     ,'.-');   xlabel('T2 index');   ylabel('T2 [ms]');
subplot(212);   plot(exp(-t./T2_axis),'k-');   xlabel('T2 index');   ylabel('Signal' ); hold on;
end;

return;



% 	if ((T2l>18) || (T2h < 30)), error(0); end;
% 
% 	t(1) = 20e-3;
% 	t(2) = 21e-3;
% 	d    = exp(-1e-3/t(2)) - exp(-1e-3/t(1));
% 
% 	for idx = (t(2)*1000+1):1:T2h
% 		prv_exp_val = exp(-1e-3/t(end));
% 		cur_exp_val = exp(-1e-3/(idx*1e-3));
% 		cur_d       = cur_exp_val - prv_exp_val;
% 		if cur_d > d
% 			t(end+1) = idx*1e-3;
% 		end;
% 	end;
% 
% 	t1=[(T2l*1e-3:1e-3:(t(1)-1e-3)), t];
% 	figure; subplot(121); plot(t1                  ,'.-'); title('t');
% 			subplot(122); plot(exp(-1e-3./t1)      ,'.-'); title('exp(-1./t)');
% 	%       subplot(133); plot(diff(exp(-1e-3./t1)),'.-'); title('diff(exp(-1./t))');
