function [t1,y1,t3,y3,t4,y4]=Run_hybrid(params,y0,Beta)
 tspan = [0 200];
 options = ['RelTol',1e-3];

 [t1,y1] = ode15s(@Model_Publication,tspan,y0,options,params);

% Run Main simulation
 params{1,1}(1,1)=1;

 [t2,y2] = ode15s(@Model_Publication,tspan,y1(end,:),options,params);

%%%%Amplification of genes
params{3}(93:end)=0.5*(y2(end,93:end)./y1(end,93:end)).^(Beta);
params{3}(params{3}<0.001)=0.001;
params{3}(params{3}>1)=1;     
% 

[t3,y3]= ode15s(@Model_Publication,tspan,y2(end,:),options,params);


params{3}(36)=0.8;
[t4,y4]= ode15s(@Model_Publication,tspan,y3(end,:),options,params);


end