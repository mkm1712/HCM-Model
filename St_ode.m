function [TT,Y]=St_ode(params,C0,eps_val,tspan,dt)
% Y is the dynamics of single switch when noise in EGF and species . 
% tau: a vector consisting time scales of species
% C0:  initial state
% eps_val: a vector indicating amplitudes of noise in stimulus and species
% T: simualtion time
% dt: time step size

TT=0:dt:(tspan+dt);  
n=length(TT)-1; 
eps_con=[eps_val*ones(86,1);5*eps_val*ones(1084,1)];    %% Magnitude of noise constant for signaling and genes
% eps_con=eps_val;
nvar=length(C0);
rand_d=normrnd(0,1,nvar+1,n)*sqrt(dt);   %% delta Wn
Y=ones(n+1,nvar);  
tau=params{1,2}';
Y(1,:)=C0';   

for j=2:n+1
    y=Y(j-1,:);   
eta=zeros(nvar,1);
for ii=1:nvar
    eta(ii)  = eta(ii)   -dt*eta(ii)   +eps_con(ii) *rand_d(ii,j-1) + 0.5*eps_con(ii)^2*(rand_d(ii,j-1).^2-dt);
end
% Excluding Model inputs
eta(1)=0;
eta(13)=0;
eta(14)=0;
eta(35)=0;
%
dy=Model_Publicationnoise(y,params);
dy=dy+eta./tau;
    temp=Y(j-1,:)'+dy*dt;
    Y(j,:)=((abs(temp)+temp)/2)';  
    
end

end
