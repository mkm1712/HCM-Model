%% Run
tic

% Loops for generating in-sicico populations
for nn=1:50
disp(nn);
% Resest parameters in each loop
params=[];y0=[];y1=[];y2=[]; y3=[];
[params,y0] = Model_Publication_loadParams();

% Setting Run parameters
dt=0.01; 
eps_val=0.05;
tspan = 200; 

 % Run initial simulation
[t1,y1]=St_ode(params,y0,eps_val,tspan,dt);

% Run Main simulation
params{1,1}(1,1)=1;
[t2,y2]=St_ode(params,y1(end,:),eps_val,tspan,dt);
y2=real(y2);

% Amplification of genes
Beta=2;
params{3}(93:end)=0.5*(y2(end,93:end)./y1(end,93:end)).^(Beta);
params{3}(params{3}<0.001)=0.001;
params{3}(params{3}>1)=1;     
[t3,y3]=St_ode(params,y2(end,:),eps_val,tspan,dt);

%ATP low simulation
params{3}(36)=0.8; % Applying ATP low
[t4,y4]=St_ode(params,y3(end,:),eps_val,tspan,dt);

Sig_treatment (:,nn)=real(y4(end,1:86))';
Sig_control (:,nn)=real(y1(end,1:86))';
Genes_treatment(:,nn)=real(y4(end,93:end))';
Genes_control (:,nn)=real(y1(end,93:end))';
HYP_C(nn)=real(y1(end,89));
HYP_T(nn)=real(y4(end,89));
APO_C(nn)=real(y1(end,92));
APO_T(nn)=real(y4(end,92));
end
toc
%%
% Data preparation
[~, TF1]=rmmissing(Genes_control,2);
[~, TF2]=rmmissing(Genes_treatment,2);
Genes_control(:,TF1==1)=nan;
Genes_control(:,TF2==1)=nan;
Genes_treatment(:,TF1==1)=nan;
Genes_treatment(:,TF2==1)=nan;
Genes_control=rmmissing(Genes_control,2);
Genes_treatment=rmmissing(Genes_treatment,2);

% Data Imputation
for ii=1:1078
Genes_control(ii,Genes_control(ii,:)==0)=mean(Genes_control(ii,:));
Genes_treatment(ii,Genes_treatment(ii,:)==0)=mean(Genes_treatment(ii,:));
end
%%
load ('Run-Main-Paper') % Model Data used for Fig 2
%%% FC Calculation
BB1=Genes_treatment./Genes_control;
BB1(BB1==inf)=nan;
for ll=1:1078
BB2= rmoutliers(BB1(ll,:),'quartiles');
BB3(ll)=mean(BB2,'omitnan');
end

% Deterministic solution
[params,y0] = Model_Publication_loadParams();
Beta=2;
[tt1,yy1,~,~,tt4,yy4]=Run_hybrid(params,y0,Beta);
SS=real(yy4(end,1:86)./yy1(end,1:86));
AA=real(yy4(end,93:end)./yy1(end,93:end));
%%
% T-test genes
P = mattest(Genes_treatment,Genes_control);
PSS= mattest(Sig_treatment,Sig_control);
% FDR
[h_fdr, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P,0.05,'pdep','yes');
[h_fdrSS, crit_pSS, adj_ci_cvrgSS, adj_pSS]=fdr_bh(PSS,0.05,'pdep','yes');
FC_Gene=log2(BB3);
% FC_Gene=log2(AA); %Value of deterministc FC will be used for analysis and validation in the paper
FC_SS=real(log2(SS)');
h_fdr=h_fdr';
adj_p=adj_p';


%%%% Volcano Plot
figure
H=h_fdr;
Pv=adj_p;
for ii=1:1078
if FC_Gene(ii)<-0.0 && H(ii)>0.5
    HH(ii,:)=[0 0 1];
elseif FC_Gene(ii)>0.0 && H(ii)>0.5
    HH(ii,:)=[1 0 0];
else HH(ii,:)=[0.5 0.5 0.5];
end
end
sz=10;
scatter (FC_Gene,-log10(Pv),sz,HH,'filled')
xlim([-0.4,1.2]);
ylim([0,18]);
xlabel('LOG2 (FC)')
ylabel('-LOG10 (adjusted p-value)')
ADJJ=adj_p';
HHFDR=h_fdr';
FCC=FC_Gene';
%% PCA
Data=[Genes_control';Genes_treatment'];
Pheno={};
for ii=1:(size(Genes_control,2)+size(Genes_treatment,2))
    if ii<=size(Genes_control,2)  
Pheno{end+1}=['Normal', num2str(ii)];
    else
Pheno{end+1}=['HCM', num2str(ii)];
    end
end
Pheno=Pheno';
mapcaplot(Data, Pheno)

%%
%%% Calculating Signal to Noise Ratio and Cofficient of Varitation
for gg=1:1078
SNR(gg)=mean(Genes_treatment(gg,:))./std(Genes_treatment(gg,:));
CV(gg)=std(Genes_treatment(gg,:))./mean(Genes_treatment(gg,:));
end
figure
plot(SNR,'.');
figure
histogram(SNR);
figure
boxplot(CV,'Labels','Model','Widths',1);
%% DEG names
for ii=1:1078
if -log10(Pv(ii))>12 && abs(FC_Gene(ii))>0.2
   text(FC_Gene(ii),-log10(Pv(ii)),params{1,4}(ii+89));
end
end
%% Violin plot
YY=Genes_treatment(90:100,:)';
[hh,L,MX,MED]=violin(YY,'xlabel',params{1,4}(90:100));

