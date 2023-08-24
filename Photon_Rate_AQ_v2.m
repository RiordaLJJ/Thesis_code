%% 
clc
clear
close all

Nspad=1; %% number of SPADs in the array [can be large]
lam0=0E6; %% photon rate when bit 0 is sent
Tsym=40E-9; %% symbol time

% Tsym=1; %% symbol time
data_count=2E4; %% number of symbols [can be large]
pulse=ones(1,data_count); %% binary pulse sequence; all ones means unmodulated signal

Td=2E-9; %% dead time

lam1=10.^(6:0.25:11); %% photon rate when bit 1 is sent

% Photon_Rate_bit1=zeros(size(lam1));
meanPC=zeros(size(lam1)); %% average photon count 
varPC=zeros(size(lam1)); %% variance of photon count
parfor n=1:1:length(lam1)
tic
detected_count_mtx=zeros(Nspad,data_count); %% A matrix which records the detected photon count for each SPAD in each symbol time
for uu=1:Nspad  %% simulate the detected photon counts for each SPAD [use parfar to speed up if needed]
       S=pho_arr_generation(pulse,Tsym,lam1(n),lam0); %% generate the photon arrival using Poisson process 
       [~,RisEdges]=detect_photon_deadtime_AQ(S,Td);  %% get the detected photon arrival considering the dead time effect
       detected_count_mtx(uu,:)=histcounts(RisEdges,0:Tsym:Tsym*length(pulse)); %% get the photon count for the SPAD in each symbol time 
end
N_sum=sum(detected_count_mtx,1); %% get the aggregated detected photon count of the SPAD array in each symbol time
N_sum=N_sum(2:end); %% remove the first symbol which is singular
meanPC(n)=mean(N_sum(pulse(2:end)==1)); % remove the first symbol which is singular
varPC(n)=var(N_sum(pulse(2:end)==1)); % remove the first symbol which is singular
toc
end

mu = lam1.*Tsym;
x = lam1.*Td;
lamda = 1./(1+x);
exact_variance_asymptotic = (lamda.^3).*(mu+(lamda.*(x.^2).*(6+4*x+x.^2))./6);
% exact_variance_asymptotic = lam1.*Tsym./((1+lam1.*Td).^3);
exact_variance = ex_variance_1(lam1,Tsym,Td);
loglog(lam1,meanPC,'bs','linewidth',1.2)
hold on
semilogx(lam1,varPC,'ro','linewidth',1.2)
grid on
xlabel('incident photon rate [cps]')
ylabel('moments of photon count')
hold on
semilogx(lam1,lam1.*Tsym./(1+lam1.*Td),'-b','linewidth',1.2)
hold on
semilogx(lam1,exact_variance,'-r','linewidth',1.2)
hold on
semilogx(lam1,exact_variance_asymptotic,'-g','linewidth',1.2)
hold on
semilogx(lam1,lam1.*Tsym,'--','linewidth',1.2);
legend('mean, simulation','variance, simulation', 'mean, analytical','variance, analytical','variance, analytical asymptotic','Ideal');
% figure;
% plot(lam1,meanPC,'bs','linewidth',1)
% hold on
% plot(lam1,varPC,'ro','linewidth',1)
% grid on
% xlabel('incident photon rate [cps]')
% ylabel('moments of photon count')
% hold on
% plot(lam1,lam1.*Tsym./(1+lam1.*Td),'-b','linewidth',1.2)
% hold on
% plot(lam1,exact_variance,'-r','linewidth',1.2)
% hold on
% plot(lam1,exact_variance_asymptotic,'-g','linewidth',1.2)
% hold on
% % plot(lam1,lam1.*Tsym,'--','linewidth',1.2);
% legend('mean, simulation','variance, simulation', 'mean, analytical','variance, analytical','variance, analytical asymptotic','Ideal');
% 
% print('validations of AQ SPAD moments','-dpng');