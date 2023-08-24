%%for distortion noise wd
Pmax=20E-3; %% maximal optical power [W] 
Pmin=0;   %% minimal optical power [W]
wave=450E-9; %% wavelength
PDE=0.35; %% PDE
Na=8192; %% number of SPADs in array
Td=2E-9; %% dead time

PB=10E-9; %% background power [W]
DCR=0.5E6; %% dark count rate
AP=0.75E-2; %% afterpulsing
CT=2.5E-2; %% crosstalk
PaveR=10.^(-9:0.25:-1); %% average received power [W]

Eph=6.6E-34*3E8/wave; %% photon energy

for kt=3:5
    kb=-kt;
delta=(Pmax-Pmin)./(kt-kb); %% scaling factor
Pbias=(Pmin.*kt-Pmax.*kb)./(kt-kb); %% bias
fkb=1./sqrt(2.*pi).*exp(-kb.^2./2); %% f(kappa b)
fkt=1./sqrt(2.*pi).*exp(-kt.^2./2); %% f(kappa t)
vB=PDE.*PB./Eph; %% background rate
PaveT=delta.*(fkb-fkt+kt.*qfunc(kt)+kb.*qfunc(-kb))+Pbias; %% average transmitted power [W]
pathloss=PaveR./PaveT; %% channel path loss
Cs=PDE.*pathloss.*(1+AP+CT)./Eph; 
Cn=(DCR+vB).*(1+AP+CT);
phi1=Cs.*delta;
phi2=Cs.*Pbias+Cn;

%% OFDM setting
Nfft=1024;  %% size of FFT operation
N_sc=Nfft/2-1; %% number of subcarrier
frame=1000;  %% OFDM frames
M=16; %% Modulation order
k=log2(M); %% bits/symbol
Ts=40E-9;  %% symobl time
pow_ctrl_fac=Nfft./(Nfft-2).*ones(N_sc,1);  %% power control factor; Uniform power allocation
symbol_rate=1/Ts; 

flag=1; %% flag indicating if do simulation or not
if flag==1
%% simulation
[x,xc,QAM_sym_input_mtx,bit_input]=OFDM_signal_generation(N_sc,frame,k,M,pow_ctrl_fac,Nfft,kt,kb,0,0); %% 
xt=delta.*xc+Pbias; %% transmitted optical power

% Opti_pow_Rx=xt.*pathloss;
% plot(Opti_pow_Rx(1:1000))
% hold on
% plot(PaveR.*ones(size(Opti_pow_Rx(1:1000))),'-k','linewidth',1.5)

ua_1_mom_est=zeros(1,length(pathloss));
ua_2_mom_est=zeros(1,length(pathloss));
alpha_est_=zeros(1,length(pathloss));
BER_est_mean=inf(1,length(pathloss));
ws_variance=zeros(1,length(pathloss));
for n=1:1:length(pathloss)
    tic
lama=Cs(n).*xt+Cn;
% ua=lama.*Ts.*exp(-lama.*Td./Na);
ua=(lama.*Ts)./(1+lama.*Td./Na);
% vara=ua-lama.^2.*Ts.*Td./Na.*exp(-2.*lama.*Td./Na).*(2-Td./Ts);
vara=(Na./(1+lama.*Td./Na).^3).*(lama.*Ts./Na+(1/6).*(1./(1+lama.*Td./Na)).*((lama.*Td./Na).^2).*(6+4*lama.*Td./Na+((lama.*Td./Na).^2)));
% vara=Na.*(ex_variance_1(lama./Na,Ts,Td))';
ws=normrnd(0,sqrt(vara));
y=ua+ws; %% SPAD output signal
alpha_est_(n)=mean(ua.*x);
ua_1_mom_est(n)=mean(ua); 
ua_2_mom_est(n)=mean(ua.^2);
% alpha_est_(n)=alpha_est(y,x);
ws_variance(n) = mean(vara);
[SNR_vec,H_est_vec,noise_var_vec,fft_sc_out_mtx]=OFDM_signal_decoding(y,QAM_sym_input_mtx,frame,Nfft,0,0); %% fft_sc_out_mtx is the signal after fft
f_vec=(1:Nfft/2-1)./Nfft.*symbol_rate; %% frequency vector
% SNR_est(n)=mean(SNR_vec);

AfterZF=fft_sc_out_mtx./kron(H_est_vec,ones(1,frame)); %% signal after zero-forcing equalization
avg_sym_amp=sqrt(mean(abs(qammod(0:M-1,M,'gray')).^2)); %% average symbol amplitude
demod_mtx=qamdemod(avg_sym_amp.*AfterZF,M,'gray'); %% demodulation
bit_output=de2bi(demod_mtx);
total_bit_num=length(bit_output(:));
error_vec=sum(bit_output~=bit_input,2); %% error vector
total_error_num=sum(error_vec); %% total error number
error_mtx=reshape(error_vec,N_sc,frame); %% error matrix
error_sc=sum(error_mtx,2); %% error number for each subcarrier
BER_est_sc=error_sc./(frame.*k); %% the BER for each subcarrier
BER_est_mean(n)=total_error_num/total_bit_num; %% estimated average BER
toc
end

% ddd = zeros(1, length(pathloss));
% ddd1 = zeros(1, length(pathloss));
% ddd2 = zeros(1, length(pathloss));
% dddd = zeros(1, length(pathloss));
for cc=1:1:length(pathloss)
fun = @(z) (1/sqrt(2*pi)).*exp(-(z.^2)/2).*(((((phi1(cc).*z+phi2(cc))).*Ts)./(1+(phi1(cc).*z+phi2(cc)).*Td./Na)).^2);
ddd(kt-2,cc)=integral(fun,kb,kt);
ddd1(kt-2,cc)=(((((phi1(cc).*kb+phi2(cc))).*Ts)./(1+(phi1(cc).*kb+phi2(cc)).*Td./Na)).^2);
ddd2(kt-2,cc)=(((((phi1(cc).*kt+phi2(cc))).*Ts)./(1+(phi1(cc).*kt+phi2(cc)).*Td./Na)).^2);
dddd(kt-2,cc)=ddd(cc)+ddd1(cc)*qfunc(-kb)+ddd2(cc)*qfunc(kt);
end
% figure;
% plot(dddd(kt-2,:),'*g');
% hold on;
% plot(dddd(2,:),'*b');
% % hold on;
% plot(dddd,'*r');
% hold on;
% plot(ua_2_mom_est);
end
var_wd_numerical(kt-2,:)=dddd(kt-2,:)-(ua_1_mom(kt-2,:).^2)-(alpha_ana(kt-2,:).^2);
SDNR(kt-2,:) = ((alpha_ana(kt-2,:).^2).*(Nfft/(Nfft-2)))./var_wd_numerical(kt-2,:);
end

figure;
semilogx(PaveR,pow2db(SDNR(1,:)),'g','linewidth',0.8);
xlabel('average received optical power $\overline{P}_{\rm Rx}$ [W]','interpreter','latex','fontsize',12)
ylabel('SDNR [dB]','interpreter','latex','fontsize',12)
hold on;
semilogx(PaveR,pow2db(SDNR(2,:)),'b','linewidth',0.8);
hold on;
semilogx(PaveR,pow2db(SDNR(3,:)),'r','linewidth',0.8);
hold on;
% fun = @(z) (1/sqrt(2*pi)).*exp(-(z.^2)/2).*(((((phi1(cc).*z+phi2(cc))).*Ts)./(1+(phi1(cc).*z+phi2(cc)).*Td./Na)).^2);
% ddd(kt-2,cc)=integral(fun,kb,kt);
% ddd1(kt-2,cc)=(((((phi1(cc).*kb+phi2(cc))).*Ts)./(1+(phi1(cc).*kb+phi2(cc)).*Td./Na)).^2);
% ddd2(kt-2,cc)=(((((phi1(cc).*kt+phi2(cc))).*Ts)./(1+(phi1(cc).*kt+phi2(cc)).*Td./Na)).^2);
% dddd(kt-2,cc)=ddd(cc)+ddd1(cc)*qfunc(-kb)+ddd2(cc)*qfunc(kt);