%%for shot noise Ws
qqq = zeros(1, length(pathloss));
qqq1 = zeros(1, length(pathloss));
qqq2 = zeros(1, length(pathloss));
qqqq = zeros(1, length(pathloss));
for cc=1:1:length(pathloss)
%fun=@(z) (1/sqrt(2*pi)).*exp(-(z.^2)/2).*(Na./(1+(phi1(cc).*z+phi2(cc)).*Td./Na).^3).*((phi1(cc).*z+phi2(cc)).*Ts./Na+(1/6).*(1./(1+(phi1(cc).*z+phi2(cc)).*Td./Na)).*(((phi1(cc).*z+phi2(cc)).*Td./Na).^2).*(6+4*(phi1(cc).*z+phi2(cc)).*Td./Na+(((phi1(cc).*z+phi2(cc)).*Td./Na).^2)));
fun = @(z) (1/sqrt(2*pi)).*exp(-(z.^2)/2).*(Na./((1+(((phi1(cc).*z+phi2(cc))./Na).*Td)).^3)).*((((phi1(cc).*z+phi2(cc))./Na).*Ts)+((((phi1(cc).*z+phi2(cc)).*Td)./Na).^2)./(1+((phi1(cc).*z+phi2(cc)).*Td)./Na)+((2/3)*((((phi1(cc).*z+phi2(cc)).*Td)./Na).^3)./(1+((phi1(cc).*z+phi2(cc)).*Td)./Na))+((1/6)*((((phi1(cc).*z+phi2(cc)).*Td)./Na).^4)./(1+((phi1(cc).*z+phi2(cc)).*Td)./Na)));
qqq(cc)=integral(fun,kb,kt);
qqq1(cc)=(Na./((1+(((phi1(cc).*kb+phi2(cc))./Na).*Td)).^3)).*((((phi1(cc).*kb+phi2(cc))./Na).*Ts)+((((phi1(cc).*kb+phi2(cc)).*Td)./Na).^2)./(1+((phi1(cc).*kb+phi2(cc)).*Td)./Na)+((2/3)*((((phi1(cc).*kb+phi2(cc)).*Td)./Na).^3)./(1+((phi1(cc).*kb+phi2(cc)).*Td)./Na))+((1/6)*((((phi1(cc).*kb+phi2(cc)).*Td)./Na).^4)./(1+((phi1(cc).*kb+phi2(cc)).*Td)./Na)));
qqq2(cc)=(Na./((1+(((phi1(cc).*kt+phi2(cc))./Na).*Td)).^3)).*((((phi1(cc).*kt+phi2(cc))./Na).*Ts)+((((phi1(cc).*kt+phi2(cc)).*Td)./Na).^2)./(1+((phi1(cc).*kt+phi2(cc)).*Td)./Na)+((2/3)*((((phi1(cc).*kt+phi2(cc)).*Td)./Na).^3)./(1+((phi1(cc).*kt+phi2(cc)).*Td)./Na))+((1/6)*((((phi1(cc).*kt+phi2(cc)).*Td)./Na).^4)./(1+((phi1(cc).*kt+phi2(cc)).*Td)./Na)));
qqqq(cc)=qqq(cc)+qqq1(cc)*qfunc(-kb)+qqq2(cc)*qfunc(kt);
ssnr_numerical_est(cc)= alpha_est(kt-1).^2.*Nfft./(Nfft-2)./qqqq(cc);
end
%alpha_est(kt-1).^2.*Nfft./(Nfft-2)./ws_mean_var(n)
figure;
plot(ssnr_numerical_est,'*b');
hold on;
plot(SSNR_est(3,:),'g');
hold on;
plot(ws_mean_var,'g');
hold on;
plot(qqqq,'r*');
