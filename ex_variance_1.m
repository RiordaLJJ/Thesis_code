%exact variance
function exact_variance_1 = ex_variance_1(lamda,Tsym_,Td_)
variance = zeros(1,size(lamda,2));
K = floor(Tsym_/Td_);
for i = 1:1:size(lamda,2)
% k = 0:1:K+1;
% Tk = lamda(i)*(Tsym_-k*Td_);
part1 = 0;
for k = 0:1:K
    part2 = 0;
    for j = 0:1:k-1
        part2 = part2 + (k-j)*(((lamda(i)*(Tsym_-k*Td_))^j)*exp(-(lamda(i)*(Tsym_-k*Td_))))/factorial(j);
    end
    part1 = part1 + ((lamda(i)*(Tsym_-k*Td_))-k+part2);
end
c = (lamda(i)*Tsym_)/(1+lamda(i)*Td_);
variance(i) = ((2/(1+lamda(i)*Td_))*part1 - c - c^2);
end
exact_variance_1 = variance;
end