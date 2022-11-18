clear all;clc;close all;
%%----------------------------------QPSK-------------------------------------
%number of symbols in simulation
Nsyms = 1e6;
% energy per symbol
Es = 1;
% energy per bit (2 bits/symbol for QPSK)
Eb = Es / 2;
% Eb/No values to simulate at, in dB
EbNo_dB = linspace(0, 10, 11);
% Eb/No values in linear scale
EbNo_lin = 10.^(EbNo_dB / 10);
% keep track of bit errors for each Eb/No point
bit_err = zeros(size(EbNo_lin));
for i=1:length(EbNo_lin)
    % generate source symbols
    syms = (1 - 2 * (randn(Nsyms,1) > 0)) + j * (1 - 2 * (randn(Nsyms, 1) > 0));
    % add noise
    syms_noisy = sqrt(Es/2) * syms + sqrt(Eb/(2*EbNo_lin(i))) * (randn(size(syms)) + j * randn(size(syms)));
    % recover symbols from each component (real and imaginary)
    syms_rec_r = sign(real(syms_noisy));
    syms_rec_i = sign(imag(syms_noisy));
    % count bit errors
    bit_err(i) = sum((syms_rec_r ~= real(syms)) + (syms_rec_i ~= imag(syms)));
end
% convert to bit error rate
bit_err = bit_err / (2 * Nsyms);

% calculate theoretical bit error rate, functionally equivalent to:
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));
figure;
semilogy(EbNo_dB, bit_err, 'rx', EbNo_dB, bit_err_theo, 'b-');
xlabel('Eb/No (dB)');
ylabel('Bit error rate');
title('QPSK bit error rate');
legend('Simulation','Theory');
grid on;

%%-------------------------------16QAM----------------------------------------------------- 
% Bit Error Rate for 16-QAM modulation using Gray modulation mapping
N = 10^5; % number of symbols
M = 16;   % constellation size
k = log2(M); % bits per symbol

% defining the real and imaginary PAM constellation
% for 16-QAM
alphaRe = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
alphaIm = [-(2*sqrt(M)/2-1):2:-1 1:2:2*sqrt(M)/2-1];
k_16QAM = 1/sqrt(10);

Eb_N0_dB  = [0:15]; % multiple Es/N0 values
Es_N0_dB  = Eb_N0_dB + 10*log10(k);

% Mapping for binary <--> Gray code conversion
ref = [0:k-1];
map = bitxor(ref,floor(ref/2));
[tt ind] = sort(map);                                

for ii = 1:length(Eb_N0_dB)
    
    % symbol generation
    % ------------------
    ipBit = rand(1,N*k,1)>0.5; % random 1's and 0's
    ipBitReshape = reshape(ipBit,k,N).';
    bin2DecMatrix = ones(N,1)*(2.^[(k/2-1):-1:0]) ; % conversion from binary to decimal
    % real
    ipBitRe =  ipBitReshape(:,[1:k/2]);
    ipDecRe = sum(ipBitRe.*bin2DecMatrix,2);
    ipGrayDecRe = bitxor(ipDecRe,floor(ipDecRe/2));
    % imaginary
    ipBitIm =  ipBitReshape(:,[k/2+1:k]);
    ipDecIm = sum(ipBitIm.*bin2DecMatrix,2);
    ipGrayDecIm = bitxor(ipDecIm,floor(ipDecIm/2)); 
    % mapping the Gray coded symbols into constellation
    modRe = alphaRe(ipGrayDecRe+1);
    modIm = alphaIm(ipGrayDecIm+1);
    % complex constellation
    mod = modRe + j*modIm;
    s = k_16QAM*mod; % normalization of transmit power to one 
    
    % noise
    % -----
    n = 1/sqrt(2)*[randn(1,N) + j*randn(1,N)]; % white guassian noise, 0dB variance 
    
    y = s + 10^(-Es_N0_dB(ii)/20)*n; % additive white gaussian noise

    % demodulation
    % ------------
    y_re = real(y)/k_16QAM; % real part
    y_im = imag(y)/k_16QAM; % imaginary part

    % rounding to the nearest alphabet
    ipHatRe = 2*floor(y_re/2)+1;
    ipHatRe(find(ipHatRe>max(alphaRe))) = max(alphaRe);
    ipHatRe(find(ipHatRe<min(alphaRe))) = min(alphaRe);
    ipHatIm = 2*floor(y_im/2)+1;
    ipHatIm(find(ipHatIm>max(alphaIm))) = max(alphaIm);
    ipHatIm(find(ipHatIm<min(alphaIm))) = min(alphaIm);

    % Constellation to Decimal conversion
    ipDecHatRe = ind(floor((ipHatRe+4)/2+1))-1; % LUT based
    ipDecHatIm = ind(floor((ipHatIm+4)/2+1))-1; % LUT based

    % converting to binary string
    ipBinHatRe = dec2bin(ipDecHatRe,k/2);
    ipBinHatIm = dec2bin(ipDecHatIm,k/2);

    % converting binary string to number
    ipBinHatRe = ipBinHatRe.';
    ipBinHatRe = ipBinHatRe(1:end).';
    ipBinHatRe = reshape(str2num(ipBinHatRe).',k/2,N).' ;
    
    ipBinHatIm = ipBinHatIm.';
    ipBinHatIm = ipBinHatIm(1:end).';
    ipBinHatIm = reshape(str2num(ipBinHatIm).',k/2,N).' ;

    % counting errors for real and imaginary
    nBitErr(ii) = size(find([ipBitRe- ipBinHatRe]),1) + size(find([ipBitIm - ipBinHatIm]),1) ;

end 
simBer = nBitErr/(N*k);
theoryBer = (1/k)*3/2*erfc(sqrt(k*0.1*(10.^(Eb_N0_dB/10))));

figure
semilogy(Eb_N0_dB,theoryBer,'bs-',Eb_N0_dB,simBer,'rx');
axis([0 15 10^-5 1])
grid on
legend('theory', 'simulation');
xlabel('Eb/No, dB')
ylabel('Bit Error Rate')
title('Bit error probability curve for 16-QAM modulation')

%%
figure();
semilogy(EbNo_dB,bit_err,'b-*',Eb_N0_dB,simBer,'r-*');
hold on;
xlabel('SNR (dB)');
ylabel('BER');
legend(' QPSK','16QAM');
