clear all;clc;close all;

M = 16;                 % Modulation order
k = log2(M);            % Bits per symbol
EbNoVec = (0:15)';      % Eb/No values (dB)
numSymPerFrame = 100;   % Number of QAM symbols per frame

berEst = zeros(size(EbNoVec));
berEst_gr = zeros(size(EbNoVec)); 
for n = 1:length(EbNoVec)
    % Convert Eb/No to SNR
    snrdB = EbNoVec(n) + 10*log10(k);
    % Reset the error and bit counters
    numErrs = 0;
    numErrs_gr = 0;
    numBits = 0;
    
    while numErrs < 200 && numBits < 1e7
        % Generate binary data and convert to symbols
        dataIn = randi([0 1],numSymPerFrame,k);
        dataSym = bi2de(dataIn);
        
        % QAM modulate using 'binary' and 'Gray' symbol mapping
        txSig = qammod(dataSym,M,'bin');
        txSig_gr = qammod(dataSym,M);
        
        % Pass through AWGN channel
        rxSig = awgn(txSig,snrdB,'measured');
        rxSig_gr = awgn(txSig_gr,snrdB,'measured');

        % Demodulate the noisy signal
        rxSym = qamdemod(rxSig,M,'bin');
        rxSym_gr = qamdemod(rxSig_gr,M);

        % Convert received symbols to bits
        dataOut = de2bi(rxSym,k);
        dataOut_gr = de2bi(rxSym_gr,k);

        % Calculate the number of bit errors
        nErrors = biterr(dataIn,dataOut);
        nErrors_gr = biterr(dataIn,dataOut_gr);

        % Increment the error and bit counters
        numErrs = numErrs + nErrors;
        numErrs_gr = numErrs_gr + nErrors_gr;
        numBits = numBits + numSymPerFrame*k;
    end
    
    % Estimate the BER
    berEst(n) = numErrs/numBits;
    berEst_gr(n) = numErrs_gr/numBits;
end

berTheory = berawgn(EbNoVec,'qam',M);

semilogy(EbNoVec,berTheory,'b-*')
hold on
semilogy(EbNoVec,berEst,'k-*')
hold on
semilogy(EbNoVec,berEst_gr,'r-*')
grid
title('BER Characteristics for 16QAM')
legend('Theoritical BER','Estimated BER without gc','Estimated BER with gc')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
