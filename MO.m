clear all;clc;close all;

N=10000; % Number of bits to be transmitted
r=randi([0,1],1,N);   % Generating random bits consisting of 1 and 0 of form 1*N
BER_QPSKgc=[]; SNR_QPSKgc=[];  % These matrix stores the SNR and BER for QPSK with Grey Coding
BER_QPSKwgc=[]; SNR_QPSKwgc=[];% These matrix stores the SNR and BER for QPSK without Grey Coding
BER_QAM=[]; SNR_QAM=[];        % These matrix stores the SNR and BER for 16-QAM

                                    %% QPSK (with gray coding)
% QPSK Modulation
QPSKgc=[];
for i=1:2:length(r)
    if r(i)==0 && r(i+1)==0
        Y=cosd(225)+1j*sind(225);
    elseif r(i)==0 && r(i+1)==1
        Y=cosd(135)+1j*sind(135);
    elseif r(i)==1 && r(i+1)==0
        Y=cosd(315)+1j*sind(315);
    elseif r(i)==1 && r(i+1)==1
        Y=cosd(45)+1j*sind(45);
    end
    
QPSKgc=[QPSKgc Y];
end

for snrdb=1:0.5:20
    QPSKawgn=awgn(complex(QPSKgc),snrdb);   %% Sending QPSK over AWGN Channel
    % Defining Threshold for hamming criteria (Minimum Distance)
    trans=[0.707+1j*0.707,-0.707+1j*0.707,-0.707-1j*0.707,0.707-1j*0.707];
    detgr=[];

    for mm=1:length(QPSKgc)
        for nn=1:length(trans)
            a=(real(QPSKawgn(mm))-real(trans(nn)))^2;
            b=(imag(QPSKawgn(mm))-imag(trans(nn)))^2;
            errorgr(nn)=sqrt(a+b);
        end
        iden=trans(find(errorgr==min(errorgr)));
        detgr=[detgr iden];
    end
    % Demodulation 
    x_capgr=[];
    for k=1:length(detgr)
    if real(detgr(k))==0.707 && imag(detgr(k))==0.707
        d=[1 1];
    elseif real(detgr(k))==-0.707 && imag(detgr(k))==0.707
        d=[0 1];
    elseif real(detgr(k))==-0.707 && imag(detgr(k))==-0.707
        d=[0 0];
    elseif real(detgr(k))==0.707 && imag(detgr(k))==-0.707
        d=[1 0];
    end
    x_capgr=[x_capgr d];
    end
% Calculation of BER
[no_of_faulty_bits BERR] = biterr(r,x_capgr);

BER_QPSKgc=[BER_QPSKgc BERR];
SNR_QPSKgc=[SNR_QPSKgc snrdb];
end

Th_BER_QPSKgc=0.5*erfc(sqrt(10.^(SNR_QPSKgc/10)));  %theoritical value
figure();
semilogy(SNR_QPSKgc,BER_QPSKgc,'k+',SNR_QPSKgc,Th_BER_QPSKgc,'p-')
xlim([0 10]);
legend('QPSK Simulation','QPSK Theoratical');

                                    %% QPSK (without gray coding)
% QPSK Modulation
QPSKwgc=[];
for i=1:2:length(r)
    if r(i)==0 && r(i+1)==0
        Y=cosd(0)+1j*sind(0);
    elseif r(i)==0 && r(i+1)==1
        Y=-cosd(0)+1j*sind(0);
    elseif r(i)==1 && r(i+1)==0
        Y=cosd(90)+1j*sind(90);
    elseif r(i)==1 && r(i+1)==1
        Y=cosd(90)-1j*sind(90);
    end
    
QPSKwgc=[QPSKwgc Y];
end
for snrdb=1:0.5:20
    QPSKawgn=awgn(complex(QPSKwgc),snrdb);   %% Sending QPSK over AWGN Channel
    % Defining Threshold for hamming criteria (Minimum Distance)
    trans=[cosd(0)+1j*sind(0),-cosd(0)+1j*sind(0),cosd(90)+1j*sind(90),cosd(90)-1j*sind(90)];

    det=[];
    for mm=1:length(QPSKwgc)
        for nn=1:length(trans)
            a=(real(QPSKawgn(mm))-real(trans(nn)))^2;
            b=(imag(QPSKawgn(mm))-imag(trans(nn)))^2;
            error(nn)=sqrt(a+b);
        end
        iden=trans(find(error==min(error)));
        det=[det iden];
    end
    % Demodulation
    x_cap=[];
    for k=1:length(det)
    if real(det(k))==cosd(90) && imag(det(k))==-sind(90)
        d=[1 1];
    elseif real(det(k))==-cosd(0) && imag(det(k))==sind(0)
        d=[0 1];
    elseif real(det(k))==cosd(0) && imag(det(k))==sind(0)
        d=[0 0];
    elseif real(det(k))==cosd(90) && imag(det(k))==sind(90)
        d=[1 0];
    end
    x_cap=[x_cap d];
    end
% Calculation of BER
[no_of_faulty_bits BERR] = biterr(r,x_cap);
BER_QPSKwgc=[BER_QPSKwgc BERR];
SNR_QPSKwgc=[SNR_QPSKwgc snrdb];
end
%% Graph of SNR vs BER for QPSK with and without grey code
figure()
semilogy(SNR_QPSKgc,BER_QPSKgc,'ko-',SNR_QPSKwgc,BER_QPSKwgc,'r-d')
title('BER Characteristics of QPSK');
xlabel('SNR(dB)');ylabel('BER');
legend('Grey Coding','Without Grey Coding');

                                    %% 16QAM
% 16-QAM Modulation
QAM=[];
for i=1:4:length(r)
    if r(i)==0 && r(i+1)==0 && r(i+2)==1 && r(i+3)==0
        Y=sqrt(1/10)*(-3+1j*3);
    elseif r(i)==0 && r(i+1)==1 && r(i+2)==1 && r(i+3)==0
        Y=sqrt(1/10)*(-1+1j*3);
    elseif r(i)==0 && r(i+1)==0 && r(i+2)==1 && r(i+3)==1
        Y=sqrt(1/10)*(-3+1j*1);
    elseif r(i)==0 && r(i+1)==1 && r(i+2)==1 && r(i+3)==1
        Y=sqrt(1/10)*(-1+1j*1);
    elseif r(i)==1 && r(i+1)==1 && r(i+2)==1 && r(i+3)==0
        Y=sqrt(1/10)*(1+1j*3);
    elseif r(i)==1 && r(i+1)==0 && r(i+2)==1 && r(i+3)==0
        Y=sqrt(1/10)*(3+1j*3);
    elseif r(i)==1 && r(i+1)==1 && r(i+2)==1 && r(i+3)==1
        Y=sqrt(1/10)*(1+1j*1);
    elseif r(i)==1 && r(i+1)==0 && r(i+2)==1 && r(i+3)==1
        Y=sqrt(1/10)*(3+1j*1);
    elseif r(i)==0 && r(i+1)==0 && r(i+2)==0 && r(i+3)==1
        Y=sqrt(1/10)*(-3-1j*1);
    elseif r(i)==0 && r(i+1)==1 && r(i+2)==0 && r(i+3)==1
        Y=sqrt(1/10)*(-1-1j*1);
    elseif r(i)==0 && r(i+1)==0 && r(i+2)==0 && r(i+3)==0
        Y=sqrt(1/10)*(-3-1j*3);
    elseif r(i)==0 && r(i+1)==1 && r(i+2)==0 && r(i+3)==0
        Y=sqrt(1/10)*(-1-1j*3);
    elseif r(i)==1 && r(i+1)==1 && r(i+2)==0 && r(i+3)==1
        Y=sqrt(1/10)*(1-1j*1);
    elseif r(i)==1 && r(i+1)==0 && r(i+2)==0 && r(i+3)==1
        Y=sqrt(1/10)*(3-1j*1);
    elseif r(i)==1 && r(i+1)==1 && r(i+2)==0 && r(i+3)==0
        Y=sqrt(1/10)*(1-1j*3);
    elseif r(i)==1 && r(i+1)==0 && r(i+2)==0 && r(i+3)==0
        Y=sqrt(1/10)*(3-1j*3);
    end
    
QAM=[QAM Y];
end

for snrdb=1:0.5:20
    QAMawgn=awgn(complex(QAM),snrdb); %% Sending 16-QAM over AWGN Channel
    % Defining Threshold for hamming criteria (Minimum Distance)
    trans=sqrt(1/10)*[(-3+1j*3),(-1+1j*3),(-3+1j*1),(-1+1j*1),(1+1j*3),(3+1j*3),(1+1j*1),(3+1j*1),(-3-1j*1),(-1-1j*1),(-3-1j*3),(-1-1j*3),(1-1j*1),(3-1j*1),(1-1j*3),(3-1j*3)];
    det_qam=[];

    for mm=1:length(QAM)
        for nn=1:length(trans)
            a=(real(QAMawgn(mm))-real(trans(nn)))^2;
            b=(imag(QAMawgn(mm))-imag(trans(nn)))^2;
            error_qam(nn)=sqrt(a+b);
        end
        iden=trans(find(error_qam==min(error_qam)));
        det_qam=[det_qam iden];
    end

    % Demodulation
    x_qam=[];
    for k=1:length(det_qam)
        if real(det_qam(k))==-3*sqrt(1/10) && imag(det_qam(k))==3*sqrt(1/10)
            d=[0 0 1 0];
        elseif real(det_qam(k))==-1*sqrt(1/10) && imag(det_qam(k))==3*sqrt(1/10)
            d=[0 1 1 0];
        elseif real(det_qam(k))==-3*sqrt(1/10) && imag(det_qam(k))==1*sqrt(1/10)
            d=[0 0 1 1];
        elseif real(det_qam(k))==-1*sqrt(1/10) && imag(det_qam(k))==1*sqrt(1/10)
            d=[0 1 1 1];
        elseif real(det_qam(k))==1*sqrt(1/10) && imag(det_qam(k))==3*sqrt(1/10)
            d=[1 1 1 0];
        elseif real(det_qam(k))==3*sqrt(1/10) && imag(det_qam(k))==3*sqrt(1/10)
            d=[1 0 1 0];
        elseif real(det_qam(k))==1*sqrt(1/10) && imag(det_qam(k))==1*sqrt(1/10)
            d=[1 1 1 1];
        elseif real(det_qam(k))==3*sqrt(1/10) && imag(det_qam(k))==1*sqrt(1/10)
            d=[1 0 1 1];
        elseif real(det_qam(k))==-3*sqrt(1/10) && imag(det_qam(k))==-1*sqrt(1/10)
            d=[0 0 0 1];
        elseif real(det_qam(k))==-1*sqrt(1/10) && imag(det_qam(k))==-1*sqrt(1/10)
            d=[0 1 0 1];
        elseif real(det_qam(k))==-3*sqrt(1/10) && imag(det_qam(k))==-3*sqrt(1/10)
            d=[0 0 0 0];
        elseif real(det_qam(k))==-1*sqrt(1/10) && imag(det_qam(k))==-3*sqrt(1/10)
            d=[0 1 0 0];
        elseif real(det_qam(k))==1*sqrt(1/10) && imag(det_qam(k))==-1*sqrt(1/10)
            d=[1 1 0 1];
        elseif real(det_qam(k))==3*sqrt(1/10) && imag(det_qam(k))==-1*sqrt(1/10)
            d=[1 0 0 1];
        elseif real(det_qam(k))==1*sqrt(1/10) && imag(det_qam(k))==-3*sqrt(1/10)
            d=[1 1 0 0];
        elseif real(det_qam(k))==3*sqrt(1/10) && imag(det_qam(k))==-3*sqrt(1/10)
            d=[1 0 0 0];
        end
        x_qam=[x_qam d];
    end
% Calculation of BER
[no_of_faulty_bits QAM_BER] = biterr(r,x_qam);
BER_QAM=[BER_QAM QAM_BER];
SNR_QAM=[SNR_QAM snrdb];
end

figure()
semilogy(SNR_QAM,BER_QAM,'k-d')
title('BER Characteristics of 16-QAM');
xlabel('SNR(dB)');ylabel('BER');

%% Comparison of Different Modulation Schemes
figure()
semilogy(SNR_QPSKwgc,BER_QPSKwgc,'bo-',SNR_QPSKgc,BER_QPSKgc,'rx-',SNR_QAM,BER_QAM,'k-d')
title('BER Characteristics of different Modulation Schemes');
xlabel('SNR(dB)');ylabel('BER');
legend('QPSK without Gray','QPSK Gray','16QAM');