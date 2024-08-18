clc;
close all;


nSym=10^4; %No. of OFDM symbols to transmit
EbN0dB = 0:2:10;

%----OFDM Parameters
N=64;%Total number of subcarriers
Nds=48;%Numberof data subcarriers
Nps=4;%Number of pilot subcarriers
ofdmBW = 20*10^6; %OFDM BW

%----Derived Parameters
deltaF=ofdmBW/N;
Tfft=1/deltaF;
Tgi=Tfft/16;
Tsignal=Tgi+Tfft;
Ncp=N*Tgi/Tfft;
Nusc=Nds+Nps;
nBitsPerSum=Nusc;

EsN0dB = EbN0dB+10*log10(Nusc/N)+10*log10(N/(Ncp+N));

for i=1:length(EsN0dB)
    errors(i)=0;
    for j=1:nSym

        s=2*round(rand(1,Nusc))-1;

        X_Freq=[zeros(1,1) s(1:Nusc/2) zeros(1,11) s(Nusc/2+1:end)];

        %IFFT block
        x_Time=N/sqrt(Nusc)*ifft(X_Freq);

        %Adding cyclic prefix
        ofdm_signal=[x_Time(N-Ncp+1:N) x_Time];

        %----AWGN Channel modelling----

        noise=1/sqrt(2)*(randn(1,length(ofdm_signal))+1i*randn(1,length(ofdm_signal)));
        r = sqrt((N+Ncp)/N)*ofdm_signal+10^(-EsN0dB(i)/20)*noise;

        %Removing cyclic prefix
        r_Parallel=r(Ncp+1:(N+Ncp));

        %FFT block
        r_Time=sqrt(Nusc)/N*(fft(r_Parallel));

        %Extracting data carriers from FFT output
        R_Freq =r_Time([(2:Nusc/2+1) (Nusc/2+13:Nusc+12)]);

        %BPSK Demodulation
        R_Freq(R_Freq>0) =+1;
        R_Freq(R_Freq<0) =-1;
        s_cap=R_Freq;

        %Error Calculation
        numErrors=sum(abs(s_cap-s)/2);
        errors(i)=errors(i)+numErrors;
    end

    %Theoritical Error
    theoriticalBER(i)=(1/2)*erfc(sqrt(10.^(EbN0dB(i)/10)));
end

simulatedBER=errors/(nSym*Nusc);
figure();
semilogy(EbN0dB,(simulatedBER),'r*');
hold on;
semilogy(EbN0dB,(theoriticalBER),'k');
grid on;
title('BER vs Eb/N0dB for OFDM with BPSK modulation over AWGN');
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('simulated','theoritical');
savefig('qq');
