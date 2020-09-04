close all
clear all
clc
Sample = 6810;
fc= 20*10^(6); %carrier frequency 
Rs=4*fc;
M = 16;
Tu=224e-6; 
T=Tu/2048;
FS=2048;
%.............................Generating RANDOM DATA...............................

x= mod(rand,30000);
y= mod(rand,30000);
uniform=[];
for i =1:Sample 
    [RandomNumber,x,y]= uniform_A(x,y);
    uniform=[uniform RandomNumber];    
end
%.............................Generating BINARY DATA....................
for i = 1:1:Sample   
    if ( uniform(i)<0.5 )    
        binary_data(i) = 0;       
    else        
        binary_data(i) = 1;   
    end
end


%.............................adding 14 zero...........................
n=15;
holdy = zeros(1,6810);
for i = 1:1:6810
     holdy(n)=  binary_data(i);
     n=n+1;  
end

%.............................binary to decimal.....................
n=1;
for i = 1:4:6824 
     qam_input(n)= 8*holdy(i)+4*holdy(i+1)+2*holdy(i+2)+1*holdy(i+3);
     n=n+1;  
end
%.............................Mapping to 16-QAM constellation.....................

for (i = 1:1:1706)
      if(qam_input(i) == 0)
          mapper(i) = -3+3i;
       elseif(qam_input(i) == 1)
          mapper(i) = -3+1i;
       elseif(qam_input(i) == 2)
          mapper(i) = -3-1i;
       elseif(qam_input(i) == 3)
          mapper(i) = -3-3i;
       elseif(qam_input(i) == 4)
          mapper(i) = -1+3i;
       elseif(qam_input(i) == 5)
          mapper(i) = -1+1i;
       elseif(qam_input(i) == 6)
          mapper(i) = -1-1i;
        elseif(qam_input(i) == 7)
          mapper(i) = -1-3i;
        elseif(qam_input(i) == 8)
          mapper(i) = 1+3i;
        elseif(qam_input(i) == 9)
          mapper(i) = 1+1i;
        elseif(qam_input(i) == 10)
          mapper(i) = 1-1i;
        elseif(qam_input(i) == 11)
           mapper(i) = 1-3i;
        elseif(qam_input(i) == 12)
           mapper(i) = 3+3i;
         elseif(qam_input(i)== 13)
           mapper(i) = 3+1i;
         elseif(qam_input(i)== 14)
           mapper(i) = 3-1i;
          elseif(qam_input(i)== 15)
           mapper(i) = 3-3i;
      end           
end

%    scatterplot(mapper);     

% .................................zero padding..........................
 k=1;
 zero_padding = zeros(1,2048);
for (i = 1:1:853)
         zero_padding(i) = mapper(k);
          k =k+1;
end
 k=854;
for (i = 1196:1:2048)
         zero_padding(i) = mapper(k);
          k =k+1;
end

% .................................IFFT..........................

ifft_data=ifft(zero_padding,2048);
real_ifft = real(ifft_data);
imag_ifft = imag(ifft_data);
ifftt_data = real_ifft+imag_ifft;

% tt=0:T/2:Tu;
% figure(1);
% subplot(211);
% stem(tt(1:20),real_ifft(1:20));
% title('carriers inphase');
% xlabel('Time(s)');
% ylabel('Amplitude');
% subplot(212);
% stem(tt(1:20),imag_ifft(1:20));
% title('carriers Quadrature');
% xlabel('Time(s)');
% ylabel('Amplitude');

% f=(2/T)*(1:(FS))/(FS);
% subplot(211);
% plot(f,abs(fft(ifftt_data ,FS))/FS);
% title('carriers FFT');
% xlabel('Frequency(Hz)');
% ylabel('Amplitude');
% subplot(212);
% pwelch(ifftt_data,[],[],[],2/T);


    EbNo = (0:1:12);
    EsNo= EbNo+10*log10(4)+10*log10(52/44);
    snr=EsNo - 10*log10(64/60);
    ratio= [];
    no_error=[];
    
for kk=1:1:length(snr)
    
    
awgn_noise = add_awgn_noise(zero_padding,snr(kk));
ifft_noise= ifft(awgn_noise);
upscale_awgn =interp(ifft_noise,10);


UPscale_real = interp(real_ifft,10);
UPscale_imag = interp(imag_ifft,10)*1i;
UPscaled = UPscale_real + UPscale_imag;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .................................MODULATION..........................
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  for(i = 1:1:20480)
    t(i) = (1 - i)/(Rs);
  end
   
   for i=1:1:20480
      mod_real(i) = UPscale_real(i)*cos(2*pi*fc*t(i));
      mod_imag(i) = UPscale_imag(i)*sin(2*pi*fc*t(i));
      % ................NOISE.MODULATION
      mod_real_awgn(i) = real(upscale_awgn(i))*cos(2*pi*fc*t(i));
      mod_imag_awgn(i) = imag(upscale_awgn(i))*sin(2*pi*fc*t(i));
   end
   
%    modulated_signal = mod_real + mod_imag;
%    f=(2/T)*(1:(FS))/(FS);
%   subplot(211);
%   plot(f,abs(fft(modulated_signal ,FS))/FS);
%   title('Modulated signal');
%   xlabel('Frequency(Hz)');
%   ylabel('Amplitude');
%   subplot(212);
%   pwelch(modulated_signal,[],[],[],2/T);
% .................................ADDING NOISE CHANNEL..........................
   for i=1:1:20480
      mod_real_noise(i) = mod_real(i)+ mod_real_awgn(i);
      mod_imag_noise(i) = mod_imag(i)+ mod_imag_awgn(i) ;
   end
    modulated = mod_real_noise + mod_imag_noise;
    tt=0:T/2:Tu;
% figure(1);
% subplot(211);
% stem(tt(1:20),mod_real_noise(1:20));
% title('carriers inphase');
% xlabel('Time(s)');
% ylabel('Amplitude');
% subplot(212);
% stem(tt(1:20),mod_imag_noise(1:20));
% title('carriers Quadrature');
% xlabel('Time(s)');
% ylabel('Amplitude');
     
%        f=(2/T)*(1:(FS))/(FS);
%   subplot(211);
%   plot(f,abs(fft(modulated ,FS))/FS);
%   title('Spectrum at Channel');
%   xlabel('Frequency(Hz)');
%   ylabel('Amplitude');
%   subplot(212);
%   pwelch(modulated,[],[],[],2/T);
    
%    [Pxx,W] = pwelch(modulated,[],[],4096,20);    
%    plot([-2048:2047]*fc/2048,10*log10(fftshift(Pxx)));
%    xlabel('frequency, MHz')
%    ylabel('power spectral density')
%    title('Transmit spectrum OFDM ');
   
% .................................DeMODULATION..........................
   
    for i=1:1:20480
      demod_real(i) = 2*(mod_real_noise(i)*cos(2*pi*fc*t(i)));
      demod_imag(i) = 2*(mod_imag_noise(i)*sin(2*pi*fc*t(i)));
    end
    
     [b,a] = butter(5,0.6);  
     real_fiter = filter(b,a,demod_real);
     imag_fiter = filter(b,a,demod_imag);
     
     Dsc_real=downsample(real_fiter,10);
     Dsc_imag=downsample(imag_fiter,10);
     Dsc = Dsc_real + Dsc_imag;
     
%      f=(2/T)*(1:(FS))/(FS);
%     subplot(211);
%     plot(f,abs(fft(Dsc ,FS))/FS);
%     title('filter Output');
%     xlabel('Frequency(Hz)');
%     ylabel('Amplitude');
%     subplot(212);
%     pwelch(Dsc,[],[],[],2/T);
     % .................................FFT..........................
     real_fft = fft(Dsc,2048);
   
     %.............................Removing zeros.....................    
    sct= zeros (1,1706);
     for i = 1:1:853
            sct(i)= real_fft(i);
     end
     
     k=854;
    for (i = 1196:1:2048)
         sct(k) = real_fft(i);
          k =k+1;
    end

%      scatterplot(sct);
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check1 = real(sct);
check2 = imag(sct);
   for i =1:1:1706;
          if(check1(i)> 2)
           real_check1(i)= 3;
          elseif(check1(i)< -2)
           real_check1(i)= -3;
          elseif(check1(i)> 0 && check1(i)< 2)
           real_check1(i)= 1;
          elseif(check1(i)< 0 && check1(i)> -2)
            real_check1(i)= -1;
          end
   end
   
      for i =1:1:1706;
          if(check2(i)> 2)
           imag_check2(i)= 3;
          elseif(check2(i)< -2)
           imag_check2(i)= -3;
          elseif(check2(i)> 0 && check2(i)< 2)
           imag_check2(i)= 1;
          elseif(check2(i)< 0 && check2(i)> -2)
            imag_check2(i)= -1;
          end
      end
   
      fft_data = real_check1 + imag_check2*1i;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      scatterplot(fft_data);

%.............................deMapping 16-QAM constellation.....................
%       demapper= zeros (1,1706);
for i = 1:1:1706
       if(fft_data(i) == -3+3i)
          demapper(i) = 0;
       elseif(fft_data(i) == -3+1i)
          demapper(i) = 1;
       elseif(fft_data(i) == -3-1i)
          demapper(i) = 2;
       elseif(fft_data(i) == -3-3i)
          demapper(i) = 3;
       elseif(fft_data(i) == -1+3i)
          demapper(i) = 4;
       elseif(fft_data(i) == -1+1i)
          demapper(i) = 5;
       elseif(fft_data(i) == -1-1i)
          demapper(i) = 6;
        elseif(fft_data(i) == -1-3i)
          demapper(i) = 7;
        elseif(fft_data(i) == 1+3i)
          demapper(i) = 8;
        elseif(fft_data(i) == 1+1i)
          demapper(i) = 9;
        elseif(fft_data(i) == 1-1i)
          demapper(i) = 10;
        elseif(fft_data(i) == 1-3i)
          demapper(i) = 11;
        elseif(fft_data(i) == 3+3i)
           demapper(i) = 12;
         elseif(fft_data(i)== 3+1i)
           demapper(i) = 13;
         elseif(fft_data(i)== 3-1i)
           demapper(i) = 14;
          elseif(fft_data(i)== 3-3i)
           demapper(i) = 15;
      end           
end
%.............................decimal 2 binary.....................
   hold2=de2bi(demapper,'left-msb');    
   k=1;
   for j = 1:1:1706;
     for i = 1:1:4;
       binary_out2(1,k) = hold2(j,i);
       k = k+1;
     end
  end
   
   n=15;
   binary_out= zeros (1,6810);
for i = 1:1:6810 
    binary_out(i)= binary_out2(n);
     n=n+1;  
end

[no_of_error(kk),ber(kk)]=biterr(binary_data , binary_out);
end

%.............................BER..............................
% figure; 
% stem(binary_out(1:50),'filled');
% title('Uncoded Received bits');
% xlabel('Binary index');
% ylabel('Binary value');
% figure; 
% stem(binary_data(1:50),'filled');
% title('Generated binary');
% xlabel('Binary index');
% ylabel('Binary value');


% berD = (1/4)*3/2*erfc(sqrt(4*0.1*(10.^(EbNo/10))));
% figure ; semilogy(EbNo,ber,'--b','linewidth',2);
% hold on;
% semilogy(EbNo,berD,'--*r','linewidth',2);
% xlabel('Eb/No (dB)');
% ylabel('Bit Error Rate')
% legend('Uncoded','Theoritical')
% grid on ;
% hold off;