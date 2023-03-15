% Kane O'BRIEN
% A1880046
% Signal Processing Applications 2023; Brian NG.
% Computer Exercise 1
% 7/3/23
clc; clear all; close all;
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Question 1:
clc; clear all; close all; fig_n=1;
% Consider a matrix of 4 random variables at three timesteps; determine the
% autocovariance manually; compare with the covariance function.
% 1.1
fprintf("The considered matrix is:")
x = [11, -2, 3; -4, -7, 13; 15, -6, 4; -1, 8, 7]

x_bar = mean(x); % columnwise by default
x_hat = x-x_bar;

% 1.2, 1.3
[N] = size(x);
xvar = zeros(N(2));
N=N(1);
for i = 1:N
    xvar = xvar + x_hat(i,:)' * x_hat(i,:);
end

% 1.4
fprintf("Via autocorrelation from first principles:\n")
x_autocov = xvar/(N-1)

% 1.5
fprintf("Via the Covariance function:\n")
x_cov = cov(x)

if x_autocov ~= x_cov
    fprintf("The matrix dont match; Houston we have a problem\n")
elseif ~ ishermitian(x_autocov)
    fprintf("The matrix is not Hermitian; there is a problem\n")
else
    fprintf("The matricies match\n")
end


% 1.6
% Generate a 100x10 matrix of random variables; determine the covariance.
% is this whats expected? How can you visualise it?
y = randn(100,10);

% 1.7
fprintf("The covariance function estimation is:\n")
y_cov = cov(y)

% First, it is the expectation that the resulting covariance matrix would 
% be of shape 10x10 - as for an MxN matrix, the covariance is NxN.
% Second; due to the diagonalisation of the matrix (in transpose-
% multiplication) we expect that the auto-covariance matrix to be Hermitian
% That said, it is typical for covariance to be real-only hence
% only real-symmetric.

if ~ishermitian(y_cov)
    fprintf("The matrix is not Hermitian/ real-symmetric; if you used the cov function you're in trouble.\n")
else
    fprintf("The matrix is Hermitian/ Real-Symmetric.\n")
end

% 1.8
% Visualising covariance is probably best done with a heatmap/ surf map.
% When the x(N) correlates with y(N) it is positive covariance; a large
% covariance leads to a clear diagonal heatmap;

figure(fig_n)
    fig_n = fig_n+1;
    fig = heatmap(y_cov);
    colormap(parula)
    fig.YDisplayData = flipud(fig.YDisplayData); 
    title("Visualisation of Autocovariance")
    xlabel("X[n]")
    ylabel("X[n]")
% Of course it will be a different distribution every time, but off the
% diagonal, the variance never seems to go above 0.2~0.3 whilst on the
% diagonal, the variance is typically above 0.8. Thresholding/ Conditional
% formatting is another way this could be visualised.

%% Question 2:
clc; clear all;
% A fundamental limitation of DFT-based spectrum estimation techniques is 
% imposed by the number of available samples. Construct N samples of the 
% signal, where N = (20,100,1000);

N = [20, 100, 1000];
F = [0.05, 0.1, 0.102];

syms x(n)
x(n) = 1*cos(2*pi*F(1)*n) + 1*cos(2*pi*F(2)*n) + 1*cos(2*pi*F(3)*n);

x1 = double( x(0:N(1)-1));
x2 = double( x(0:N(2)-1));
x3 = double( x(0:N(3)-1));

X1 = abs(fft(x1)).^2;
X2 = abs(fft(x2)).^2;
X3 = abs(fft(x3)).^2;

figure(fig_n)
  fig_n = fig_n+1;
  sgtitle("DFT$(x(n) = 1\cdot\cos(2\pi0.05n) + 1\cdot\cos(2\pi0.1n) + " + ...
      "1\cdot\cos(2\pi0.102n))$")
  subplot(3,3,[1,2,3])
    plot(0:N(1)-1, X1)
    str= "N = " + num2str(N(1)); 
    title(str)
    xlim([0 (N(1)-1)])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,3,[4,5,6])
    plot(0:N(2)-1, X2)
    str= "N = " + num2str(N(2));
    title(str)
    xlim([0 N(2)-1])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,4,[9,10])
    plot(0:N(3)-1, X3)
    str= "N = " + num2str(N(3));
    title(str)
    xlim([0 (N(3)-1)])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,4,[11,12])
    plot(0:N(3)-1, X3)
    title("Zoomed "+str)
    xlabel("Freq Bin (n)")
    ylabel("Mag")
    xlim([48 104])
    xticks(48:2:104)


% Due to use of Cos, the resulting signal is real and hence two symmetric
% bins are produced in the DFT. 

% Within the first subplot where N=20, the frequency resolution is very
% poor. Whilst there is three superimposed signals - their frequency
% spacing given by the equation (f0=0.05rad/n, f1= 0.1rad/n, f2=0.102rad/n)
% is too fine that they cannot be resolved individually. The first plot
% shows two of the signals (f1,f2) are entirely contained within the same
% frequency bin #2,18 and f0's energy can be seen in the adjacent bins
% #1,19 however cannot be truely resolved.

% Within the second subplot, N=100, the frequency resolution is improved
% such that f0 can be resolved from f1,f2, found in frequency bins #5,95
% and #10,90 respectively. f1 and f2 cannot be resolved individually.

% Finally within the last subplots, N=1000, the frequency resolution is
% improved again such that f1 and f2 can now be resolved. f0 can be seen in
% bins #50,950, f1 is seen in #100,900 and f2 in #102,898

% This shows that the resolution rule \Delta f = 1/N, where \Delta f is
% integer. 

% Verifying the first frequency bins (Symmetric)
F1a = floor(F(1) * N(1));
F2a = floor(F(2) * N(1));
F3a = floor(F(3) * N(1));

F1b = floor(F(1) * N(2));
F2b = floor(F(2) * N(2));
F3b = floor(F(3) * N(2));

F1c = floor(F(1) * N(3));
F2c = floor(F(2) * N(3));
F3c = floor(F(3) * N(3));

F_bins = [F1a, F1b, F1c; F2a, F2b, F2c; F3a, F3b, F3c]
% We cannot resolve frequencies when they exist in the same bin (or even in
% the adjacent bin) that  F1, F2, F3 show.


% 2.1
% Repeat the above using Hamming windowed DFT
w1 = hamming(N(1));
w2 = hamming(N(2));
w3 = hamming(N(3));

XW1 = abs(fft(x1.*w1')).^2;
XW2 = abs(fft(x2.*w2')).^2;
XW3 = abs(fft(x3.*w3')).^2;

figure(fig_n)
  fig_n = fig_n+1;
    sgtitle("DFT$(x(n)\times w(n))$ - Hamming window DFT")
  subplot(3,3,[1,2,3])
    plot(0:N(1)-1, XW1)
    str= "N = " + num2str(N(1)); 
    title(str)
    xlim([0 (N(1)-1)])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,3,[4,5,6])
    plot(0:N(2)-1, XW2)
    str= "N = " + num2str(N(2));
    title(str)
    xlim([0 N(2)-1])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,4,[9,10])
    plot(0:N(3)-1, XW3)
    str= "N = " + num2str(N(3));
    title(str)
    xlim([0 (N(3)-1)])
    xlabel("Freq Bin (n)")
    ylabel("Mag")
  subplot(3,4,[11,12])
    plot(0:N(3)-1, XW3)
    title("Zoomed "+str)
    xlabel("Freq Bin (n)")
    ylabel("Mag")
    xlim([48 104])
    xticks(48:2:104)

figure(fig_n)
  fig_n = fig_n+1;
    plot(0:N(3)-1, X3)
    hold on
    plot(0:N(3)-1, XW3)
    title("Zoomed "+str)
    xlabel("Freq Bin (n)")
    ylabel("Mag")
    xlim([48 104])
    ylim([0 280000])
    xticks(48:2:104)
    title("Unwindowed DFT vs. Hamming Window DFT")
    legend("Un-Windowed","Hamming Window")

% When we use the windowing function on real sinewaves ona N-point DFT it 
% is of limited value as the energy is lost to the window, and spectral 
% spreading occurs. N-point DFT does not have ~ periodocity issues with 
% continuity errors as the FFT has. Windows become very useful for the 
% signals that are non-periodic.

%% Question 3:
clc; clear all;
% Compare the following spectral estimation methods covered in the 
% prescribed reading:
% - Welch (w/ Parameters)
% - Blackman-Tukey
% - AR (Burg method)
% - MUSIC
% Upon a N-Samples triple tone in complex,white noise

sigma = 0.1;
F = [0.3 0.8 1.9]/2;
A = [1 1 1];

for N = [20,200]
  [X_welch,btwindowlen,X_BT,ar_order,X_AR,music_order,X_music] = specdestimate(N,sigma,A,F);
  figure(fig_n)
    fig_n = fig_n+1;
      subplot(4,1,1)
        plot(0:(2)/length(X_welch):2-((2)/length(X_welch)),X_welch)
        title("Welch's Method")
      subplot(4,1,2)
        plot(0:(2)/length(X_BT):2-((2)/length(X_BT)),abs(X_BT))
        title("Blackman-Tukey : m = "+ btwindowlen)
      subplot(4,1,3)
        plot(0:(2)/length(X_AR):2-((2)/length(X_AR)),abs(X_AR))
        title("AutoRegression (Burg Method): order= "+ ar_order)
      subplot(4,1,4)
        plot(0:(2)/length(X_music):2-((2)/length(X_music)),abs(X_music))
        title("MUSIC: order= "+ music_order)
      sgtitle("Power Spectral Density Estimates: N= "+ num2str(N))
        for ii = 1:4
            subplot(4,1,ii)
            xlim([0 2])
            xticks(0:0.1:2)
            grid on
            xlabel("Normalised Frequency ($\times\pi$ rad/Sample)")
            ylabel("$|Mag|$")
        end
end

% What we see is that DFT based spectrum estimation requires a high number
% of samples to estimate accurately. This is demonstrated with N=20 in
% figure 5, in which Welch and Blackman-Tukey do not find peaks, or 
% the frequency bin of the energy. Parametric based spectrum estimation
% requires some heuristic tuning (testing variable model orders) however we
% can see for the provided parameters, Burg-3 and MUSIC-8, the frequency
% bins are estimated accurately, however the amplitudes and hence energy is
% not.
% Figure 6 shows the same three signals in noise with N=200; the DFT-based
% estimators in this case function very well, and show a matched amplitude
% (as created the signals). In this case, the model based estimators still
% show peaks in the 0.3 and 0.8 frequency bins, however their peak is very
% weak. This again, likely falls into some heuristic tuning and changing of
% model order.

%% Q4
% pkg load signal # Octave funnyness
load('SPA2023-computer-exercises1-data/fm.mat');

t= length(x)/fs;
fprintf("The sample runs for %3.2f seconds\n",t)
n =12
N = 2^n;
window = N; Nfft = N
overlap = .5; % Default is 0.5, play with parameter sometimes.

X  = 10*log10(pwelch(x,window,overlap,Nfft).^2);
% Octave ignores un-full Nfft windows and doesn't zero pad default;
% Double check matlab behaviour

str= strcat("Power spectral density of the FM Broadcast Channel, N= 2^ ", num2str(n));
figure(fig_n)
  fig_n = fig_n+1;
  clf
  subplot(2,3,[1 2 3])
    plot(0 :1/N:(1)-(1/N), X)
    title(str)
    xlim([0 1])
    ylim([0 120])
    xticks([0:0.05:2])
    yticks([0:5:120])
    grid on
    ylabel("Log-Mag (dB)")
    xlabel("Normalised Frequency ($\times 2\pi$ rad/Sample)")
  subplot(2,3,[4])
    plot(0 :1/N:(1)-(1/N), X)
    title("Zoomed section of the FM Broadcast Channel")
    xlim([0.35 0.45])
    ylim([40 120])
    xticks([0.35:0.01:0.45])
    grid on
    ylabel("Log-Mag (dB)")
    xlabel("Normalised Frequency ($\times 2\pi$ rad/Sample)")
  subplot(2,3,[5 6])
    plot(0 :1/N:(1)-(1/N), X)
    title("Zoomed 'Symmetric' component of FM Channel")
    xlim([0.4 0.44])
    ylim([40 120])
    xticks([0.4:0.005:0.44])
    yticks([40:5:120])
    grid on
    ylabel("Log-Mag (dB)")
    xlabel("Normalised Frequency ($\times 2\pi$ rad/Sample)")
    hold on
    text(0.41,90, "19kHz Pilot")
    text(0.401,118, "Carrier @ 800kHz offset")
    text(0.428, 50, "OOB/ RDS Subchannels")


% FM Carriers are distinct due to their Bessel-curve-with-horns appearance;
% The main peak is carrier / mixing oscillator; the "Horns" are due to the 19kHz
% Pilot Tones
% Broadcast FM is powerful in urban environments,  these will be the largest
% signals present. In the given samples, there appears to be some LO spurs or IMD
% present on the spectrum.  Samples look to be 8-bit IQ, and best guess is
% capture from RTLSDR(RTL2832U) or derivative SDR.

% Commerical FM channels are known to use Wideband-Frequency Modulation. This is
% a typical 200 kHz wide channel to allow Pilot carriers (tones to allow cheap
% IF), SUM (Mono, L+R) and DIFFERERNCE (stereo, L-R) channels and out of band
% carriers (RDS, etc). Different modulation schemes are used for the OOB, notably
% the DSBSC for the difference channel (38kHz) and BPSK for the RDS (57kHz)
% See [FM Broadcast Radio](https://www.sigidwiki.com/wiki/FM_Broadcast_Radio)
% for more detail.

% I suspect the spectrum may look more familiar if we use the fft shift. This
% figure now shows the "curtains" / digital filter roll off from the receiver
% itself; further, its easier to see the noise floor and spurs from this view
figure(fig_n)
  fig_n = fig_n+1;
  clf
    plot(0 :1/N:(1)-(1/N), fftshift(X))
    title(str)
    xlim([0 1])
    ylim([0 120])
    xticks([0:0.05:2])
    yticks([0:5:120])
    hold on
    grid on
    ylabel("Log-Mag (dB)")
    xlabel("Normalised Frequency ($\times 2\pi$ rad/Sample)")
    plot([0 1],[47.5,47.5],'r')
    text(0.5,108, "Local Oscillator")
    text(0.8,108, "WFM Broadcast")
    text(0.3,45,"Noise Floor",'color','red')
    text(0.13,65, "IMD/ Spurs")
    text(0.06,25, "Filter roll-off 'Curtains'")
    
    
%% Q5

%% Q6

%% Functions

function [X_welch,btwindowlen,X_BT,ar_order,X_AR,music_order,X_music] = specdestimate(N,sigma,Amps,Fns)
  % Complex noise; Gaussian Distributed
    w  = (sigma/sqrt(2)) * (randn(1,N)+1j*randn(1,N));
    
  % Symbolic, the lazy mans way.
    syms x(a,fn,n)
    x(a,fn,n) = exp(1j*2*pi*fn*n);
  % Create three signals
    x1 = double( x(Amps(1), Fns(1), 0:(N-1)));
    x2 = double( x(Amps(2), Fns(2), 0:(N-1)));
    x3 = double( x(Amps(3), Fns(3), 0:(N-1)));
  % Combine signals and noise
    x = x1+x2+x3+w;
    
  % 3.1 - Welch Method
    X_welch = pwelch(x);
    
  % 3.2 - Blackman-Tukey Method
    %   Windowed autocorrelation with Bartlett window (Ham,Hann do not satisfy BT)
    btwindowlen = N/4; % Tunable
    w_bt = bartlett(btwindowlen);
    x_cov = xcorr(x);
    w_bt = [w_bt', zeros([1,length(x_cov)-length(w_bt)])]; % Zero packed (matrix length match)
    X_BT = abs(fft(x_cov.*w_bt));
    
  % 3.3 - Burg AR Method
    ar_order = 3; %Tuneable
    X_AR = pburg(x,ar_order);
    
  % 3.4 - MUlitple SIgnal Classification
    music_order = 8;  % Tuneable
    X_music = pmusic(x, music_order);
end
