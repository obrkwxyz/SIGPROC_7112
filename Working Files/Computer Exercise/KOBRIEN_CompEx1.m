% Kane O'BRIEN
% A1880046
% Signal Processing Applications 2023; Brian NG.
% Computer Exercise 1
% 7/3/23

set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');

%% Question 1:
clc; clear all; close all
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

figure(1)
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
clc; clear all; close all
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

figure(2)
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


% Through use of Cos, the resulting signal is real and hence two symmetric
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
% the adjacent bin) that  F1,F2,F3 show.


% 2.1
% Repeat the above using Hamming windowed DFT
w1 = hamming(N(1));
w2 = hamming(N(2));
w3 = hamming(N(3));

XW1 = abs(fft(x1.*w1')).^2;
XW2 = abs(fft(x2.*w2')).^2;
XW3 = abs(fft(x3.*w3')).^2;

figure(3)
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

figure(4)
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
clc; clear all; close all
% Compare the following spectral estimation methods covered in the 
% prescribed reading:
% - Welch (Default Parameters)
% - Blackman-Tukey
% - AR (Burg method)
% - MUSIC
% Upon a N-Samples triple tone in white noise

N = [20, 200];
f = [0.04 -0.05 0.5];
A = [1 1 1];
sigma = 0.1;

syms x(a,fn,n)
x(a,fn,n) = exp(1j*2*pi*fn*n);

% Random noise, variance of sigma
w  = (sigma/sqrt(2)) * (randn(1,N(1))+1j*randn(1,N(1)));

% Create the signals
x1 = double( x(A(1), f(1), 0:(N(1)-1)));
x2 = double( x(A(2), f(2), 0:(N(1)-1)));
x3 = double( x(A(3), f(3), 0:(N(1)-1)));

% Combine Signals
x = x1+x2+x3+w;

% Use builtin pwelch()
X_welch = pwelch(x);

% Perform Blackman-Tukey Method
%   Windowed autocorrelation with Bartlett window
%   Hamming / Hanning windows do not satisfy BT method.
window_length = N(1)/4; % tunable parameter

x_cov = xcorr(x);
w_bt = bartlett(window_length);
w_bt = [w_bt', zeros([1,length(x_cov)-length(w_bt)])];
X_BT = abs(fft(x_cov.*w_bt));

% Use builtin arburg()
ar_order = 2;
X_AR = pburg(x,ar_order);


% Via MUSIC
pmusic_order =8;
X_music = pmusic(x, pmusic_order);


figure(5)
  subplot(4,1,1)
    plot(X_welch)
    xlim([0 length(X_welch)])
    title("Power Spectral Estimate via pwelch()")
  subplot(4,1,2)
    plot(real(X_BT))
    xlim([1 length(X_BT)])
    title("Power Spectral Estimate via Blackman-Tukey, m= "+window_length)
  subplot(4,1,3)
    plot(real(X_AR))
    xlim([0 length(X_AR)])
    title("Power Spectral Estimate via AutoRegression (Burg Method) - pburg(), order= "+ar_order)
  subplot(4,1,4)
    plot(real(X_music))
    xlim([0 length(X_music)])
    title("Power Spectral Estimate via MUSIC - pmusic(), order= "+ pmusic_order)
 


%     hold on
%     plot(imag(X_BT)


    