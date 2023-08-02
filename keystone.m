% keystone
%
% Demo of keystone formatting for correcting range-Doppler
% measurements for range migration.
%
% This code closely follows the equations in the tech memo "The Keystone
% Transformation for Correcting Range Migration in Range-Doppler
% Processing" by Mark A. Richards, Mar. 2014, available at www.radarsp.com.
%
% Mark Richards
%
% March 2014
clear all
close all
c = 3e8; % speed of light
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 128; % fast time dimension, samples
M = 101; % slow-time dimension, samples; keep it odd.
% K_L and K_M must be even to avoid labeling problems later
K_L = 2^(nextpow2(128)+1); % fast-time DFT size for interpolation and shifting
K_M = 2^(nextpow2(512)+1); % slow-time DFT size
% Lref is the range bin # of the target, on a 0:L-1 scale, at the center of the CPI (middle pulse)
% Lref = round(L/2); % puts target at middle range bin on the middle pulse
Lref = 100;
F0 = 1000e6; % RF (Hz)
B = 200e6; % 波形带宽 (Hz)
v = 440; % 朝向雷达的速度 velocity in m/s 
% 采样间隔和速率
Fsft = 2.3*B;
PRF = 10e3;
% sinc滤波器阶数
Nsinc = 11;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END USER INPUT SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some derived parameters
m_end = (M-1)/2;
ms = (-m_end:m_end); % slow time index labels
Fd = 2*v*F0/c; % 多普勒频移, Hz
Tft = 1/Fsft; % 快时间采样间隔
dr = c*Tft/2; % range bin 距离门
Tst = 1/PRF; % 慢时间采样间隔

Dfd = 1/M; % Rayleigh Doppler resolution in cycles/sample
Drb = (1/B)/Tft; % Rayleigh fast-time resolution in range bins 快速采样时间段分辨率
if (PRF < Fd/2)
    fprintf('\nWarning: PRF < Fd/2. PRF = %g, Fd = %g.\n',PRF,Fd)
end
% 计算并报告在距离范围内所产生的距离门徙动 ？？
RM = v*Tst/dr; % amount of range migration per pulse in range bins
RMtot = M*RM; % total range migration over the dwell, in range bins
fprintf('\nTotal range migration = %g m = %g range bins.\n',RMtot*dr,RMtot)
if ( (floor(Lref-RMtot/2) < 0) || (ceil(Lref+RMtot/2) > L-1) )
    fprintf(['\nWarning: Target will migrate out of range.', ...
    'RMtot = %g range bins, Lref = %g, L = %g range bins.\n'],RMtot,Lref,L)
end
% 计算归一化多普勒并wrap it
fd = Fd*Tst;
fdn = mod(fd + 0.5,1) - 0.5; % alias it back into [-0.5,+0.5]
amb_num = round(fd - fdn); % number of Doppler wraps
% Define corners of a box centered on the expected target coordinates, and
% one Rayleigh width wide in each direction and each dimension (i.e., a
% null-to-null resolution box for a well-formed sinc spectrum). Will use
% this to draw such a resolution box on some of the figures.
L1 = Lref - Drb;
L2 = Lref + Drb;
fd1 = fdn - Dfd;
fd2 = fdn + Dfd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create synthetic data. First compute pulse-to-pulse phase shift. Then set
% up a matrix. Loop over pulses, computing current range. The range
% profile of the compressed data for a single pulse is assumed to be a sinc
% function with a zero spacing equal to 1/B seconds, and a phase shift of
% the usual -(4*pi*F0/c)*R, where R = Rref - v*Tst*m and m = pulse number.
% Don't worry about amplitude. Also don't bother with the
% -(4*pi*F0/c)*Rref phase term, it is the same for all pulses.
del_phi = -4*pi*(F0/c)*v*Tst; % pulse-to-pulse phase increment due to range change
y = zeros(L,M);
% 慢时间做频率压缩
for m = 1:M
    mm = ms(m); % counts from -(M-1/2) to +(M-1)/2
    y(:,m) = exp(-1i*del_phi*mm)*sinc( B*Tft*((0:L-1)-Lref+v*Tst*mm/dr) );
end
% Now examine the data. Then compute the range-Doppler matrix with a
% slow-time DFT and look at that.
figure
imagesc(ms,0:L-1,abs(y))
grid
colormap(flipud(gray))
xlabel('pulse number')
ylabel('range bin')
title('Raw Fast-time/Slow-time Data Pattern')
% 存在疑问
Y_rD = fft(y,K_M,2);
Y_rD = fftshift(Y_rD,2);
Y_rD_dB = db(abs(Y_rD),'voltage');


Y_rD_dB = Y_rD_dB - max(Y_rD_dB(:)); % normalize to 0 dB max
Y_rD_dB(:) = max(-40,Y_rD_dB(:)); % limit dynamic range for plot purposes
fD = (-K_M/2:K_M/2-1)/K_M; % this only works correctly if K_M is even

% 原来存在注释的
figure
mesh(fD,0:L-1,Y_rD_dB)
xlabel('normalized Doppler')
ylabel('range bin')
title('Raw Range-Doppler Matrix')

% hvline--------注释
figure
imagesc(fD,0:L-1,Y_rD_dB)
% % hline(Lref,':w'); vline(fdn,':w') % mark the correct spectrum center
% line([fd1 fd1 fd2 fd2 fd1],[L1 L2 L2 L1 L1],'Color','w','LineWidth',2) % resolution box
title(['Rng Migration = ',num2str(RMtot),...
 ' bins over CPI; Normalized Doppler = ',num2str(fdn),' cyc/samp'])
xlabel('normalized Doppler')
ylabel('range bin')
colorbar
shg

% ------- 另外多的两幅图
% It is convenient to look at the fast time DFT of the raw data as well. We
% will need this product as the starting point for keystoning a little
% further down. Apply a fast-time DFT. fftshift in fast-time freq dimension
% to center the origin. This will be the frequency corresponding to F0.
% Also work out axis label in Hz.
% Y_Rd是什么
Y_Rd = fftshift(fft(y,K_L,1),1);
Fl = (-K_L/2:K_L/2-1)/K_L*Fsft;
figure
subplot(1,2,1)
imagesc(ms,F0+Fl,abs(Y_Rd))
xlabel('slow time')
ylabel('fast time frequency (Hz)')
title('Magnitude after fast-time FT')
colorbar
subplot(1,2,2)
imagesc(ms,F0+Fl,angle(Y_Rd))
xlabel('slow time')
ylabel('fast time frequency (Hz)')
title('Unwrapped phase after fast-time FT')
colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now do compensation by explicit shifting and interpolation in range for
% each pulse to see that that does indeed work (for a single target with
% known velocity) and verify formulas
%现在在范围内通过显式移位和插值进行补偿%的脉冲，看看是否确实有效(对于单个目标已知速度)并验证公式
% 显示移位操作，在慢时间维度上
yc = zeros(size(y)); % this will be the compensated fast-time/slow-time data
wl = 2*pi*(-K_L/2:K_L/2-1)'/K_L; % normalized fast time freq in radians/sample
for m = 1:M % loop over pulses
    mm = ms(m); % counts from -(M-1/2) to +(M-1)/2
    Lm = v*Tst*mm/dr; % # of range bins of shift needed
    Y_shift = fftshift(fft(y(:,m),K_L)).*exp(-1i*wl*Lm);
    y_shift_temp = ifft(ifftshift(Y_shift));
    y_shift(:,m) = y_shift_temp(1:L);
end % of loop over pulses

figure
imagesc(ms,0:L-1,abs(y_shift))
colormap(flipud(gray))
grid
xlabel('pulse number')
ylabel('range bin')
title('Data Pattern after Compensation by Shifting')
Y_shift = fft(y_shift,K_M,2);
Y_shift = fftshift(Y_shift,2);
Y_shift_dB = db(abs(Y_shift),'voltage');
Y_shift_dB = Y_shift_dB - max(Y_shift_dB(:)); % normalize to 0 dB max
Y_shift_dB(:) = max(-40,Y_shift_dB(:)); % limit to 40 dB range
figure
mesh(fD,0:L-1,Y_shift_dB)
xlabel('normalized Doppler')
ylabel('range bin')
title('Range-Doppler Matrix after Compensation by Shifting')

% hvline--------注释
% figure
% imagesc(fD,0:L-1,Y_shift_dB)
% hline(Lref,':w'); vline(fdn,':w') % mark the correct spectrum center
% line([fd1 fd1 fd2 fd2 fd1],[L1 L2 L2 L1 L1],'Color','w','LineWidth',2) % resolution  box
% xlabel('normalized Doppler')
% ylabel('range bin')
% title('Range-Doppler Matrix after Compensation by Shifting')
% shg
% colorbar
% Let's also look at the fast-time DFT of this data. It will show what

% 移位--------画彩图
% we're trying to get to with the keystoning in the next section.
% Y_Rd_shift = fftshift(fft(y_shift,K_L,1),1);
% figure
% subplot(1,2,1)
% imagesc(ms,F0+Fl,abs(Y_Rd_shift))
% xlabel('slow time')
% ylabel('fast time frequency (Hz)')
% title('Magnitude after fast-time FT')
% colorbar
% subplot(1,2,2)
% imagesc(ms,F0+Fl,angle(Y_Rd_shift))
% xlabel('slow time')
% ylabel('fast time frequency (Hz)')
% title('Unwrapped phase after fast-time FT')
% colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 从Y_Rd开始，即在快时间内而不是在慢时间内进行数据DFT。
% 对于每个快时频率仓，计算一个新的插值慢时序列。
% 使用现有的函数sinc_interp为带宽有限，汉明加权插值做的工作
% Now start over and correct by keystoning. Begin with Y_Rd, i.e. data
% DFT'ed in fast time but not in slow time. For each fast-time frequency
% bin, compute a new, interpolated slow-time sequence. Use the existing
% sinc_interp function for bandlimited, Hamming-weighted interpolation to
% do the work.
Y_Rd_key = zeros(size(Y_Rd));

% 在快时间内做数据插值
for k = 1:K_L
    [y_temp,mm_i] = sinc_interp(Y_Rd(k,:),ms,(F0/(F0+Fl(k)))*ms,Nsinc,1);
 % y_temp = interp1(ms,Y_Rd(k,:),(F0/(F0+Fl(k)))*ms,'spline');
 % Mi will always be odd the way I'm setting up the problem. Also, Mi <= M.
    Mi = length(y_temp);
    dM = M - Mi; % dM will be even so long as M and Mi are odd
    Y_Rd_key(k,1+dM/2:1+dM/2+Mi-1) = y_temp; % center the interpolated data in slow time
end
% Now correct the modified spectrum for the ambiguity number of the
% Doppler. This code uses the ambiguity number of the first target. So if
% the other targets have a different ambiguity number it won't be correct
% for them. The first version of the correction corresponds to the Li et al
% paper and is consistent with my memo. The second (commented out)
% corresponds to the Perry et al paper and can be obtained from the first
% using a binomial expansion approximation to (F0/(F0+Fl)). Either one
% works if done either after the keystone correction, as is done here, or
% before.
for mp = 1:M
    for k = 1:K_L
        mmp = ms(mp); % counts from -(M-1/2) to +(M-1)/2
        Y_Rd_key(k,mp) = Y_Rd_key(k,mp)*exp(1i*2*pi*amb_num(1)*mmp*(F0/(F0+Fl(k))));
% Y_Rd_key(k,mp) = Y_Rd_key(k,mp)*exp(-1i*2*pi*amb_num(1)*mmp*(Fl(k)/F0));
    end
end
% Now IDFT in fast-time and DFT in slow time to get range-Doppler matrix
y_temp_key = ifft( ifftshift(Y_Rd_key,1),K_L,1 );
y_rd_key = y_temp_key(1:L,:);
Y_rD_key = fftshift( fft(y_rd_key,K_M,2),2 );
Y_rD_key_dB = db(abs(Y_rD_key),'voltage');
Y_rD_key_dB = Y_rD_key_dB - max(Y_rD_key_dB(:)); % normalize to 0 dB max
Y_rD_key_dB(:) = max(-40,Y_rD_key_dB(:)); % limit to 40 dB range
figure
imagesc(ms,0:L-1,abs(y_rd_key))
grid
colormap(flipud(gray))
xlabel('slow-time')
ylabel('fast time')
title('Keystoned Fast-time/Slow-time Data Pattern')
figure
mesh(fD,0:L-1,Y_rD_key_dB)
xlabel('normalized Doppler')
ylabel('fast time')
title('Keystoned Range-Doppler Matrix')

% hvline--------注释
% figure
% imagesc(fD,0:L-1,Y_rD_key_dB)
% hline(Lref,':w'); vline(fdn,':w') % mark the correct spectrum center
% line([fd1 fd1 fd2 fd2 fd1],[L1 L2 L2 L1 L1],'Color','w','LineWidth',2) % resolution box
% xlabel('normalized Doppler')
% ylabel('range bin')
% title('Keystoned Range-Doppler Matrix')
% shg
% colorbar

figure
 subplot(1,2,1)
imagesc(ms,F0+Fl,abs(Y_Rd_key))
xlabel('slow time')
ylabel('fast time frequency (Hz)')
title('Magnitude after interpolation')
colorbar
subplot(1,2,2)
imagesc(ms,F0+Fl,angle(Y_Rd_key))
xlabel('slow time')
ylabel('fast time frequency (Hz)')
title('Unwrapped phase after interpolation')
colorbar
