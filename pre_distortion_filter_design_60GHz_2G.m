%% pre-distortion filter method 1

%% run Cal5GHz_r1_wb_jk.
% BB_BW (baseband bandwidth) in line 164 is now 6(GHz).
% Select an element.
% output: 'jk_mean_ir' of size [81880x1] and of 40GHz BW
% 'jk_mean_ir' corresponds to the uncalibrated impulse response.
Cal5GHz_r1_wb_jk_60GHz_2Gbps;   % select an element, check BB_BW (2GHz? 6GHz?)
% npn: 81880
% spb: 40
% 81880 => 40GHzdfadsfgasgdsag
mean_ir = reshape(jk_mean_ir,npn,1);

mean_ir = reshape(jk_mean_ir,npn,1); %%% lkjlkj
%% resample 'jk_mean_ir'
% 40GHz to 12GHz
% interpolate jk_mean_ir with the interpolation factor of 12.
% decimate interp(jk_mean_ir,40) with the decimation factor of 40.
% output: 'jk_mean_ir_12GHz'
interp_factor = 12;
decim_factor = 40;
interp_mean_ir = interp([0;mean_ir(1:end-1)],interp_factor);
decim_mean_ir = decimate(interp_mean_ir,decim_factor);asdfadf

% Debug
if (1)
    % time domain
    org_ts = 1/spb;
    org_t_range = (1:npn)*org_ts;
    interp_ts = 1/spb/interp_factor;
    interp_t_range = (1:numel(interp_mean_ir))*interp_ts;
    decim_ts = 1/spb/interp_factor*decim_factor;
    decim_t_range = (1:numel(decim_mean_ir))*decim_ts;

    [~,max_index] = max(abs(mean_ir));
    max_time = org_t_range(max_index);
    xlim_min = floor(max_time)-10;
    xlim_max = floor(max_time)+50;
    
    figure(119);
    subplot(311);
    plot(org_t_range,10*log10(abs(real(mean_ir))),'bx:');
    hold on;
    plot(interp_t_range,10*log10(abs(real(interp_mean_ir))),'r-');
    plot(decim_t_range,10*log10(abs(real(decim_mean_ir))),'k.:');
    hold off;
    grid on;
    title('real part');
    legend('org','interp','decim','Location','best');    
    xlim([xlim_min xlim_max]);
    xlabel('time(ns)');
    ylabel('real part value(dB)');
    
    subplot(312);
    plot(org_t_range,10*log10(abs(imag(mean_ir))),'bx:');
    hold on;
    plot(interp_t_range,10*log10(abs(imag(interp_mean_ir))),'r-');
    plot(decim_t_range,10*log10(abs(imag(decim_mean_ir))),'k.:');
    hold off;
    grid on;
    title('imaginary part');
    legend('org','interp','decim','Location','best');
    xlim([xlim_min xlim_max]);
    xlabel('time(ns)');
    ylabel('imag part value(dB)');
    
    subplot(313);
    plot(org_t_range,10*log10(abs(mean_ir)),'bx:');
    hold on;
    plot(interp_t_range,10*log10(abs(interp_mean_ir)),'r-');
    plot(decim_t_range,10*log10(abs(decim_mean_ir)),'k.:');
    hold off;
    grid on;
    title('abs');
    legend('org','interp','decim','Location','best');
    xlim([xlim_min xlim_max]);
    xlabel('time(ns)');
    ylabel('magnitude(dB)');

    % freq domain
    org_f_range = [(-npn/2):(npn/2-1)]/npn*spb;
    org_fft = fftshift(fft(mean_ir));

    ninterp = npn*interp_factor;
    interp_f_range = [(-ninterp/2):(ninterp/2-1)]/ninterp*spb*interp_factor;
    interp_fft = fftshift(fft(interp_mean_ir));

    ndecim = ninterp/decim_factor;
    decim_f_range = [(-ndecim/2):(ndecim/2-1)]/ndecim*spb*interp_factor/decim_factor;
    decim_fft = fftshift(fft(decim_mean_ir));

    figure(120);
    subplot(211);
    plot(org_f_range,10*log10(abs(org_fft)),'bx:');
    hold on;
    plot(interp_f_range,10*log10(abs(interp_fft)),'r-');
    plot(decim_f_range,10*log10(abs(decim_fft)),'k.:');
    hold off;
    grid on;
    title('|fft|');
    legend('org','interp','decim','Location','best');
    xlim([-4 4]);
    ylim([-100 -30]);
    xlabel('freq(GHz)');
    ylabel('|fft|(dB)');
    
    subplot(212);
    plot(org_f_range,angle(org_fft),'bx:');
    hold on;
    plot(interp_f_range,angle(interp_fft),'r-');
    plot(decim_f_range,angle(decim_fft),'k.:');
    hold off;
    grid on;
    title('angle(fft)');
    legend('org','interp','decim','Location','best');
    xlim([-4 4]);
    xlabel('freq(GHz)');
    ylabel('phase');
end

jk_mean_ir_12GHz = decim_mean_ir;   % [2047x12 1]
jk_mean_ir_12GHz = reshape(jk_mean_ir_12GHz,numel(jk_mean_ir_12GHz),1);
[~,jk_mean_ir_12GHz_max_idx]=max(abs(jk_mean_ir_12GHz));
%{  
% fft-based sample shift instead of the simple circshift
jk_mean_ir_12GHz_fft = fftshift(fft(jk_mean_ir_12GHz));
jk_nn = numel(jk_mean_ir_12GHz_fft);
jk_ir_freq_vec = ((-jk_nn/2):(jk_nn/2-1))'/jk_nn;
jk_mean_ir_12GHz_fft_based_shift = jk_mean_ir_12GHz_fft...
    .*exp(-1j*2*pi*jk_ir_freq_vec*(-jk_mean_ir_12GHz_max_idx+1));
jk_mean_ir_12GHz_time_after_fft_based_shift = ifft(ifftshift(jk_mean_ir_12GHz_fft_based_shift));
%}
jk_mean_ir_12GHz = circshift(jk_mean_ir_12GHz,-jk_mean_ir_12GHz_max_idx+1);

%% run awgM8190_3GHzIF_12GSPS_0121.
% target: 'final' x_bb_long of size [2047x12x64]
% 1. get x_bb_long and scale it to make the total energy to be 64 (64
% repetition).
% 2. get ideal impulse response using the target signal and scale the ideal
% impulse response; make the total energy to be 1.
% 2. get the pre-distortion filter using the ideal impulse response and
% 'jk_mean_ir_12GHz'.
awgM8190_3GHzIF_12GSPS_0121_60GHz; % Line 10 (clear all;) is removed temporarily.
% x_bb_long: (2047x12)x64
% Samp_mkr_repetition: 64
% fs_awg: 12 (Gsps)
% fsymb_awg: 1 (Gsps)
% samp_per_symb = fs_awg/fsymb_awg: 12

x_bb_before_mixer = x_bb_long;
[pn_theory_at_rx,T] = pn(11,[11,8,5,2],samp_per_symb); % 2047x12 samples
pn_theory_at_rx = round(pn_theory_at_rx);
pn_theory_at_rx = reshape(pn_theory_at_rx,numel(pn_theory_at_rx),1);    % column vector
pn_theory_at_rx_long = repmat(pn_theory_at_rx,Samp_mkr_repetition,1);   % [2047x12x64 1]
% The next line has been modified from 'pre_distortion_filter_design_60GHz_v2.m' to remove the DC component.
ideal_ir_fft = fft(x_bb_before_mixer-mean(x_bb_before_mixer)).*conj(fft(pn_theory_at_rx_long));
% FILTER CALCULATION INPUT 1: ideal_ir_fft
ideal_ir_time = ifft(ideal_ir_fft);
%*** average power of ideal_ir_time = mean(abs(ideal_ir_time).^2)
avgPW_ideal_ir_time = mean(abs(ideal_ir_time).^2);

% Debug
if (1)
    ts_ideal_ir = 1/samp_per_symb;  % ns
    n_ideal_ir = numel(ideal_ir_time);
    ideal_ir_t_range = (1:n_ideal_ir)*ts_ideal_ir;
    ideal_ir_f_range = [(-n_ideal_ir/2):(n_ideal_ir/2-1)]/n_ideal_ir*samp_per_symb; % GHz
    figure(121);
    subplot(211);
    plot(ideal_ir_t_range,10*log10(abs(ideal_ir_time)));
    xlabel('time(ns)');
    ylabel('mag(ideal impulse response) (dB)');
    subplot(212);
    plot(ideal_ir_f_range,10*log10(fftshift(abs(ideal_ir_fft))));
    xlabel('freq(GHz)');
    ylabel('|fft(ideal impulse response)| (dB)');
end

%% fft filter calculation
jk_mean_ir_12GHz_long = repmat(jk_mean_ir_12GHz,Samp_mkr_repetition,1); % [2047x12x64 1]
%*** average power of jk_mean_ir_12GHz_long = mean(abs(jk_mean_ir_12GHz_long).^2)
avgPW_jk_mean_ir_12GHz_long = mean(abs(jk_mean_ir_12GHz_long).^2);
jk_mean_ir_12GHz_long = jk_mean_ir_12GHz_long*sqrt(avgPW_ideal_ir_time/avgPW_jk_mean_ir_12GHz_long);
jk_mean_ir_fft_long = fft(jk_mean_ir_12GHz_long);
% replace zeros with the next min value.
jk_mean_ir_fft_long(jk_mean_ir_fft_long==0) = 1e5;
min_jk_mean_ir_fft_long = min(jk_mean_ir_fft_long);
jk_mean_ir_fft_long(jk_mean_ir_fft_long==1e5) = min_jk_mean_ir_fft_long;
% FILTER CALCULATION INPUT 2: jk_mean_ir_fft_long

test_filter_fft = ideal_ir_fft./jk_mean_ir_fft_long;
% replace zeros with the next min value.
test_filter_fft(test_filter_fft==0) = 1e5;
min_test_filter_fft = min(test_filter_fft);
test_filter_fft(test_filter_fft==1e5) = min_test_filter_fft;

test_result_fft = jk_mean_ir_fft_long.*test_filter_fft;
test_result_time = ifft(test_result_fft);

%pre_filter = test_filter_fft;
%save_str = ['pre_filter_M16_',datestr(now,'mmdd')];
%save(save_str,'pre_filter');

% Debug
if (1)
    figure(122);
    subplot(221);
    plot(ideal_ir_t_range,10*log10(abs(ideal_ir_time)));
    hold on;
    plot(ideal_ir_t_range,10*log10(abs(test_result_time)),'r.:');
    hold off;
    grid on;
    xlabel('time(ns)');
    ylabel('magnitude (dB)');
    legend('ideal IR','test IR');
    
    subplot(223);
    plot(ideal_ir_t_range,angle(ideal_ir_time));
    hold on;
    plot(ideal_ir_t_range,angle(test_result_time),'r.:');
    hold off;
    grid on;
    xlabel('time(ns)');
    ylabel('phase');
    legend('ideal IR','test IR');
    
    subplot(222);
    plot(ideal_ir_f_range,10*log10(fftshift(abs(ideal_ir_fft))));
    hold on;
    plot(ideal_ir_f_range,10*log10(fftshift(abs(test_result_fft))),'r.:');
    hold off;
    xlim([-2 2]);
    %ylim([-10 100]);
    xlabel('freq(GHz)');
    ylabel('|fft| (dB)');
    legend('ideal IR','test IR');
    
    subplot(224);
    plot(ideal_ir_f_range,angle(fftshift(abs(ideal_ir_fft))));
    hold on;
    plot(ideal_ir_f_range,angle(fftshift(abs(test_result_fft))),'r.:');
    hold off;
    xlim([-2 2]);
    xlabel('freq(GHz)');
    ylabel('phase');
    legend('ideal IR','test IR');
end

%%
% Debug
% original transmitted signal 'x_bb_long'
% vs. 
% distorted signal x_bb_distorted
x_bb_distorted_fft = fft(x_bb_before_mixer).*test_filter_fft;

figure(123);
subplot(211);
plot(ideal_ir_t_range,(real(x_bb_before_mixer)));
hold on;
plot(ideal_ir_t_range,(real(ifft(x_bb_distorted_fft))),'r-');
hold off;
legend('org','distort');
xlim([270 360]);

x_distorted_mixer = ifft(x_bb_distorted_fft).*x_local;
x_distorted_mixer = real(x_distorted_mixer);
figure(123);
subplot(212);
plot(ideal_ir_t_range,x_mixer);
hold on;
plot(ideal_ir_t_range,x_distorted_mixer,'r-');
hold off;
legend('org','distort');
xlim([270 360]);