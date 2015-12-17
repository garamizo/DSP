classdef myDSP  
methods(Static)

function [SIGNAL, nowin] = reshape(signal, win_size, overlap_p)
%RESHAPE Reshape vector in window blocks, with overlaps
% FRF Estimation of simulated, real-looking data
%     overlap = .25;
%     Fs = 4096;
%     bsize = 4096;
%     len = 409600;
% 
%     % Simulate data ------------------------------
%     sys = tf(1, [1e-2 5e-1 3e5]);
%     sx = 1e-2;
%     sy = 1e-2 * dcgain(sys);
% 
%     t = (0:len-1).' / Fs;
%     x = filter(ones(10,1)/10, 1, reshape(repmat(randn(len/4,1), [1 4]).', [len 1]));
%     y = lsim(sys, x, t) + sy*randn(size(x));
%     x = x + sx*randn(size(x));
%     % --------------------------------------------
% 
%     blocks = floor(len*(1+overlap)/bsize);
% 
%     for k = 1 : blocks
%         idx = (k-1)*blocks + (1:bsize);
%         xb(:,k) = x(idx);
%         yb(:,k) = y(idx);
%     end
% 
%     win = window(@hann, bsize);
% 
%     GX = fft(bsxfun(@times, xb, win));
%     GX = bsxfun(@times, GX, [1 2*ones(1, bsize-1)]'/(mean(win)*bsize));
% 
%     GY = fft(bsxfun(@times, yb, win));
%     GY = bsxfun(@times, GY, [1 2*ones(1, bsize-1)]'/(mean(win)*bsize));
% 
%     Gxx = mean(conj(GX) .* GX, 2);
%     Gyy = mean(conj(GY) .* GY, 2);
%     Gxy = mean(conj(GX) .* GY, 2);
%     Gyx = mean(conj(GY) .* GX, 2);
% 
%     H1 = Gxy ./ Gxx;
%     H2 = Gyy ./ Gyx;
%     coh = (Gxy .* Gyx) ./ (Gxx .* Gyy);
% 
%     H0 = mean(GY ./ GX, 2);
% 
%     f = linspace(0, Fs, bsize);
%     [mag, pha] = bode(sys, f*2*pi);
% 
%     figure(1)
%     subplot(311); semilogy(f, [squeeze(mag) abs([H1 H2 H0])]); xlim([0 Fs/2])
%     legend('Nominal', 'H1', 'H2', 'H0')
%     subplot(312); plot(f, [deg2rad(squeeze(pha)) angle([H1 H2 H0])]); xlim([0 Fs/2])
%     subplot(313); plot(f, coh); xlim([0 Fs/2]); ylim([0 1])
%     set(gcf, 'Position', [59 82 645 856])
% 
%     figure(2)
%     subplot(211); semilogy(f, sqrt(Gxx)); xlim([0 Fs/2]); ylabel('input')
%     subplot(212); semilogy(f, sqrt(Gyy)); xlim([0 Fs/2]); ylabel('output')
%     set(gcf, 'Position', [14 46 476 451])


    len = size(signal, 1); % signal length
    nocha = size(signal, 2); % number of channels

    if win_size > len
        error('Window size is greater than signal size')
    end
    if overlap_p > 1 || overlap_p < 0
        error('Overlap should range between 0 and 1')
    end

    overlap = round(overlap_p*win_size);
    nowin = floor((len-win_size)/(win_size-overlap) + 1); % # of windows

    SIGNAL = zeros(win_size, nowin, nocha);
    for k = 1 : nocha
        % make overlap multiple of win_size, percentual to absolute
        idx = repmat((1:win_size).', [1 nowin]) + ...
            repmat((0:nowin-1)*(win_size-overlap), [win_size 1]);
        SIGNAL(:,:,k) = reshape(signal(idx,k), [win_size, nowin]);
    end
end 


function Y = discretize( U, N, V )
    % Discretize vector U on the range +-V into N bits

    U = uencode( U, N, V, 'signed');

    U = cast( U, 'double' );
    N = cast( N, 'double' );
    V = cast( V, 'double' );

    Y = interp1( [ -2^(N-1) 2^(N-1)-1 ], [-V V], U );
end

function [blocksize, noc] = optimize_blocksize(F0, Fs, noc_range, Fmin)
    
    min_time = 2/Fmin;
    min_blocksize = min_time/Fs;
    
    nocs = noc_range(1) : noc_range(2);
    approx_blocksizes = nocs * Fs/F0;
    % disregard too small block sizes
    approx_blocksizes(approx_blocksizes < min_blocksize) = [];
    if isempty(approx_blocksizes)
        error('Slower component do not fit window choices')
    end
    blocksizes = round(approx_blocksizes);
    
    [err, idx] = min(abs(approx_blocksizes - blocksizes));
    blocksize = blocksizes(idx(1));
    noc = blocksize * F0/Fs;
    
    if err > 1e-3
        warning(['Bad block size. Round error: ' num2str(err)])
    end
end

function rpm = speed_from_tach(time, tach)
   
    len = length(tach);
    tach = reshape(tach, [len 1]);
    
    thres = (prctile(tach, 75) - prctile(tach, 25))*0.6;
    trigger = [false; tach(2:end) - tach(1:end-1) > thres];
    time_rising = time(trigger);
    
    bad = zeros(size(time_rising)) > 0;
    avg_dt = time_rising(2) - time_rising(1);
    for k = 2 : length(time_rising)
        if time_rising(k)-time_rising(k-1) < 0.7*avg_dt
            bad(k) = true;
        else
            avg_dt = avg_dt*0.5 + 0.5*(time_rising(k)-time_rising(k-1));
        end
    end
    triggerb = find(trigger);
    time_rising = time(triggerb(~bad));
    
    w = 2*pi./(time_rising(2:end) - time_rising(1:end-1));
    w = [w; w(end)];
    
%     plot(time, tach, time(trigger), tach(trigger), 'o', time(triggerb(~bad)), tach(triggerb(~bad)), '+')
    rpm = interp1(time(triggerb(~bad)), w, time, 'linear', 'extrap') * 60/(2*pi);
end

function colorMap(order, rpm, Gyy)
    
    [ORDER, RPM] = meshgrid(order, rpm);
    ACC = sqrt(Gyy);

    figure;
    surf(ORDER, RPM, ACC, 'EdgeColor','None');
    view(2); ylim([min(rpm) max(rpm)]); xlim([min(order) max(order)])
    colorbar; xlabel('Order [1]'); ylabel('shaft speed [RPM]')
end

function [o, rpmOut, Gyy] = TimeTVDFTOrderTracker(signalF, Fs, rpmF)

    if ~iscolumn(rpmF) || ~iscolumn(signalF)
        error('Wrong shapes')
    end
    
    % window parameters
    winsize = 500;                  % window size on time
    overlap = 0;                    % window overlap
    win = window(@hann, winsize);   % time window
    
    % reshape signal and rpm into block sizes
    signal = bsxfun(@times, win, myDSP.reshape(signalF, winsize, overlap));
    rpm = myDSP.reshape(rpmF, winsize, overlap);
    
    % compute time variables
    o = linspace(0, 25, 500);
    
    % ang is the cos/sin argument. a/b is the cos/sin coefficient
    ang = cumsum(bsxfun(@times, rpm, permute(o, [1 3 2])), 1)*2*pi/(60*Fs);
    a = squeeze(mean(bsxfun(@times, signal, cos(ang)), 1))*2/rms(win);
    b = squeeze(mean(bsxfun(@times, signal, sin(ang)), 1))*2/rms(win);
    
    GY = a + 1j*b;                  % linear spectrum
    Gyy = conj(GY) .* GY;           % auto-power
    
    figure
    semilogy(o, sqrt(Gyy))
    xlabel('Order [1]')
    xlim([0 o(end)/2])
    
    rpmOut = mean(rpm, 1);
end

function [o, rpmOut, Gyy] = OrderTVDFTOrderTracker(signalF, Fs, rpmF)

    if ~iscolumn(rpmF) || ~iscolumn(signalF)
        error('Wrong shapes')
    end
    
    % compute order variables
    do = 1/10;
    Omax = 25;
    o = linspace(0, Omax, round(Omax/do));       % order vector
    revs = cumsum(rpmF)/(60*Fs);
     
    a = zeros(floor(revs(end)*do), length(o));
    b = zeros(floor(revs(end)*do), length(o));
    rpmOut = zeros(floor(revs(end)*do), 1);
    for k = 1 : floor(revs(end)*do)
        idx = (k-1)/do < revs & revs < k/do;
        
        win = window(@hann, sum(idx));
        rpm = rpmF(idx);
        signal = signalF(idx) .* win * 2 / rms(win);
        
        ang = cumsum(bsxfun(@times, rpm, o), 1) * 2*pi/(60*Fs);
        a(k,:) = squeeze(mean(bsxfun(@times, signal, cos(ang)), 1)).';
        b(k,:) = squeeze(mean(bsxfun(@times, signal, sin(ang)), 1)).';
        
        rpmOut(k) = mean(rpm);
    end
    
    GY = a + 1j*b;                  % linear spectrum
    Gyy = conj(GY) .* GY;           % auto-power
    
    figure
    semilogy(o, sqrt(Gyy))
    xlabel('Order [1]')
    xlim([0 o(end)/2])
end


function [o, rpm, Gyy] = OrderResamplingOrderTracker(signal, Fs, rpmFull)
    
    len = size(signal, 1);
    
    revs = cumsum(rpmFull/60, 1)/Fs;
    y = interp1(revs, signal, linspace(revs(1), revs(end), len).');
    
    Os = len / revs(end);
    
    [f, rpm, Gyy] = myDSP.TimeResamplingOrderTracker(y, Os, rpmFull);
    xlabel('Order [1]')
    
    o = f * Os/Fs;
end


function [f, rpm, Gyy] = TimeResamplingOrderTracker(signal, Fs, rpmFull)
    
    if ~iscolumn(signal) || ~iscolumn(rpmFull)
        error('input column vectors')
    end
    
    winsize = 500;
    overlap = 0.3;
    
    [yreshaped, nowin] = myDSP.reshape(signal, winsize, overlap);
    
    tmp = myDSP.reshape(rpmFull, winsize, overlap);
    rpm = mean(tmp, 1);

    win = window(@hann, winsize);
    gain = [1; 2*ones(winsize-1, 1)] / (rms(win)*winsize);

    GY = bsxfun(@times, fft(bsxfun(@times, yreshaped, win)), gain);
    Gyy = conj(GY) .* GY;

    f = linspace(0, Fs, winsize).';
    
    F = repmat(f, [1 nowin]);
    RPM = repmat(rpm, [winsize 1]);
    
    figure; plot3(F(:,1:10:end), RPM(:,1:10:end), sqrt(Gyy(:,1:10:end)))
    xlabel('Frequency [Hz]'); ylabel('Shaft Speed [RPM]'); zlabel('Amplitude [m/s^2]')
    xlim([0 Fs/2]); grid on
    set(gca, 'ZScale', 'log')
end

end
end

