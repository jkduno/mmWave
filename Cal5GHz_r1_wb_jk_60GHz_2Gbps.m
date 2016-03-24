%% Program to read and process back to back calibration data
%%Read digitized IF  and data, downconvert and correlate, claculate cal
%%factors, plot and store

%% original flie: Cal5GHz_r1_wb_jk_60GHz.m
%% file location: work and study\mmWave\New System\ReportAndShare\Pre Distortion Filter\PreDistFilter60GHz

%% modified file: Cal5GHz_r1_wb_jk_60GHz_2Gbps.m
%% file location: work and study\mmWave\New System\ReportAndShare\Pre Distortion Filter\PreDistFilter60GHz(2Gbps)
%% Just two things changed from org file: 
% 1. pathname for wbcal018 and wbcal022
% 2. Use gran(Nfile) in wbcal018hdr (or wbcal022hdr), instead of 1.9778
% which was used only for wbcal007.raw.

%% enter fixed parameters
sr=40e9;    %default Guzic sample rate
dt=1/sr;    %sample per sec
br=2.0e9;   %bit rate
pn_bits=2047;       %number of bits per pn word
spb=sr/br;          %samples/bit
npn=pn_bits*spb; %sample points per pn code word
Fif=3e9;    % JK: 2016-3-14
%Fif=20e9;
Fs=sr;
dFs=sr/npn;
Fscale=1e-9*(-sr/2:dFs:(Fs/2)-dFs); %frequency scale in GHz
%pathname= 'c:\Users\pbp\Documents\MATLAB\NISTMFiles\RawData\BacktoBackCalData\';
%pathname= 'C:\Users\jnc3\Documents\MATLAB\BackToBack\';
pathname= 'C:\JK\from Desktop\work and study\mmWave\New System\ReportAndShare\Pre Distortion Filter\PreDistFilter60GHz(2Gbps)\org_files_from_boulder\'; 
%cd pathname
FilterSpec=[pathname '*.mat'];
%% open hdr file or enter data parameters
%hinput=input('(R)ead parameters from hdr file or E(nter),R or E?','s');
hinput='R';
if strcmp(hinput,'R') || strcmp(hinput,'r')
    display('Find hdr file');
    [Filename,PathName,FilterIndex] = uigetfile(FilterSpec,'hdr calibration file');
    load([PathName Filename])
end
if strcmp(hinput,'E') || strcmp(hinput,'e')
    Nfiles=input('Enter the number of files(power levels) to read, Nfiles =');
    FilterSpec=[pathname '*.RAW'];
    display('SELECT first data file')
    [filename,PathName,FilterIndex] = uigetfile(FilterSpec,'raw data file');
    bfn=str2num(filename(4:6)); %#ok<ST2NM>
    for Nfile=1:Nfiles
        display(['Enter data for file = ' num2str(Nfile)]);
        gran(Nfile)=input(['enter granularity(mv/lsb) for file' num2str(Nfile) ' = ']);
        Pin(Nfile)=input(['enter input power, Pin(dBm) for ' num2str(Nfile) ' = ']);
    end
    hrdfile=input('enter name for hdr information file','s');
    save ([PathName hrdfile 'hdr.mat'], 'Pin', 'gran' ,'Nfiles');
end

%targetNfiles = 6;   % jk
targetNfiles = 1;   % JK: 2016-3-14
%targetNfiles = 5;   % jk
%targetNfiles = 8;   % jk
%% open ifraw 5 GHz IF data file
FilterSpec=[pathname '*.raw'];
%for Nfile=1:Nfiles     % jk, rm
for Nfile=targetNfiles:targetNfiles     % jk
    if Nfile == 1
        %display('enter first data file name')
        %[filename,PathName,FilterIndex] = uigetfile(FilterSpec,'raw calibration data files');
        filename=[Filename(1:8) '.raw'];
        %filename='wbcal847.raw';    % M1, Pin(dBm) of -57.56
        fid=fopen([PathName filename],'r');
        [ifraw,count]=fread(fid,inf,'uint8');
        bf=filename(6:8);
        bfn=str2num(bf);
        display(['Reading ' filename]);
        %data(Nfile,:)=ifraw;   % jk, rm
        jk_data = ifraw;        % jk
        %ns(Nfile)=count;       % jk, rm
        jk_ns = count;          % jk
        clear ifraw count;
        fclose(fid);
        %fname(Nfile,:)=filename;   % jk, rm
        jk_fname = filename;        % jk
    end
    if Nfile > 1
        filename=[Filename(1:8) '.raw'];    % jk
        bf=filename(6:8);   % jk
        bfn=str2num(bf);    % jk
        fnum=bfn+Nfile-1;
        if fnum < 10;
        nfilename=[filename(1:6) num2str(fnum) '.raw'];
        fid=fopen([PathName nfilename],'r');
        [ifraw,count]=fread(fid,inf,'uint8');
        display(['Reading ' nfilename]);
        elseif fnum >= 10 && fnum < 100;
        nfilename=[filename(1:5) num2str(fnum) '.raw'];
        fid=fopen([PathName nfilename],'r');
        [ifraw,count]=fread(fid,inf,'uint8');
        display(['Reading ' nfilename]);
        elseif fnum >=100 && fnum <1000
        nfilename=[filename(1:5) num2str(fnum) '.raw'];
        fid=fopen([PathName nfilename],'r');
        [ifraw,count]=fread(fid,inf,'uint8');
        display(['Reading ' nfilename]);
        else
            display(['error fnum>1000 fnum= ' num2str(fnum)]);
        end
        %data(Nfile,:)=ifraw;   % jk, rm
        jk_data = ifraw;        % jk
        %jk_data = circshift(ifraw,-48290);  % jk M1
        %jk_data = circshift(ifraw,-52406);  % jk M2
        %jk_data = circshift(ifraw,-11215);  % jk M3
        %ns(Nfile)=count;       % jk, rm
        jk_ns = count;          % jk
        clear ifraw count;
        fclose(fid);
        %fname(Nfile,:)=nfilename;  % jk, rm
        jk_fname = nfilename;          % jk
    end
end

% if strcmp(hinput,'E') || strcmp(hinput,'e')
% Nfiles=input('Enter the number of files(power levels) to read, Nfiles =');
%     for Nfile=1:Nfiles
%         display(['Enter data for file = ' num2str(Nfile)]);
%         gran(Nfile)=input(['enter granularity(mv/lsb) for file' num2str(Nfile) ' = ']);
%         Pin(Nfile)=input(['enter input power, Pin(dBm) for ' filename ' = ']);
%     end
%     save ([PathName filename(1:6) 'hdr.mat'], 'Pin', 'gran' ,'Nfiles');
% end

%%Determine Processing Parameters
% NImp=number of impulses(code words)/file or data segment
%Nimpulses=round(min(ns)/npn);  % jk, rm
Nimpulses=floor(jk_ns/npn);     % jk
display(['Nimpulses= ' num2str(Nimpulses)]);
% find largest 2^n divisor of Nimp
for  pofn=1:10
    if 2^pofn<=Nimpulses
        nmax=pofn;
    end
end
display(['largest 2^n divisor of Nimpulses(nmax) = ' num2str(nmax) ', Nimpuses= ' num2str(2^nmax)]);
if_V=zeros(2^nmax,npn);         %dimension if_V
ir=zeros(2^nmax,npn);           %dimension ir
pdp=zeros(2^nmax,npn);           %dimension pdp
%apdp=zeros(Nfiles,nmax,npn);   % jk, rm
%sigPW=zeros(Nfiles,nmax);      % jk, rm
%sigdB=zeros(Nfiles,nmax);      % jk, rm
%pdp_iod=zeros(2^nmax,npn);           %dimension pdp % jk, rm
%apdp_iod=zeros(Nfiles,nmax,npn);       % jk, rm
%sigPW_iod=zeros(Nfiles,nmax);          % jk, rm
%sigdB__pdp=zeros(Nfiles,nmax);         % jk, rm

%% process data to uncalibrated PDP'd and APDP's this is for multiple file option add multiple segments in one file later
%downconvert file by file or segment by segment
%[pn_theory,T] = pn(9,[9 4],spc/dfac);
[pn_theory,T]=pn(11,[11,8,5,2],spb); 
%zerofill=zeros(1,npn/2);  %for fft ifv % jk, rm
%apdp_iods=zeros(Nfiles,nmax,npn);  % jk, rm
%apdp=zeros(Nfiles,nmax,npn);       % jk, rm
pw_pn_meas = zeros(2^nmax,1);   % jk
%for Nfile=1:Nfiles                     % jk, rm
clear z
for Nfile=targetNfiles:targetNfiles     % jk
    %for segNo=1:Nseg;
    for impNo = 1:2^nmax     %index raw data into impulses records for one channel from the currnet file/segment/burst
        idx1=1+npn*(impNo-1);
        idx2=npn+npn*(impNo-1);
        %if_V(impNo,:) = (data(Nfile,idx1:idx2) - mean(data(Nfile,idx1:idx2)))*gran(Nfile)*1e-3;    % jk, rm
        if_V(impNo,:) = (jk_data(idx1:idx2) - mean(jk_data(idx1:idx2)))*gran(Nfile)*1e-3;
        %if_V(impNo,:) = (jk_data(idx1:idx2) - mean(jk_data(idx1:idx2)))*1.9778*1e-3;    % JK: 2016-3-14
        z(impNo,:) = if_V(impNo,:).*exp(-1i*(0:(npn-1))*2*pi*Fif*dt);
    end %for impNo
    % Take fourier transform of baseband signal and zero 2nd lobe (indices 512-1533)
    Zfft = (sqrt(2)/npn)*fft(z,[],2);
    %zerofill=zeros(1,npn/2);  %for Nip     % jk, rm
    %zerofill=repmat(zerofill,2^nmax,1);    % jk, rm
    %Zfft(:,npn/4+1:3*npn/4) = zerofill;    % jk, rm
    
    Zfft = fftshift(Zfft,2);    % jk
    dc_index = npn/2+1;         % jk
    BB_BW = 20;              % jk
    srGHz = sr/1e9;             % jk
    BB_freq_indices = dc_index + [(-npn/srGHz*BB_BW/2):(npn/srGHz*BB_BW/2-1)];    % jk
    tmp_Zfft = zeros(size(Zfft));   % jk
    tmp_Zfft(:,BB_freq_indices) = Zfft(:,BB_freq_indices);  % jk
    Zfft = ifftshift(tmp_Zfft,2);    % jk
    
    pn_meas = ifft(Zfft,[],2);
    %pn_meas = ifft(Zfft,[],2)*length(Zfft); % 3-17-2016
    
    for impNo=1:2^nmax       %impulse processing one burst one channel at a time
        pw_pn_meas(impNo) = mean(abs(pn_meas(impNo,:)).^2);
        ir(impNo,:) = ifft(fft(pn_meas(impNo,:)-mean(pn_meas(impNo,:))).*conj(fft(round(pn_theory'))));
        %ir(impNo,:) = ifft(fft(pn_meas(impNo,:)-mean(pn_meas(impNo,:)))/length(pn_meas).*conj(fft(round(pn_theory'))/length(pn_theory)))/numel(pn_theory);
    end   %end for impNo=1:Nimpulses
    pdp_tst=abs(ir).^2;
    pdp=abs(ir.^2);
    %% calculare apdp receive power and calibration factor for all power levels using differnt avegaging 
    % nmax is the max n which give the max 2^n divisor of Nimpulses
    display(['Nfile='  num2str(Nfile) ' Pin(dBm)=' num2str(Pin(Nfile))]);
end     % jk

jk_mean_ir = mean(ir,1);
[~,maxid]=max(10*log10(abs(jk_mean_ir)))
10*log10(mean(mean(pdp))*1000)

