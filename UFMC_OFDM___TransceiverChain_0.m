% =======================================================================
% MatLab Script “UFMC_OFDM___TransceiverChain.m” Terms of Use
% Version: 2014-03-21
%
% This MatLab script was created and is made available free of charge
% by Alcatel-Lucent Deutschland AG, Lorenzstraße 10, 70435 Stuttgart
% (“Alcatel-Lucent”).
%
% You may use this script for any purpose, and you may modify the
% contents of this script and incorporate the script into your own
% scripts in whole or in part.
%
% If you use this script without modification, you must preserve the
% copyright notice identifying Alcatel-Lucent. If you incorporate this
% script into your own scripts, you must include into them a notice
% stating that portions of the script were created by Alcatel-Lucent.
%
% Alcatel-Lucent does not provide any technical support for the script
% or for MatLab itself. For enquiries related to the above, you may
% contact
%     Frank Schaich, e-mail frank.schaich@acatel-lucent.com, or
%     Thorsten Wild, e-mail thorsten.wild@alcatel-lucent.com.
%
% Alcatel-Lucent does not provide any representation or warranty with
% regard to the functionality of this script, and Alcatel-Lucent does
% does not assume any liability for the functionality of this script.
% Furthermore, Alcatel-Lucent does not represent, warrant, guarantee
% or otherwise assume any liability for the fitness of the script for
% any particular purpose and for any consequences the use of this
% script may have.
%
% By making available the script, Alcatel-Lucent does not provide an
% express or implied license as to any of its or its related companies’ 
% patents or patent applications . The granting of rights embodied
% in this notice relates only to the script itself.
%
% MatLab is software licensed separately by The MathWorks, Inc.
%
% ========================================================================
%
% Further technical details on UFMC can be found in the following 
% references (and references therein):
% [1] F. Schaich, T. Wild, Y. Chen , “Waveform contenders for 5G – 
%     suitability for short packet and low latency transmissions”, 
%     accepted for IEEE VTCs’14, Seoul, Korea, April 2014
% [2] V. Vakilian, T. Wild, F. Schaich, S.t. Brink, J.-F. Frigon, 
%     "Universal-Filtered Multi-Carrier Technique for Wireless Systems 
%     Beyond LTE", 9th International Workshop on Broadband Wireless Access 
%    (BWA) @ IEEE Globecom'13, Atlanta, GA, USA, December 2013.




clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simplified UFMC chain:                                           %%%
%%% single-user, no delay, AWGN, BPSK/QPSK, ZF/MF/MMSE and FFT based %%%
%%% detection, UFMC and CP-OFDM                                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Parameter settings %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   burst placing, modulation   %%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PAR.blockShift = 28;% frequency shift of lowest-frequency block in subcarriers
PAR.nPRB = 10;% Allocation width in number of subbands (sub-band width defined in PAR.blockSize below)   
PAR.DataSource = 'QPSK'; % 'BPSK', 'QPSK'

PAR.Constellation{1} = [-1 1 ; 1 0];              % BPSK signal constellation 
PAR.Constellation{2} = [(+1+1i)/sqrt(2) 0 0 ; ... % QPSK signal constellation (Gray mapping)
                        (+1-1i)/sqrt(2) 0 1 ; ...
                        (-1+1i)/sqrt(2) 1 0 ; ...
                        (-1-1i)/sqrt(2) 1 1 ];

PAR.Tx.Flag_UndoFilterResponse = 1; %if flag set to 1: filter response in pass-band is undone in Tx to uniform power distribution beteen subcarriers.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Parameteres for Rx  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PAR.Rx.ZF = 0; %ZF active?  1 --> active, 0 --> not active
PAR.Rx.MF = 0; %MF active?  1 --> active, 0 --> not active
PAR.Rx.MMSE = 0; %MMSE active?   1 --> active, 0 --> not active
PAR.Rx.FFTbasedRx = 1; %FFT based detection active?   1 --> active, 0 --> not active
PAR.Rx.ChanEst = 'viaKnownSymbs'; %'viaKnownSymbs' (ideal yet, i.e. all symbols used as pilots)
PAR.Rx.flag_CFOcomp_on = 1; % CFO compensation (time domain) done in receiver: 1 = yes; 0 = no

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Settings regarding synchronization misalignments per alloc and per layer  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PAR.rCFO = 0.0;%relative carrier frequency offset in subcarrier spacings  

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   General settings  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
PAR.SNR_dB = [0:10];
PAR.NsymbsperTTI = 14;% number of Multi-carrier symbols per TTI (i.e. per simulation drop)
PAR.FFTsize = 1024;
PAR.lFIR = 74; % filter length: 1 means OFDM, >1 uses a Dolph-Chebychev FIR filter
PAR.FilterPar_dB = 40; % sideband attenuation (design parameter of Dolph-Chebychev filters)
PAR.blockSize = 12; % width of subband in number of subcarriers (needs to match to Filterbandwidth)
PAR.CPlength = 73; % length of OFDM CP in samples  
PAR.NTTIs = 10000; % number of TTIs/drops

if (PAR.Tx.Flag_UndoFilterResponse && (PAR.Rx.ZF || PAR.Rx.MF || PAR.Rx.MMSE))
    error('PAR.Tx.Flag_UndoFilterResponse and at least one of the linear receivers is active. PAR.Tx.Flag_UndoFilterResponse is only applicable to FFT based Rx!')
end

%% Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  basic parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% number of samples per multicarrier symbol
if (PAR.lFIR == 1) % OFDM
    PARderived.lMCsym = PAR.FFTsize + PAR.CPlength;
else
    PARderived.lMCsym = PAR.FFTsize + PAR.lFIR -1;
end

switch PAR.DataSource
        case {'BPSK'}
            PARderived.Bit_per_Symbol = 1;
        case {'QPSK'}
            PARderived.Bit_per_Symbol = 2;
end

PARderived.nSNR = length(PAR.SNR_dB);

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Signal generation  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% OFDM CP and CP removal
% matrix for cyclic prefix addition
CPadd = zeros(PAR.CPlength+PAR.FFTsize,PAR.FFTsize);
CPadd(1:PAR.CPlength,(PAR.FFTsize-PAR.CPlength+1):end) = diag(ones(1,PAR.CPlength));
CPadd((PAR.CPlength+1):end,:) = diag(ones(1,PAR.FFTsize));
PARderived.Tx.CPadd=CPadd;
% matrix for cyclic prefix removal
CPrem = zeros(PAR.FFTsize,PAR.CPlength+PAR.FFTsize);
CPrem(:,(PAR.CPlength+1):end) = diag(ones(1,PAR.FFTsize));
PARderived.Rx.CPrem=CPrem;

% Allocation widths
PARderived.nUsedCarr = PAR.nPRB*PAR.blockSize;
% allocated subcarriers
PARderived.allocatedSubcarriers = [1 : PARderived.nUsedCarr] + PAR.blockShift;

% Generation of IDFT spreading matrices carrying the relevant columns of the IDFT matrice with size PAR.FFTsize
% Dimension of the matrices: [PAR.FFTsize x PARderived.nUsedCarr]
PARderived.V = zeros(PAR.FFTsize,PARderived.nUsedCarr);    
for c = 1:PARderived.nUsedCarr   %loop through all allocated subcarriers     
    SubcarrierIndex=PARderived.allocatedSubcarriers(c);
    PARderived.V([1:PAR.FFTsize],c) = exp(2*pi*1i*([1:PAR.FFTsize]-1)*SubcarrierIndex/PAR.FFTsize); %generation of the IDFT vector
end

% final multicarrier modulation matrix T
if PAR.lFIR == 1  % OFDM
    V=PARderived.V;
    T=(1/norm(V))*CPadd*V; %CP-OFDM = IDFT spreading matrices plus CP addition
    PARderived.T = T;
else % UFMC
    f = chebwin(PAR.lFIR,PAR.FilterPar_dB); %Dolph-Chebyshev
    % initialize helper matrices
    F_all = [];
    V_all = zeros(PAR.FFTsize*PAR.nPRB,PARderived.nUsedCarr);
    for iPRB = 1:PAR.nPRB
        % shift to center carrier
        blockShift = PARderived.allocatedSubcarriers(1)-1; %edge of the allocation
        carrierind = blockShift + (PAR.blockSize+1)/2 + (iPRB-1)*PAR.blockSize; % center carrier
        centerFshift = zeros(PAR.lFIR,1);
        for k = 1:PAR.lFIR
            centerFshift(k) = exp(2*pi*1i*(k-1)*carrierind/PAR.FFTsize);
        end
        % frequency-shifted FIR window
        f1 = f.*centerFshift;
        PARderived.Filterresponse_shifted{iPRB}=f1;
        % generate Toeplitz matrix for convolution
        F{iPRB} = toeplitz([f1;zeros(PAR.FFTsize-1,1)],[f1(1),zeros(1,PAR.FFTsize-1)]);
        % stacked Toeplitz matrices implement multicarrier modulation
        F_all = [F_all F{iPRB}];
        % generate expanded IDFT matrix
        V_all( (1+(iPRB-1)*PAR.FFTsize):(iPRB*PAR.FFTsize), ...
        (1+(iPRB-1)*PAR.blockSize):(iPRB*PAR.blockSize)) = ...
                                PARderived.V(:,(1+(iPRB-1)*PAR.blockSize):(iPRB*PAR.blockSize));
    end
    T = F_all*V_all;
    % Final normalized multicarrier modulation matrix    
    TimeDomainSig=T*ones(PARderived.nUsedCarr,1);
    T=T/sqrt(mean(abs(TimeDomainSig).^2)/PARderived.nUsedCarr*PAR.FFTsize);   
    PARderived.T = T;
    %determine FreqResp in pass-band
    TimeDomainSig=T*ones(PARderived.nUsedCarr,1);
    FreqDomSig_oversampled=fft([TimeDomainSig.' zeros(1,2*PAR.FFTsize-length(TimeDomainSig))])/sqrt(PAR.FFTsize);
    FreqDomSig=FreqDomSig_oversampled(1:2:end);
    if PAR.Tx.Flag_UndoFilterResponse
        PARderived.PredistortionResponse=(FreqDomSig(PARderived.allocatedSubcarriers+1)./(mean(abs(FreqDomSig(PARderived.allocatedSubcarriers+1))))); 
    else
        PARderived.PredistortionResponse=ones(1,PARderived.nUsedCarr);
    end
end


%CFO generation (matrix-multiplication in time domain), per Alloc and per Layer
T = PARderived.T; %Modulation matrix, perfectly time and frequency alligned
Gamma = diag(exp((1j*2*pi*PAR.rCFO*(1:PARderived.lMCsym))/PAR.FFTsize));
PARderived.Gamma = Gamma;
PARderived.GT = Gamma*T; % Modulation matrix (Signal) and CFO matrix combined


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Signal detection  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% static receive filters (for AWGN)

V = PARderived.V;
Gamma = PARderived.Gamma;
T = PARderived.T;
GT = PARderived.GT;

if PAR.Rx.ZF %ZF active
    if ((PAR.CPlength ~= 0)&(PAR.lFIR == 1)) % OFDM w/ CP
        if PAR.Rx.flag_CFOcomp_on == 0
            w_ZF = pinv((1/norm(V))*V);
        else
            w_ZF = pinv((1/norm(Gamma(PAR.CPlength+1:end,PAR.CPlength+1:end)*V))*Gamma(PAR.CPlength+1:end,PAR.CPlength+1:end)*V);
        end
        PARderived.w_ZF = w_ZF;  
    else % UFMC
        if PAR.Rx.flag_CFOcomp_on == 0
            w_ZF = pinv(T);
        else
            w_ZF = pinv(GT);
        end
        PARderived.w_ZF = w_ZF; 
    end
end
if PAR.Rx.MF % MF
    % OFDM with CP
    if ((PAR.CPlength ~= 0)&(PAR.lFIR == 1))
        if PAR.Rx.flag_CFOcomp_on == 0
            w_MF = (1/norm(V))*(V)';
        else
            w_MF = ((1/norm(Gamma(PAR.CPlength+1:end,PAR.CPlength+1:end)*V))*Gamma(PAR.CPlength+1:end,PAR.CPlength+1:end)*V)';
        end
        D = zeros(PARderived.nUsedCarr,PARderived.nUsedCarr);
        for ii=1:PARderived.nUsedCarr
            D(ii,ii) = 1/(w_MF(ii,:)*w_MF(ii,:)');
        end
        PARderived.w_MF = D*w_MF;
     else%UFMC
        if PAR.Rx.flag_CFOcomp_on == 0
            w_MF = T';
        else
            w_MF = GT';
        end
        % normalize it
        %row-wise - each carrier must be normalized
        D = zeros(PARderived.nUsedCarr,PARderived.nUsedCarr);
        for ii=1:PARderived.nUsedCarr
            D(ii,ii) = 1/(w_MF(ii,:)*w_MF(ii,:)');
        end
        PARderived.w_MF = D*w_MF;
    end
end
if PAR.Rx.MMSE %MMSE 
    if ((PAR.CPlength ~= 0)&(PAR.lFIR == 1)) % OFDM with CP
        V_normalized=1/norm(V)*V;
        GV = Gamma(PAR.CPlength+1:end,PAR.CPlength+1:end)*V_normalized;        
        for isnr = 1:PARderived.nSNR
            nvar = 1/(10^(0.1*PAR.SNR_dB(isnr))); %perfect knowledge of noise power assumed
            if PAR.Rx.flag_CFOcomp_on == 0
                w_MMSE{isnr} = inv(V_normalized'*V_normalized + nvar*diag(ones(size(V_normalized,2),1)))*V_normalized';
            else                
                w_MMSE{isnr} = inv(GV'*GV + nvar*diag(ones(size(GV,2),1)))*GV';
            end
            PARderived.w_MMSE.SNR{isnr} = w_MMSE{isnr};
        end
    else %UFMC    
        for isnr = 1:PARderived.nSNR
            nvar = 1/(10^(0.1*PAR.SNR_dB(isnr))); %perfect knowledge of noise power assumed
            if PAR.Rx.flag_CFOcomp_on == 0
                w_MMSE{isnr} = inv(T'*T + nvar*diag(ones(size(T,2),1)))*T';
            else
                w_MMSE{isnr} = inv(GT'*GT + nvar*diag(ones(size(GT,2),1)))*GT';
            end
            PARderived.w_MMSE.SNR{isnr} = w_MMSE{isnr};
        end
    end
end

%% Begin of main simulation loops
% main simulation loops over SNR-values and TTIs
for isnr = 1:PARderived.nSNR
    sqrt_nvar = 1/sqrt(10^(0.1*PAR.SNR_dB(isnr)));

%% Tx
    %%%%%%%%%%%%%%%%%%%%%
    % Signal generation %
    %%%%%%%%%%%%%%%%%%%%%
    for iTTI = 1:PAR.NTTIs        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Symbol vector generation %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        switch PAR.DataSource
            case {'BPSK'}  %BPSK
                s_orig = sign(randn(PARderived.nUsedCarr,PAR.NsymbsperTTI));
            case {'QPSK'}  %QPSK
                s_orig = (1/sqrt(2))*(sign(randn(PARderived.nUsedCarr,PAR.NsymbsperTTI))+j*sign(randn(PARderived.nUsedCarr,PAR.NsymbsperTTI)));
        end
        %s_pilots=s_orig; % ideal chanest so far, each data symbol known and used as pilot
        
        %Predistortion to undo FilterResponse
        
        if (PAR.lFIR ~= 1 && PAR.Tx.Flag_UndoFilterResponse && ~(PAR.Rx.ZF || PAR.Rx.MF || PAR.Rx.MMSE))%UFMC
            s = s_orig./repmat(PARderived.PredistortionResponse.',1,PAR.NsymbsperTTI);
            s_pilots = s_orig./repmat(PARderived.PredistortionResponse.',1,PAR.NsymbsperTTI);
        else
            s=s_orig;
            s_pilots=s_orig;
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Transformation to time domain %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                      %delay raus  
        x = PARderived.GT*s; %signal of user of interest including CFO (G)
        x_pilots = PARderived.T*s_pilots; %pilots of user of interest (for ideal chanest, w/o CFO as CFO compensated separately in time domain)
        
        % add noise
        n = sqrt_nvar*(1/sqrt(2))*(randn(PARderived.lMCsym,PAR.NsymbsperTTI)+j*randn(PARderived.lMCsym,PAR.NsymbsperTTI));
        y = x + n; %superimpose layers and add noise
        y_pilots = x_pilots;
        

%% Rx
        % CP removal
        if (PAR.lFIR == 1) % OFDM
            y = PARderived.Rx.CPrem*y;
        end

        %% obtain symbol estimates
        for iNsymbsperTTI=1:PAR.NsymbsperTTI
            % ZF
            if PAR.Rx.ZF %ZF active?
                w_ZF = PARderived.w_ZF;
                s_est_ZF(:,iNsymbsperTTI) = w_ZF*y(:,iNsymbsperTTI);
                % hard decision
                switch PARderived.Bit_per_Symbol
                    case {1}
                        s_HD_ZF(:,iNsymbsperTTI) = sign(real(s_est_ZF(:,iNsymbsperTTI)));
                    case {2}
                        s_HD_ZF(:,iNsymbsperTTI) = (1/sqrt(2))*(sign(real(s_est_ZF(:,iNsymbsperTTI))) + j*sign(imag(s_est_ZF(:,iNsymbsperTTI))));
                end
            end
            % MF
            if PAR.Rx.MF %MF active?
                w_MF = PARderived.w_MF;
                s_est_MF(:,iNsymbsperTTI) = w_MF*y(:,iNsymbsperTTI);
                % hard decision
                switch PARderived.Bit_per_Symbol
                    case {1}
                        s_HD_MF(:,iNsymbsperTTI) = sign(real(s_est_MF(:,iNsymbsperTTI)));
                    case {2}
                        s_HD_MF(:,iNsymbsperTTI) = (1/sqrt(2))*(sign(real(s_est_MF(:,iNsymbsperTTI))) + j*sign(imag(s_est_MF(:,iNsymbsperTTI))));
                end
            end
            % MMSE
            if PAR.Rx.MMSE %MMSE active?
                w_MMSE{isnr} = PARderived.w_MMSE.SNR{isnr};
                s_est_MMSE(:,iNsymbsperTTI) = w_MMSE{isnr}*y(:,iNsymbsperTTI);
                % hard decision
                switch PARderived.Bit_per_Symbol
                    case {1}
                        s_HD_MMSE(:,iNsymbsperTTI) = sign(real(s_est_MMSE(:,iNsymbsperTTI)));
                    case {2}
                        s_HD_MMSE(:,iNsymbsperTTI) = (1/sqrt(2))*(sign(real(s_est_MMSE(:,iNsymbsperTTI))) + j*sign(imag(s_est_MMSE(:,iNsymbsperTTI))));
                end
            end
            % UFMC FFT based detection
            if (PAR.lFIR ~= 1 && PAR.Rx.FFTbasedRx) %UFMC with FFT based detection active
                if PAR.Rx.flag_CFOcomp_on == 1 % undo CFO
                    y(:,iNsymbsperTTI)  = y(:,iNsymbsperTTI).*exp(-(1j*2*pi*PAR.rCFO*(1:PARderived.lMCsym))/PAR.FFTsize).';
                end
                switch PAR.Rx.ChanEst %Frequency domain estimation. Estimation of the frequency response of the filters in pass-band.
                    case {'viaKnownSymbs'} %all QAM symbols used as pilots
                        y_pilots_padded=[y_pilots(:,iNsymbsperTTI).' zeros(1,2048-length(y_pilots))];
                        H_oversampled = fft(y_pilots_padded)/sqrt(PAR.FFTsize);
                        H = H_oversampled(1:2:end);
                        H_Alloc = H(PARderived.allocatedSubcarriers+1)./s_orig(:,iNsymbsperTTI).';
                end

                y_padded=[y(:,iNsymbsperTTI).' zeros(1,2048-length(y(:,iNsymbsperTTI)))];
                s_oversampled=fft(y_padded)/sqrt(PAR.FFTsize);
                s_Rx=s_oversampled(1:2:end); %all subcarriers unequalized
                s_est_FFT(:,iNsymbsperTTI) = s_Rx(PARderived.allocatedSubcarriers+1)./H_Alloc;
                switch PARderived.Bit_per_Symbol
                    case {1} %BPSK
                        s_HD_FFT(:,iNsymbsperTTI) = sign(real(s_est_FFT(:,iNsymbsperTTI)));
                    case {2} %QPSK
                        s_HD_FFT(:,iNsymbsperTTI) = (1/sqrt(2))*(sign(real(s_est_FFT(:,iNsymbsperTTI))) + j*sign(imag(s_est_FFT(:,iNsymbsperTTI))));
                end
                % OFDM FFT based
            elseif (PAR.lFIR == 1 && PAR.Rx.FFTbasedRx) %OFDM with FFT based detection active
                if PAR.Rx.flag_CFOcomp_on == 1 % undo CFO
                    y(:,iNsymbsperTTI)  = y(:,iNsymbsperTTI).*exp(-(1j*2*pi*PAR.rCFO*(PAR.CPlength+1:PARderived.lMCsym))/PAR.FFTsize).';
                end
                s_Rx=(fft(y(:,iNsymbsperTTI))/sqrt(PAR.FFTsize)).';
                %generate 8 times oversampled frequency sig.
                sH_FFT_tmp(:,iTTI)=(fft([y(:,iNsymbsperTTI).' zeros(1,7*PAR.FFTsize)])/sqrt(PAR.FFTsize)).';
                %
                s_est_FFT(:,iNsymbsperTTI) = s_Rx(PARderived.allocatedSubcarriers+1);
                switch PARderived.Bit_per_Symbol
                    case {1} %BPSK
                        s_HD_FFT(:,iNsymbsperTTI) = sign(real(s_est_FFT(:,iNsymbsperTTI)));
                    case {2} %QPSK
                        s_HD_FFT(:,iNsymbsperTTI) = (1/sqrt(2))*(sign(real(s_est_FFT(:,iNsymbsperTTI))) + j*sign(imag(s_est_FFT(:,iNsymbsperTTI))));
                end
            end
        end
        % track MSE and raw SER
        if (PAR.Rx.ZF)
            SE_ZF(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = abs(s_orig-s_est_ZF).^2;
            symErr_ZF(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = sum(s_HD_ZF~=s_orig);
        end
        if (PAR.Rx.MF)
            SE_MF(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = abs(s_orig-s_est_MF).^2;
            symErr_MF(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = sum(s_HD_MF~=s_orig);
        end
        if (PAR.Rx.MMSE)
            SE_MMSE(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = abs(s_orig-s_est_MMSE).^2;
            symErr_MMSE(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = sum(s_HD_MMSE~=s_orig);
        end
        if (PAR.Rx.FFTbasedRx) %UFMC with FFT based detection active
            SE_FFT(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = abs(s_orig-s_est_FFT).^2;
            symErr_FFT(:,(iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = sum(s_HD_FFT~=s_orig);
        end
        iEb((iTTI-1)*PAR.NsymbsperTTI+1:iTTI*PAR.NsymbsperTTI,isnr) = diag(x'*x/(PARderived.nUsedCarr*PARderived.Bit_per_Symbol));
    end
end
%% compute simulation results (SER, EbN0 and MSE)
for isnr = 1:PARderived.nSNR
    % MSE
    if  (PAR.Rx.ZF)
        RES.MSE_ZF(isnr) = mean(mean(squeeze(SE_ZF(:,:,isnr))));
    end
    if  (PAR.Rx.MF)
        RES.MSE_MF(isnr) = mean(mean(squeeze(SE_MF(:,:,isnr))));
    end
    if  (PAR.Rx.MMSE)
        RES.MSE_MMSE(isnr) = mean(mean(squeeze(SE_MMSE(:,:,isnr))));
    end
    if (PAR.Rx.FFTbasedRx)
       RES.MSE_FFT(isnr) = mean(mean(squeeze(SE_FFT(:,:,isnr))));
    end
    % Eb / N0   
    RES.EbN0_dB(isnr) = 10*log10(mean(iEb(:,isnr)))+PAR.SNR_dB(isnr); % SNR is in fact noise level
    % symbol error rate
    if  (PAR.Rx.ZF)
        RES.SER_ZF(isnr)= sum(sum(squeeze(symErr_ZF(:,:,isnr)))) / (PAR.NTTIs*PAR.NsymbsperTTI*PARderived.nUsedCarr);
    end
    if  (PAR.Rx.MF)
        RES.SER_MF(isnr)= sum(sum(squeeze(symErr_MF(:,:,isnr)))) / (PAR.NTTIs*PAR.NsymbsperTTI*PARderived.nUsedCarr);
    end
    if  (PAR.Rx.MMSE)
        RES.SER_MMSE(isnr)= sum(sum(squeeze(symErr_MMSE(:,:,isnr)))) / (PAR.NTTIs*PAR.NsymbsperTTI*PARderived.nUsedCarr);
    end
    if (PAR.Rx.FFTbasedRx)
        RES.SER_FFT(isnr)= sum(sum(squeeze(symErr_FFT(:,:,isnr)))) / (PAR.NTTIs*PAR.NsymbsperTTI*PARderived.nUsedCarr);
    end
end
