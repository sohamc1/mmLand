% Last version finalized by Thomas H King, reorganized by Jizheng (Daniel)
% He, jizheng4@illinois.edu, 11/10/2022.

%% Micro-displacement Extraction Pipeline.
% - Parameters:
%   - start_and_end_delay: integer, number of frames elapsed where
%           the tag is held still at the start & end of one measurement.
%   - exptDir: directory of all experiment folder
%   - saveDir: directory for saving results. Leave empty for plotting.
% - Returns: TODO
%distances(experimentno, fileNo, tagNo, sampleNo) of estimated distances moved. To get a single
%estimate, call mean(distances, 2)
%GTErrors1 is distanes - errors
%GTErrors2 is distances + errors, use either depending on motion stage direction

%% Organization of Experiments/Data files
%  exptDir          //  Outside directory for all experiments
%  |- experiment1   //  Each experiment has its own subfolder
%  |  |- dataFile1  //  and within are separate datafiles with slight
%  |  |- dataFile2  //  differences in displacement magnitude/direction.
%  |  |- ...
%  |  |- GTFile     //  There is also a ground truth file.
%  |- experiment2 
%  |  |- ...

function [distances, GTErrors1, GTErrors2, unwrappedDist] = microDisp_pipeline(start_and_end_delay, exptDir, saveDir)
arguments
    start_and_end_delay (1,1) double {mustBeInteger} = 50;
    exptDir (1,1) string {mustBeNonzeroLengthText} = "/Data/20230219-veldiff/";
    saveDir (1,1) string = "";
%     exptIndices (1,:) 
end

%% Directory Loads
% Include necessary library directories
CurPath = pwd();
addpath([CurPath, '/DemoRadBasic/DemoRadUsb']);
addpath([CurPath, '/DemoRadBasic/Class']);
FuSca               =   7.5989e-06;

% Experiment files directory
exptDir             =   CurPath + exptDir;
exptFolders         =   dir(exptDir);       %   All subfolders of main expt dir
exptFolders         =   exptFolders(3:end); %   Leave out '.' and '..'
maxRangeList        =   zeros([1 size(exptFolders, 1)]);
minRangeList        =   zeros([1 size(exptFolders, 1)]);

% Parameters
switchPer           =   (625*10^-6);        %   Tag switching rate 
bkndSubtract        =   1;                  %   Enable background subtraction or not
c0                  =   physconst('LightSpeed');

% TODO CLEAN THIS PART
vels = [0.2 0.4 0.6 0.8 1.0 2.0 3.0 4.0 5.0 6.0 8.0 10.0];
for i = 1 : size(exptFolders, 1)
    maxRangeList(i) = 25;
    minRangeList(i) = 15;
end

unwrappedDist = {};
GTErrors1 = {};
GTErrors2 = {};
distances = {};

%% Main Processing Loop
for exptIdx         =   1:size(exptFolders, 1) 
    loadGT          =   1;
    
    exptName        =   exptFolders(exptIdx).name;
    RMax            =   maxRangeList(exptIdx);

    if exptName(end-3:end) == ".mat"
        continue;
    end
    dataFiles       =   dir(exptDir + exptName + "/*.mat");
    dataFileCnt     =   size(dataFiles, 1) - 1; % -1 to leave out GT file

    startFrame      =   1;
    modDuty         =   50;

    selRAll         =   [];
    velocityAll     =   [];
    corrAll         =   [];
    rangeAll        =   [];
    
    if loadGT
        fileMax     =   size(dataFiles, 1);
        GTAll       =   load(exptDir + exptName + "/" + dataFiles(fileMax).name);
    end

    for dataIdx     =   1 : dataFileCnt
        realFFT = []; imagFFT = []; absFFT = []; allFFT = [];
        tagImagAll = []; tagRealAll = []; tagMagAll = []; tagAllAll = [];

        %% Basic parameter retrieval and creation
        %  See visualize.m for config parameter details
        disp("Currently Processing: " + dataFiles(dataIdx).name);
        RadData     =   load(exptDir + exptName + "/" + dataFiles(dataIdx).name);
        DataAll     =   RadData.Data;
        RadData.dtime.TimeZone = 'America/New_York';
        numFrame    =   size(DataAll, 2);

        Cfg         =   RadData.Cfg;
        if (length(Cfg.Seq) == 1)
            numTx   =   1;
        else
            numTx   =   2;
        end

        N           =   floor(RadData.N);
        fs          =   RadData.fs; % sampling freq
        numRx       =   RadData.NrChn;
        NrChn       =   numRx * numTx;
        CalDat      =   RadData.CalDat;
        DataAll(1:N:(size(DataAll,1)), :, :) = []; % Remove the chirp number
        modF        =   1 ./ (2*switchPer);

        Win2D           =   repmat(hanning(N-1), 1, Cfg.FrmMeasSiz, NrChn);
        ScaWin          =   sum(Win2D(:,1,1));
        NFFT            =   2^10;
        NFFTVel         =   2^8;

        kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;
        fc              =   (Cfg.fStop + Cfg.fStrt)/2;
        lambda          =   c0/fc;
        vRange          =   (0:NFFT-1).'./NFFT.*fs.*c0/(2.*kf);

        RMin            =   minRangeList(exptIdx);
        RMax            =   min(RMax, ((Cfg.N/Cfg.TRampUp)*c0) / (4*kf));
%         RMin = 0.2; RMax = 15; % Temp Values for Debugging Purposes

        [~, RMinIdx]    =   min(abs(vRange - RMin));
        [~, RMaxIdx]    =   min(abs(vRange - RMax));
        vRangeExt       =   vRange(RMinIdx:RMaxIdx);

        WinVel          =   hanning(Cfg.FrmMeasSiz);
        ScaWinVel       =   sum(WinVel);
        WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

        % Calibration data
        mCalData        =   permute(repmat(CalDat(1:NrChn), 1, Cfg.FrmMeasSiz, N-1), [3 2 1]);

        tagRangeWidth   =   1.0; % Length of one side of a tag's signature in a range domain in m (approximate for now)
        vRangeWidth     =   round(tagRangeWidth*NFFT/fs/c0*(2*kf) * 125/N); % Convert above measurement into a number of bins

        RPExtAll = zeros(numFrame, numel(vRangeExt), Cfg.FrmMeasSiz, NrChn);
        RDAll = zeros(numFrame, numel(vRangeExt), NFFTVel, NrChn);
    
        %% Initial Processing and FFTs
        for MeasIdx = 1 : numFrame - startFrame+1
            Data        =  reshape(squeeze(DataAll(:,MeasIdx,:)), N-1, [], NrChn);
            
            % Background subtraction with first chirp of each frame
            if bkndSubtract
                Data    =   Data - Data(:, 1, :); %#ok<UNRCH> 
            end

            % Calculate range profile including calibration
            RP          =   2*fft(Data.*Win2D.*mCalData, NFFT, 1).*FuSca/ScaWin;
            RPExt       =   RP(RMinIdx:RMaxIdx, :, :);
            RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RPExtAll(MeasIdx-startFrame+1, :, :, :) = RPExt;
            RDAll(MeasIdx-startFrame+1, :, :, :) = RD;
        end

        %% Find Range Bins/Do Matched Filtering       
        for MeasIdx = 1 : numFrame-1+startFrame
            [selR, range, velocity, corr] = matched_filtering(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, squeeze(RDAll(MeasIdx,:,:,:)));
            selRAll(dataIdx, MeasIdx-startFrame+1, :) = selR; 
            rangeAll(dataIdx, MeasIdx-startFrame+1, :) = range;
            velocityAll(dataIdx, MeasIdx-startFrame+1, :) = velocity;
            corrAll(dataIdx, MeasIdx-startFrame+1, :) = corr;
        end

        selROverall(dataIdx, :) = round(median(selRAll(dataIdx, :, :), 2));
        velocityOverall(dataIdx, :) = mean(velocityAll(dataIdx, :, :), 2);
        
        %% Find Angle from Radar
        %  Temporarily deleted, see microDisp_pipeline_uncleaned.m

        %% Calculate Tag FFTs for Phase Extraction
        for tagNum = 1 : length(modF)
            expDopFinals = calculateFilter_wMobility(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc(1), velocityOverall(dataIdx, tagNum));

            [~, indxAll] = max(expDopFinals(1, 1:floor(NFFTVel/2), :), [], 2);
            indx = squeeze(indxAll);
            for MeasIdx = 1:1:numFrame-startFrame+1
                %disp(MeasIdx);
                indsAll = selROverall(dataIdx,tagNum) + (-vRangeWidth:vRangeWidth);
                indsUse = indsAll(indsAll >= 1 & indsAll <= size(RPExtAll, 2));
                RPExt = squeeze(RPExtAll(MeasIdx-startFrame+1, indsUse, :,  :));
                RPExt = permute(RPExt, [2 1 3]);

                realFFT(:, :, :) = fft(real(RPExt), NFFTVel, 1);
                imagFFT(:, :, :) = fft(imag(RPExt), NFFTVel, 1);
                absFFT(:, :, :) = fft(abs(RPExt), NFFTVel, 1);
                allFFT(:, :, :) = fft((RPExt), NFFTVel, 1);

                nB = 2;
                tagImagAll(dataIdx, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(imagFFT(max(indx-nB, 1):min(indx+nB, end), :, :), [], 1));
                tagRealAll(dataIdx, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(realFFT(max(indx-nB, 1):min(indx+nB, end), :, :), [], 1));
                tagMagAll(dataIdx, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(absFFT (max(indx-nB, 1):min(indx+nB, end), :, :), [], 1));
                tagAllAll(dataIdx, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(allFFT(max(indx-nB, 1):min(indx+nB, end), :, :), [], 1));
            end
        end

        %% Extract Tag Phase from FFTs
        for tagNum = 1 : length(modF)
            measIndices = 1:numFrame-startFrame+1;
            realFFTAll = squeeze((tagRealAll(dataIdx, measIndices, tagNum, 1:size(tagRealAll, 4), :)));
            imagFFTAll = squeeze((tagImagAll(dataIdx, measIndices, tagNum, 1:size(tagRealAll, 4), :)));
            magFFTAll  = squeeze((tagAllAll(dataIdx, measIndices, tagNum, 1:size(tagRealAll, 4), :)));

            scaledReal = abs(realFFTAll) ./ abs(magFFTAll);
            scaledImag = abs(imagFFTAll) ./ abs(magFFTAll);

            phi = (0 : pi/4000 : pi); % 43 = size of distance btwn two peaks
            x = 1 : size(magFFTAll, 2);

            % Quick and dirty frequency extraction
            freqIndex = [];

            for MeasIdx = 1 : floor((numFrame-startFrame+1)/10) : numFrame-startFrame+1 % select 10 representative samples
                for rxAntennaIdx = 2:3 % across the two central antennas
                    [~, fIdx1] = findpeaks(scaledReal(MeasIdx, :, rxAntennaIdx), 'MinPeakProminence', 0.5);
                    [~, fIdx2] = findpeaks(scaledImag(MeasIdx, :, rxAntennaIdx), 'MinPeakProminence', 0.5);
                    per1temp = diff(fIdx1);
                    per2temp = diff(fIdx2);
                    freqIndex = [freqIndex, per1temp, per2temp];
                end
            end
            period = trimmean(freqIndex, 10);

            triWave = repmat(triang(length(x))', [size(scaledReal, 1)  1 4]);
            %for(MeasIdx = 1:250)%numFrame-startFrame+1)
            autoCorr2a = [];
            autoCorr2b = [];
            for phiInd = 1 : length(phi)
                sinW = repmat(abs(sin(x*pi/period + phi(phiInd))), [size(scaledReal, 1) 1 4]);
                cosW = repmat(abs(cos(x*pi/period + phi(phiInd))), [size(scaledReal, 1) 1 4]);

                autoCorr2a(:, phiInd, :) = (sum(((scaledReal(:, :, :))) .* sinW .* triWave, 2));
                autoCorr2b(:, phiInd, :) = (sum(((scaledImag(:, :, :))) .* cosW .* triWave, 2));
            end

            autocorr = autoCorr2a(:, :, :)+autoCorr2b(:, :, :);
            [~, I]   = max(autocorr, [], 2);
            angles   = squeeze(phi(I));

            % Do some comparison of angle as well?
            unwrappedAngle2 = unwrapPi(angles);
            unwrappedDist{exptIdx}{dataIdx}{tagNum} = 1*lambda/(4*pi) * unwrappedAngle2;
        end

        %% Compare displacements to GT
        numSamples = 50;

        if (loadGT)
            for tagNum = 1:length(modF)
                gtError1 = zeros(numSamples, 1);
                gtError2 = zeros(numSamples, 1);
                distance = zeros(numSamples, 1);
                for sample = 1:numSamples
                    randindsStart = randperm(start_and_end_delay,round(start_and_end_delay/2))+5;
                    randindsEnd = numFrame + 1 - randperm(start_and_end_delay,round(start_and_end_delay/2));

                    pp    = (squeeze(unwrappedDist{exptIdx}{dataIdx}{tagNum}));
                    dist2 = squeeze(mean(trimmean(pp(randindsEnd, :), 10) - mean(pp(randindsStart, :), 1), 2));

                    % GTErrors1 or 2 depends on the direction of your motion stage - I'm not sure how your setup works.
                    gtError1(sample) = 1000*dist2 + GTAll.posn(dataIdx);
                    gtError2(sample) = 1000*dist2 - GTAll.posn(dataIdx);
                    distance(sample) = 1000*dist2; %*3.7778; %Some cnst scale
%                     gtdists1(exptIdx, dataIdx, tagNum, sample) = GTAll.posn(dataIdx);
                end
                GTErrors1{exptIdx}{dataIdx}{tagNum} = gtError1;
                GTErrors2{exptIdx}{dataIdx}{tagNum} = gtError2;
                distances{exptIdx}{dataIdx}{tagNum} = distance;
            end
        end
    end

    % If indicated save path, save processing results; otherwise plot all
    % extracted microdisplacements
    if saveDir.strlength > 0
        save(saveDir + string(datetime("now"), "yyMMdd-HHmmss") + "-" + exptName + "-processedData.mat", ...
             "unwrappedDist", "GTErrors1", "GTErrors2");
%         save(saveDir + string(datetime("now"), "yyMMdd-HHmmss") + "-" + exptName + "-processedData.mat", ...
%              "unwrappedDist", "GTErrors1", "GTErrors2", "distances", "gtdists1");
    else
        for dataIdx = 1 : dataFileCnt
            figure(exptIdx)
            if (dataIdx == 1)
                title("Experiment " + exptIdx)
            end
            subplot(1, dataFileCnt, dataIdx)
            for tagNum = 1 : length(modF)
                plot(squeeze(unwrappedDist{exptIdx}{dataIdx}{tagNum}))
            end
            title("jj =  " + dataIdx + " : GT " + (1/100*GTAll.posn(dataIdx)))
        end
    end
end

end

%% Sinc Template Filter Function Generator
% function template = calculateFilter_wMobility(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, vel)
% c0 = physconst('LightSpeed');
% 
% % Time vector, oversampled to get rid of noise
% t = (0 : Cfg.Perd : Cfg.FrmMeasSiz*Cfg.Perd); 
% 
% % Calculate doppler shift in frequency domain
% doppShift = fft(cos(2*pi*(2*vel*fc/c0).*t), NFFTVel, 2);
% doppShift = doppShift(1:NFFTVel/2) ;
% doppShift = doppShift / max(doppShift);
% doppShift(doppShift < 0.6) = 0;
% 
% % Emulate tag modulation square wave in freq. domain, normalized such that
% % each filter has same "amount" of correlation strength
% tagMod = fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2);
% tagMod = abs(tagMod);
% tagMod = tagMod ./ max(tagMod);
% tagMod(tagMod<0.1) = 0;
% tagMod = tagMod ./ sum(tagMod)*100;
% 
% % Convolute sq. wave and doppler shift in freq. domain. Flip doppler shift
% % based on velocity direction.
% if vel > 0
%     doppShiftPos = [doppShift zeros(1, NFFTVel/2)];
%     template = conv(doppShiftPos, tagMod);
%     template = abs(template(1 : NFFTVel));
% elseif vel < 0
%     doppShiftNeg = flip([doppShift  zeros(1, NFFTVel/2)]);
%     template = conv(doppShiftNeg, tagMod);
%     template = abs(template(NFFTVel : end));
% else
%     template = abs(tagMod);
% end
% 
% template = template ./ max(template);
% template(template<0.1) = 0;
% % template(template < 0.3) = 0;
% template = repmat(template, length(vRangeExt), 1);
% end

function template = calculateFilter_wMobility(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, vel)
    c0 = 2.998e8;

    t = (0 : Cfg.Perd : 20*Cfg.Perd*10);
    doppShift           = fft(cos(2*pi*(2*vel*fc/c0).*t), NFFTVel, 2); % Doppler shift of sq. wave

    undopp_square=fft(square(2*pi*(modF.').*t(1:end-1),modDuty), NFFTVel, 2);
    if (vel>0)
        doppShiftPos = [doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)];
        sq_wav = conv(doppShiftPos, undopp_square);
        expdopt = abs(sq_wav(1:NFFTVel));
    elseif (vel < 0)
        doppShiftNeg = flip([doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)]);
        sq_wav = conv(doppShiftNeg, undopp_square);
        expdopt = abs(sq_wav(NFFTVel:end));
    else
        expdopt = abs(undopp_square);
    end

    expdopt = expdopt ./ max(expdopt);
    expdopt(expdopt<0.1) = 0;
    template = repmat(expdopt, length(vRangeExt), 1);
end

%% Correlation-based Match Filtering
function [selR, range, velocity, corr] = matched_filtering(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, RDA)
% Define the signal for correlation matching
% sig dimensions: (frames) * (chirps) * (channels)
sig = squeeze(abs(RDA));

% Normalize the signal and sum 1~3 channels (TODO why?)
normA = sig - min(sig(:));
normSig = normA ./ max(normA(:));
normSig = sum(normSig(:,:,1:3), 3);
% normSig dim: (frames) * (chirps)

% Calculate background noise across frequency bins
doppNoise = sum(sum(normSig, 3), 1);
convW = 17;
doppNoise = conv(doppNoise, triang(convW));
doppNoise = doppNoise(ceil(convW/2) : end-floor(convW/2));
doppNoise = normalize(doppNoise, 2, 'range');

% Default uses values 10 and -0.2
% Mobility case  using values 400*(-0.03)
sigmoid(:) = 1 ./ (1+exp(10*(-0.2+doppNoise)));

vel = -0.2 : 0.05 : 0.2; % for no mobility cases

% Container for results. Dimension (mod freq) * (velocities) * (range bins)
normSumAll = zeros(length(modF), length(vel), length(vRangeExt));

% Sweep each filter across results
for vIdx = 1 : length(vel)
    for fIdx = 1 : length(modF)
        template = calculateFilter_wMobility(Cfg, modF(fIdx), modDuty, NFFTVel, vRangeExt, fc, vel(vIdx));
        template(isnan(template)) = 0;
        template = template ./ (sum(template(1, :))) * 3; % Want to bias away from a matched filter with more peaks/just larger
        template = template .* sigmoid;

        normSum = squeeze(sum(abs(normSig) .* repmat(abs(template), 1, 1, size(normSig, 3)), 2));
        normSumAll(fIdx, vIdx, :) = normSum;
    end
end

corr = zeros(length(modF), 1);
velocity = zeros(length(modF), 1);
selR = zeros(length(modF), 1);

% Find maximum value and record
for tagNum = 1:length(modF)
    [peakIdxs, corr(tagNum)] = findPeaks2D(squeeze(normSumAll(tagNum, :, :)), 1);
    velocity(tagNum) = vel(peakIdxs(1));
    selR(tagNum) = peakIdxs(2);
end

range = vRangeExt(selR);
end