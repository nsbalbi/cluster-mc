function [aIdxRecord] = bundleMCOneSide(nA, N, b, deltaHAA, tetherGridSpacing, options)
%BUNDLEMCONESIDE Monte Carlo simulation of a polymer brush with sticky
% ends. See associated publication SI for algorithm details.
%  
%   BUNDLEMCONESIDE() runs the model with default sample parameters and 
%   saves the resulting video as an .avi. Returns linear indeces of A
%   groups over time.
%
%   BUNDLEMCONESIDE(nA, N, b, deltaHAA, tetherGridSpacing) runs the
%   model with nA A chain ends, N Kuhn segments, Kuhn length b, A-A
%   interaction strength deltaHAA, and a distance between tether points of 
%   tetherGridSpacing.
%
%   BUNDLEMCONESIDE(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
%   name-value pairs that can be used to change default function
%   parameters. See "arguments" for details.

%% Sample Starting Parameters
% 12 kDa PS, DAP-Thy pairs, Data from Santos, 2018
% options.CellLength = 1;  % ~size of a Thy/DAP group
% T = 298;  % RT
% b = 1.8;  % from Rubenstein and Colby
% N = 12000/104.15*1.8/9.5/b = 21.8/b;  % n*l/b with l from b/C_inf from R&C
% deltaHAA < (1/8)*deltaHAB = -475;  % since eight A-A bonds can form per A
% tetherGridSpacing = 2;  % density at interface = density on surf*(r1/r2)^2

arguments
    % Chain Parameters
    nA = 1024;  % number of A groups (must be a square number)
    N double = 10;  % polymer number of Kuhn segments (21.8/b for PS)
    b double = 1.8;  % polymer Kuhn length in nm
    deltaHAA double = -300;  % energy bonus of each A-A neighboring (over kB)
    tetherGridSpacing double = 4;  % number of simulation cells separating tether points

    % Simulation Setup Parameters
    options.NIterations = 300000;  % number of iterations in model
    options.SimSize = 199;  % sim is 2 x SimSize x SimSize, pick an odd number
    options.CellLength = 1;  % length of a cell in nm
    options.T double = 300;  % system temp in K
    options.TetherPlacementType = 'grid';  % 'random' (nA, nB any #) or 'grid' (nA, nB square numbers)
    options.EndPlacementType = 'aboveTether';  % 'aboveTether', 'random', 'packed', or 'bundled'
    options.InitBundleSize = 25;  % intialized bundle size for 'bunded' EndPlacementType

    % Iteration Parameters
    options.IterationMethod = 'random';  
    % 'random' - update a random end each time, 
    % 'all' - update all ends in a random order

    % Equilibration Parameters
    options.TempProfileType = 'const';  % 'const', 'ramp', 'jump', 'anneal', 'spikeTemp', 'spikeEnthalpy'
    options.TempEquibIterations = 500;  % 'ramp'/'jump' - # iterations before T changes, 'anneal'/'spikeTemp' - n iterations after ramp/spikes
    options.TempRampRate double = 0.1;  % 'ramp' - increase in temp per iteration in K
    options.Temp2 double = 2000;  % 'jump' - second temp in K, 'anneal' - starting temp, 'spikeTemp' - spike temp
    options.SpikeFrequency = 100;  % number of iterations between start of temp spikes
    options.SpikeDuty = 10;  % number of frames temp is spiked for
    
    % Plot Parameters
    options.MakeVideo logical = true;  % make a video and save it
    options.MakeTrajectoryPlot logical = false;  % make plot showing chain end trails
    options.VideoTempLabel logical = false;  % show temp at current iteration in corner of vid

    % Export Parameters
    options.MakeCSV logical = false;  % export A positions as a csv
    options.ExportInterval = 1000;  % export A positions every x iterations
end

%% System Variables

density = 1/(tetherGridSpacing^2 * options.CellLength^2);  % chains/nm^2
interfaceSep = 1.78*(density)^(1/3)*N;  % distance between binding interface and tether interface in nm
% MODEL IS INDEPENDENT OF THIS

%% Initialization

tic

aArray = zeros(options.SimSize);  % array in which A's are placed, 1+ = A, 0 = empty
% Array that returns real space coords at a given grid location
coordArray = ones(options.SimSize,options.SimSize,3);  % x = (:,:,1), y = (:,:,2), z = (:,:,3)
[coordArray(:,:,1), coordArray(:,:,2)] = meshgrid(options.CellLength*((1:options.SimSize)-(options.SimSize+1)/2), options.CellLength*((1:options.SimSize)-(options.SimSize+1)/2));
coordArray(:,:,3) = interfaceSep*coordArray(:,:,3);

% Assign tether locations in real space
if (strcmp(options.TetherPlacementType, 'random'))
    
    % Places tether locations randomly with 10% padding on each edge
    % (different for A and B)

    aTetherLocs = 0.8 * options.CellLength * options.SimSize * (rand(nA, 2) - 0.5);
    aTetherLocs = [aTetherLocs, zeros(nA, 1)];

elseif (strcmp(options.TetherPlacementType, 'grid'))
    % Places tether points on an evenly spaced grid (same for A and B)

    tetherGridSize = sqrt(nA);
    [tetherAX, tetherAY] = meshgrid(tetherGridSpacing*options.CellLength*((1:tetherGridSize)-(tetherGridSize+1)/2), ...
                                    tetherGridSpacing*options.CellLength*((1:tetherGridSize)-(tetherGridSize+1)/2));
    tetherAZ = zeros(tetherGridSize);
    aTetherLocs = cat(3, tetherAX, tetherAY, tetherAZ);
    aTetherLocs = reshape(aTetherLocs,[nA 3]);
   
end

% Place A groups in array
if strcmp(options.EndPlacementType, 'random')
    simIndeces = randperm(options.SimSize^2);  % shuffle array of indices
    aIdx = simIndeces(1:nA);  % nth entry is Idx of A_n (A with ID n)
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's at the first nA random indices
    end
elseif strcmp(options.EndPlacementType, 'aboveTether')
    % Approximate tether points as spaces in the sim grid
    [AX, AY] = meshgrid(tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2, ...
                                    tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2);
    aIdx = sub2ind([options.SimSize, options.SimSize],round(AY),round(AX));  % convert to indices
    aIdx = reshape(aIdx,1,[]);
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's at the first nA random indices
    end
elseif strcmp(options.EndPlacementType, 'packed')
    % Approximate tether points as spaces in the sim grid
    [AX, AY] = meshgrid(((1:tetherGridSize)) + 0.5*(options.SimSize - sqrt(nA)) - 0.5, ...
                                    ((1:tetherGridSize)) + 0.5*(options.SimSize - sqrt(nA)) - 0.5);
    aIdx = sub2ind([options.SimSize, options.SimSize],round(AY),round(AX));  % convert to indices
    aIdx = reshape(aIdx,1,[]);
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's at the first nA random indices
    end
elseif strcmp(options.EndPlacementType, 'bundled')
    iBSDim = round(sqrt(options.InitBundleSize));  % initial bundle size dimension
    iBSOptions = [iBSDim, iBSDim;
                  iBSDim - 1, iBSDim;
                  iBSDim - 1, iBSDim - 1;
                  iBSDim - 1, iBSDim + 1;
                  iBSDim + 1, iBSDim;
                  iBSDim + 1, iBSDim + 1];  % list of bundle configuraitons near a square
    % Chooses bundle size closest to desired initial bundle size
    [~, iBSChoice] = min(abs(options.InitBundleSize - (iBSOptions(:,1).*iBSOptions(:,2))),[],"all");
    iBSDimX = iBSOptions(iBSChoice,1);  % initial bundle size x dimension
    iBSDimY = iBSOptions(iBSChoice,2);  % initial bundle size y dim

    iBSRemainderX = mod(sqrt(nA),iBSDimX);  % number of SMBs in remainder in X dim
    iBSRemainderY = mod(sqrt(nA),iBSDimY); 
    iBSNDivsX = floor(sqrt(nA)/iBSDimX);  % number of divisions along X dim
    iBSNDivsY = floor(sqrt(nA)/iBSDimY);
    if (iBSRemainderX == 0)  % move final chunk into the remained if it fits perfectly
        iBSNDivsX = iBSNDivsX - 1; 
        iBSRemainderX = iBSDimX;
    end
    if (iBSRemainderY == 0)
        iBSNDivsY = iBSNDivsY - 1; 
        iBSRemainderY = iBSDimY;
    end

    % Calculating number of grid cells defining tether area
    [AXTemp, ~] = meshgrid(tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2, ...
                                    tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2);
    tetherZoneSize = max(round(AXTemp),[],'all') - min(round(AXTemp),[],'all');

    % Calculate size of spaces between bundles in each dim
    iBSDivSizeX = floor((tetherZoneSize - sqrt(nA))/iBSNDivsX);  
    iBSDivSizeY = floor((tetherZoneSize - sqrt(nA))/iBSNDivsY);

    % Calculate remaining spaces, they are tacked on the end
    iBSRemainderDivX = mod(tetherZoneSize - sqrt(nA),iBSNDivsX);  
    iBSRemainderDivY = mod(tetherZoneSize - sqrt(nA),iBSNDivsY);

    % Array counters
    aArray = zeros(tetherZoneSize);
    idxCount = 1;
    IDCount = 1;

    % Read and write columns and rows recursively into a grid following
    % area divisions calculated above
    for iY = 1:iBSNDivsY
        for jY = 1:iBSDimY
            for xI = 1:iBSNDivsX
                for xJ = 1:iBSDimX
                    aArray(idxCount) = IDCount;
                    idxCount = idxCount + 1;
                    IDCount = IDCount + 1;
                end
                for xJ = 1:iBSDivSizeX
                    idxCount = idxCount + 1;
                end
            end
    
            for xJ = 1:iBSRemainderX
                aArray(idxCount) = IDCount;
                idxCount = idxCount + 1;
                IDCount = IDCount + 1;
            end
            for xJ = 1:iBSRemainderDivX
                idxCount = idxCount + 1;
            end
        end
        for jY = 1:iBSDivSizeY
            idxCount = idxCount + tetherZoneSize;
        end
    end
    
    for iY = 1:iBSRemainderY
        for xI = 1:iBSNDivsX
            for xJ = 1:iBSDimX
                aArray(idxCount) = IDCount;
                idxCount = idxCount + 1;
                IDCount = IDCount + 1;
            end
            for xJ = 1:iBSDivSizeX
                idxCount = idxCount + 1;
            end
        end

        for xJ = 1:iBSRemainderX
            aArray(idxCount) = IDCount;
            idxCount = idxCount + 1;
            IDCount = IDCount + 1;
        end
        for xJ = 1:iBSRemainderDivX
            idxCount = idxCount + 1;
        end
    end
    for iY = 1:iBSRemainderDivY
        idxCount = idxCount + tetherZoneSize;
    end

    % Pad array above to fit it into the middle of the sim
    aArray = padarray(aArray,[floor((options.SimSize - tetherZoneSize)/2), floor((options.SimSize - tetherZoneSize)/2)],'both');
    if (mod(tetherZoneSize,2) == 0)
        aArray = padarray(aArray,[1 1],'pre');
    end
    aIdx = find(aArray)';  % gets A idx from A locations in aArray

end

%% Utility Matrices

% Coord shift array; Center, N, then clockwise
ss = options.SimSize;
shiftArray = [0; -1; -1+ss; ss; 1+ss; 1; 1-ss; -ss; -1-ss];  % index shift array
% Neighboring locations for each shift location (index shift)
pairArrayV2 = [
    -1; -1+ss; ss; 1+ss; 1; 1-ss; -ss; -1-ss;  % center (8)
    -2; -2+ss; -1+ss; ss; -ss; -1-ss; -2-ss;  % N (7)
    -2+ss; -2+2*ss; -1+2*ss; 2*ss; ss; -1; -2;  % NE (7)
    -1+ss; -1+2*ss; 2*ss; 1+2*ss; 1+ss; 1; -1;  % E (7)
    ss; 2*ss; 1+2*ss; 2+2*ss; 2+ss; 2; 1;  % SE (7)
    ss; 1+ss; 2+ss; 2; 2-ss; 1-ss; -ss;  % S (7)
    -ss; 1; 2; 2-ss; 2-2*ss; 1-2*ss; -2*ss; % SW (7)
    -1-ss; -1; 1; 1-ss; 1-2*ss; -2*ss; -1-2*ss;  % W (7)
    -2-ss; -2; -1; -ss; -2*ss; -1-2*ss; -2-2*ss];  % NW (7)

%% Iterate

T1 = options.T;  % initial temp for annealing
HAA1 = deltaHAA;

% Only save once every options.ExportInterval iterations
iterRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, 1);  % record of recorded iterations
iterRecord(1) = 0;
aIdxRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);  % record of group positions over time
aIdxRecord(1,:) = aIdx;
tempRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, 1);
tempRecord(1) = options.T;

progressBarNums = floor(linspace(1,options.NIterations,11));

% Iterate the model NIterations times
for i = 1:options.NIterations
    [aIdx, aArray] = iterate(aIdx, aArray);  % iterate once

    % Print Progess
    if (ismember(i, progressBarNums))
        fprintf('Iteration %i/%i\n',i,options.NIterations);
    end

    % Data Recording
    if (mod(i, options.ExportInterval) == 0)  % record every ExportInterval iterations
        iterRecord(i/options.ExportInterval+1) = i;
        aIdxRecord(i/options.ExportInterval+1,:) = aIdx;  % record new positions
        tempRecord(i/options.ExportInterval+1) = options.T;
    end

    % Temperature Control
    if (strcmp(options.TempProfileType, 'jump'))
        % Change temp to [Temp2] after [TempEquibIterations]
        if (i > options.TempEquibIterations)  
            options.T = options.Temp2;  % change T
        end
    elseif (strcmp(options.TempProfileType, 'ramp'))
        % Ramp temp over time with [TempRampRate] after [TempEquibIterations]
        if (i > options.TempEquibIterations)  
            options.T = options.T + options.TempRampRate;  % incrementT
        end
    elseif (strcmp(options.TempProfileType, 'anneal'))
        % Immediately start ramping temp down from [Temp2] to [T], then
        % equilibriate for [TempEquibIterations]
        options.TempRampRate = (T1 - options.Temp2)/(options.NIterations - options.TempEquibIterations);
        if (i == 1)
            options.T = options.Temp2;
        elseif (i <= options.NIterations - options.TempEquibIterations) 
            options.T = options.T + options.TempRampRate;  % incrementT
        end
    elseif (strcmp(options.TempProfileType, 'spikeTemp'))
        % Spikes to Temp2 for [SpikeDuty] iterations every
        % [SpikeFrequency] iterations up until [TempEquibIterations]
        if ((mod(i, options.SpikeFrequency) < options.SpikeDuty) && i <= (options.NIterations-options.TempEquibIterations))
            options.T = options.Temp2;
        else
            options.T = T1;
        end
    elseif (strcmp(options.TempProfileType, 'spikeEnthalpy'))
        % Spikes enthalpy to 0 for [SpikeDuty] iterations every
        % [SpikeFrequency] iterations up until [TempEquibIterations]
        if ((mod(i, options.SpikeFrequency) < options.SpikeDuty) && i <= (options.NIterations-options.TempEquibIterations))
            deltaHAA = 0;
        else
            deltaHAA = HAA1;
        end
    end

end

toc

%% Plotting

% Video Plot

if (options.MakeVideo)

    figVideo = figure();

    % Vary array size depending on iteration methods
    aXRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);
    aYRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);
    
    for i = 1:nA  % for each particle, covert idx record to coords
        aXRecord(:, i) = coordArray(aIdxRecord(:,i));
        aYRecord(:, i) = coordArray(aIdxRecord(:,i) + options.SimSize^2);
    end
    
    % Plot for video, updated each frame
    hold on
    videoPlotA = scatter(aXRecord(1,:), aYRecord(1,:), 10,  'r', "filled", 'MarkerFaceAlpha',.5);  % scatter A points
    if (options.VideoTempLabel)  % add temperature label if desired
        videoTempText = text(0.65*options.CellLength*options.SimSize/2, 0.8*options.CellLength*options.SimSize/2, ['T: ' num2str(tempRecord(1))] );
    end
    hold off
    xlim([min(coordArray(:,:,1),[],'all'), max(coordArray(:,:,1),[],'all')])
    ylim([min(coordArray(:,:,2),[],'all'), max(coordArray(:,:,2),[],'all')])
    xlabel('X [nm]');
    ylabel('Y [nm]');
    
    writerObj = VideoWriter(['MCOneSide_H' num2str(deltaHAA) '_' options.EndPlacementType],'MPEG-4');  % videoWrite object
    writerObj.FrameRate = 60;
    open(writerObj)  % start writing video
    for i = 1:size(aXRecord, 1)
        % update plot data
        set(videoPlotA, 'XData', aXRecord(i,:), 'YData', aYRecord(i,:));
        if (options.VideoTempLabel)
            set(videoTempText, 'String', ['T: ' num2str(tempRecord(i))]);
        end
        F = getframe(figVideo);  % capture frame
        writeVideo(writerObj,F);  % record frame
    end
    close(writerObj);  % stop writing video

end

%% Export Data

if (options.MakeCSV)
    writematrix(aIdxRecord, ['AIdxRecord_' num2str(nA) '_' num2str(abs(deltaHAA)) '_' num2str(options.T) '_' num2str(N) '_' num2str(b) '_' num2str(tetherGridSpacing) '_' num2str(options.InitBundleSize) '_' num2str(options.ExportInterval) '.csv']);  % Export record of A idxs
end

%% Functions

    % Iteration shell function, handles different iteration methods
    function [aIdxOut, aArrayOut] = iterate(aIdxIn, aArrayIn)
    
        % Record input variables
        aIdxOut = aIdxIn;
        aArrayOut = aArrayIn;

        if (strcmp(options.IterationMethod, 'random'))                

            randA = randi(nA);  % current A ID (idx from OG list)
            [aIdxOut, aArrayOut] = iterateMain(aIdxOut, aArrayOut, randA);

        elseif (strcmp(options.IterationMethod, 'all'))

            % Shuffled list of As and Bs
            randAs = randperm(length(aIdxOut));

            % For each A in random order, update state
            for j = 1:nA
        
                randA = randAs(j);  % current A (idx from OG list)
                [aIdxOut, aArrayOut] = iterateMain(aIdxOut, aArrayOut, randA);  
    
            end

        end
        
    end

    % Function that actually does the iterating
    function [aIdxOut, aArrayOut] = iterateMain(aIdxIn, aArrayIn, aID)
        
        % Record input variables
        aIdxOut = aIdxIn;
        aArrayOut = aArrayIn;

        %% Calculating Streching Contribution
                        
        % Coords of A shifted in all directions
        aShift = repmat(aIdxOut(aID), [9 1]) + shiftArray;
        aShift(or(aShift < 1, aShift > options.SimSize^2)) = NaN;  % replace out of bounds with NaN
    
        % Convert rows and columns to indices
        aShiftXIdxs = aShift;
        aShiftYIdxs = aShiftXIdxs + options.SimSize^2;
        aShiftZIdxs = aShiftYIdxs + options.SimSize^2;
        aShiftOOB = isnan(aShiftXIdxs);  % out of bounds A shifts
    
        % Replace NaN with 1 to avoid errors
        aShiftXIdxs(aShiftOOB) = 1; aShiftYIdxs(aShiftOOB) = 1; aShiftZIdxs(aShiftOOB) = 1;
    
        % Get real space coords of each A shift
        aShiftLocs = [coordArray(aShiftXIdxs), coordArray(aShiftYIdxs), coordArray(aShiftZIdxs)];
        aShiftLocs(repmat(aShiftOOB, [1, 3])) = NaN;  % set OOB values to NaN

        % x and y distance between each location and it's tether point
        aShiftDelta = aShiftLocs(:,1:2) - repmat(aTetherLocs(aID,1:2),[9, 1]);
        % G = 3/2 * (DeltaR)^2 / Nb^2
        deltaGAStretch = 3/2*sum(aShiftDelta.^2,2)/N/(b^2);
    
        %% Calculating Collision

        % Evaluate A array at shift locations to check for collision
        aCollision = aArrayOut(aShiftXIdxs);
        aCollision = [1; ~aCollision(2:end)];  % 0 if collision, else 1

        % Default select lowest G move (w/out collision) if GStretch is
        % significantly large (to avoid floating point errors)
        if sum(deltaGAStretch > 300) == 9  % if all are very large
          
            [~, stateChosen] = min(deltaGAStretch - 1000*aCollision);
            % Calculate new idx of A given selected move
            aNewIdx = aIdxOut(aID) + shiftArray(stateChosen);
        
            % Update A array
            aArrayOut(aIdxOut(aID)) = 0;  % set old position to 0
            aArrayOut(aNewIdx) = aID;  % set new position to A idx
            aIdxOut(aID) = aNewIdx;  % set new A idx into master list

            return;

        end

        %% A-A Interactions
        
        % Look for other A connections with each potential A shift
        aShiftPairIdxs = repmat(aIdxOut(aID), [64 1]) + pairArrayV2;
        aShiftPairIdxs(or(aShiftPairIdxs < 1, aShiftPairIdxs > options.SimSize^2)) = NaN;  % replace out of bounds with NaN

        aShiftPairOOB = isnan(aShiftPairIdxs);  % out of bounds A shifts
        % Replace NaN with 1 to avoid errors
        aShiftPairIdxs(aShiftPairOOB) = 1;
        % Set value to 1 if there's another A present
        aShiftPairVals = double(aArrayOut(aShiftPairIdxs) > 0);
        aShiftPairVals(aShiftPairOOB) = 0;  % set OOB values to 0

        % Group values by shift direction
        aShiftPairSum = [
            sum(aShiftPairVals(1:8)); sum(aShiftPairVals(9:15)); sum(aShiftPairVals(16:22));
            sum(aShiftPairVals(23:29)); sum(aShiftPairVals(30:36)); sum(aShiftPairVals(37:43));
            sum(aShiftPairVals(44:50)); sum(aShiftPairVals(51:57)); sum(aShiftPairVals(58:64))];
        % Calculate Delta G contribution of A-A pairing based on well depth
        deltaGAA = deltaHAA/options.T * aShiftPairSum;

        %% Boltzmann Factor
    
        % Calculate Boltzmann weight of each potential state change
        % Multiplies by collision matrix to eliminate any states that have collisions
        boltzmannWeights = aCollision.*exp(-(deltaGAStretch + deltaGAA));
        boltzmannWeights(isnan(boltzmannWeights)) = 0;
    
        % Cumulative sum of Boltzmann weights
        boltSum = cumsum(boltzmannWeights);
        boltSum = boltSum/boltSum(end);
    
        % Pick state randomly using weighted distribution
        randVal = rand();  % random value between 0 and 1
        stateChosen = 1;  % default chosen state
        if randVal > boltSum(1)
            % set chosen state baed on randvals position along the cumsum
            stateChosen = max(find(boltSum < randVal,1,'last'))+1;
        end

        % Calculate new idx of A given selected move
        aNewIdx = aIdxOut(aID) + shiftArray(stateChosen);
    
        % Update A array
        aArrayOut(aIdxOut(aID)) = 0;  % set old position to 0
        aArrayOut(aNewIdx) = aID;  % set new position to A idx
        aIdxOut(aID) = aNewIdx;  % set new A idx into master list

    end

    % Helper function for abShiftPairIdxs
    function val = randValue(arrayIn)
        % Returns a random value from an input array, or NaN if the array is
        % empty
        if (isempty(arrayIn))
            val = NaN;
        else
            val = arrayIn(randi(length(arrayIn), 1));
        end
    end

end

