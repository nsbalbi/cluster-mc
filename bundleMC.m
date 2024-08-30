function [aIdxRecord, bIdxRecord] = bundleMC(nAB, N, b, deltaHAABB, deltaHAB, tetherGridSpacing, options)
%BUNDLEMC Monte Carlo model simulating chain end interactions at the 
% interface between two end-functionalized polymer brushes. See associated
% publication SI for algorithm details.
% 
%   BUNDLEMC() runs the model with default sample parameters and saves the
%   resulting video as an .avi. Returns linear indeces of A and B groups
%   over time.
%
%   BUNDLEMC(nAB, N, b, deltaHAABB, deltaHAB, tetherGridSpacing) runs the
%   model with nAB A and B chain ends, N Kuhn segments, Kuhn length b, A-A
%   and B-B interaction strength deltaHAABB, A-B interaction strength
%   deltaHAB, and a distance between tether points of tetherGridSpacing.
%
%   BUNDLEMC(...,PARAM1,VAL1,PARAM2,VAL2,...) are additional
%   name-value pairs that can be used to change default function
%   parameters. See "arguments" for details.

%% Sample Starting Parameters

% 12 kDa PS, DAP-Thy pairs, Data from Santos, JACS 2018
% cellLength = 1;  % ~size of a Thy/DAP group
% interfaceSep = 10;  % derived from SI exp/sim param table
% T = 298;  % RT
% b = 1.8;  % from Rubenstein and Colby
% N = 12000/104.15*1.8/9.5/b = 21.8/b;  % n*l/b with l from b/C_inf from R&C
% deltaHAB = -31.93 (kJ/mol) *1000/1.38E-23/6.022E23 = -3800;  % from Santos, 2018
% deltaHAA < (1/8)*deltaHAB > -475;  % since eight A-A bonds can form per A
% tetherGridSpacing = 2;  % density at interface = density on surf*(r1/r2)^2

arguments
    % Chain Parameters
    nAB = 400;  % number of A groups (square number if grid)
    N double = 12.1  % polymer number of Kuhn segments (21.8/b for PS)
    b double = 1.8;  % polymer Kuhn length in nm
    deltaHAABB double = -300;  % energy bonus of each A-A and B-B neighboring (over kB)
    deltaHAB double = -100;  % energy bonus of each A-B neighboring (over kB)
    tetherGridSpacing double = 4;  % number of simulation cells separating tether points

    % Simulation Setup Parameters
    options.NIterations = 250;  % number of iterations in model
    options.SimSize = 149;  % sim is 2 x SimSize x SimSize, pick an odd number
    options.CellLength = 1;  % length of a cell in nm
    options.T double = 300  % system temp in K
    options.TetherPlacementType = 'grid';  % 'random' (nA, nB any #) or 'grid' (nA, nB square numbers)
    options.EndPlacementType = 'aboveTether';  % 'aboveTether', 'random', 'packed', or 'bundled'
    options.InitBundleSize = 6;  % intialized bundle size for 'bunded' EndPlacementType

    % Iteration Parameters
    options.IterationMethod = 'sequential';
    % 'sequential' - only one A-B bond per A/B, when updating, all As are updated, then all Bs
    % 'random' - A-B bonds are shared like A-A/B-B, and a random chain end is updated each iteration

    % Temperature Ramp Parameters
    options.TempProfileType = 'const';  % 'const', 'ramp', 'jump', 'anneal'
    options.TempEquibIterations = 500;  % 'ramp'/'jump' - # iterations before T changes, 'anneal' - n iterations after ramp
    options.TempRampRate double = 0.1;  % 'ramp' - increase in temp per iteration in K
    options.Temp2 double = 498;  % 'jump' - second temp in K, 'anneal' - starting temp
    
    % Plot Parameters
    options.MakeVideo logical = true;  % make a video and save it
    options.VideoTempLabel logical = false;  % show temp at current iteration in corner of vid
    
    % Export Parameters
    options.makeCSV logical = false;  % export A,B positions as a csv
    options.ExportInterval = 1;  % export A positions every x iterations
end

%% System Variables

% Chain Parameters
nA = nAB;  % number of A groups
nB = nAB;  % number of B groups

% Interaction Parameters
deltaHAA = deltaHAABB;  % energy bonus of each A-A neighboring (over kB)
deltaHBB = deltaHAABB;  % energy bonus of each B-B neighboring (over kB)

density = 1/(tetherGridSpacing^2 * options.CellLength^2);  % chains/nm^2
interfaceSep = 1.78*(density)^(1/3)*N;  % distance between binding interface and tether interface in nm
% MODEL IS INDEPENDENT OF THIS

%% Initialization

tic

aArray = zeros(options.SimSize); % array in which A's are placed, 1+ = A, 0 = empty
bArray = zeros(options.SimSize); % array in which B's are placed, 1+ = B, 0 = empty
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

    bTetherLocs = 0.8 * options.CellLength * options.SimSize * (rand(nB, 2) - 0.5);
    bTetherLocs = [bTetherLocs, zeros(nB, 1)];

elseif (strcmp(options.TetherPlacementType, 'grid'))

    % Places tether points on an evenly spaced grid (same for A and B)

    tetherGridSize = sqrt(nA);
    [tetherAX, tetherAY] = meshgrid(tetherGridSpacing*options.CellLength*((1:tetherGridSize)-(tetherGridSize+1)/2), ...
                                    tetherGridSpacing*options.CellLength*((1:tetherGridSize)-(tetherGridSize+1)/2));
    tetherAZ = zeros(tetherGridSize);
    aTetherLocs = cat(3, tetherAX, tetherAY, tetherAZ);
    aTetherLocs = reshape(aTetherLocs,[nA 3]);

    bTetherLocs = aTetherLocs;
   
end

% Place A and B groups in array
if strcmp(options.EndPlacementType, 'random')
    simIndecesA = randperm(options.SimSize^2);  % shuffle array of indeces
    simIndecesB = randperm(options.SimSize^2);  % shuffle array of indeces
    aIdx = simIndecesA(1:nA);  % nth entry is ID of A_n
    bIdx = simIndecesB(1:nB);  
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's at the first nA random indeces
        bArray(bIdx(i)) = i;  % place B's at the first nB random indeces
    end
elseif strcmp(options.EndPlacementType, 'aboveTether')
    % Approximate tether points as spaces in the sim grid
    [AX, AY] = meshgrid(tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2, ...
                                    tetherGridSpacing*((1:tetherGridSize)) + 0.5*(options.SimSize - tetherGridSpacing*sqrt(nA)) - tetherGridSpacing/2);
    aIdx = sub2ind([options.SimSize, options.SimSize],round(AY),round(AX));  % convert to indeces
    aIdx = reshape(aIdx,1,[]);
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's above tether locations
    end
    bIdx = aIdx;
    bArray = aArray;
elseif strcmp(options.EndPlacementType, 'packed')
    % Approximate tether points as spaces in the sim grid
    [AX, AY] = meshgrid(((1:tetherGridSize)) + 0.5*(options.SimSize - sqrt(nA)) - 0.5, ...
                                    ((1:tetherGridSize)) + 0.5*(options.SimSize - sqrt(nA)) - 0.5);
    aIdx = sub2ind([options.SimSize, options.SimSize],round(AY),round(AX));  % convert to indeces
    aIdx = reshape(aIdx,1,[]);
    for i = 1:nA
        aArray(aIdx(i)) = i;  % place A's in packed config
    end
    bIdx = aIdx;
    bArray = aArray;
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
    % COULD TACK ON HALF ON EACH SIDE
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
    aIdx = find(aArray);  % gets a IDs from a locations in aArray

    bArray = aArray;
    bIdx = aIdx;

end

%% Utility Matrices

% Coord shift array; Center, N, then clockwise
% shiftArrayRowCol = [0 0; -1 0; -1 1; 0 1; 1 1; 1 0; 1 -1; 0 -1; -1 -1];
ss = options.SimSize;
shiftArray = [0; -1; -1+ss; ss; 1+ss; 1; 1-ss; -ss; -1-ss];  % index shift array
% Neighboring locations for each shift location (index shift) NOT INCLUDING 
% location moved from
pairArray = [
    -1; -1+ss; ss; 1+ss; 1; 1-ss; -ss; -1-ss;  % center (8)
    -2; -2+ss; -1+ss; ss; -ss; -1-ss; -2-ss;  % N (7)
    -2+ss; -2+2*ss; -1+2*ss; 2*ss; ss; -1; -2;  % NE (7)
    -1+ss; -1+2*ss; 2*ss; 1+2*ss; 1+ss; 1; -1;  % E (7)
    ss; 2*ss; 1+2*ss; 2+2*ss; 2+ss; 2; 1;  % SE (7)
    ss; 1+ss; 2+ss; 2; 2-ss; 1-ss; -ss;  % S (7)
    -ss; 1; 2; 2-ss; 2-2*ss; 1-2*ss; -2*ss; % SW (7)
    -1-ss; -1; 1; 1-ss; 1-2*ss; -2*ss; -1-2*ss;  % W (7)
    -2-ss; -2; -1; -ss; -2*ss; -1-2*ss; -2-2*ss];  % NW (7)
% Neighboring locations for each shift location (index shift) INCLUDING
% location moved from AND current location
pairArrayAB = [
    0; -1; -1+ss; ss; 1+ss; 1; 1-ss; -ss; -1-ss;  % center (9)
    -1; -2; -2+ss; -1+ss; ss; 0; -ss; -1-ss; -2-ss;  % N (9)
    -1+ss; -2+ss; -2+2*ss; -1+2*ss; 2*ss; ss; 0; -1; -2;  % NE (9)
    ss; -1+ss; -1+2*ss; 2*ss; 1+2*ss; 1+ss; 1; 0; -1;  % E (9)
    1+ss; ss; 2*ss; 1+2*ss; 2+2*ss; 2+ss; 2; 1; 0;  % SE (9)
    1; 0; ss; 1+ss; 2+ss; 2; 2-ss; 1-ss; -ss;  % S (9)
    1-ss; -ss; 0; 1; 2; 2-ss; 2-2*ss; 1-2*ss; -2*ss; % SW (9)
    -ss; -1-ss; -1; 0; 1; 1-ss; 1-2*ss; -2*ss; -1-2*ss;  % W (9)
    -1-ss; -2-ss; -2; -1; 0; -ss; -2*ss; -1-2*ss; -2-2*ss];  % NW (9)

%% Iterate

T1 = options.T;  % initial temp for annealing

% Only save once every options.ExportInterval iterations
iterRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, 1);  % record of recorded iterations
iterRecord(1) = 0;
aIdxRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);  % record of group positions over time
aIdxRecord(1,:) = aIdx;
bIdxRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nB);  % record of group positions over time
bIdxRecord(1,:) = bIdx;
tempRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, 1);
tempRecord(1) = options.T;

progressBarNums = floor(linspace(1,options.NIterations,11));
if strcmp(options.IterationMethod, 'sequential')
    fprintf('Sequential Iteration\n');
elseif strcmp(options.IterationMethod, 'random')
    fprintf('Random Iteration\n');
end

for i = 1:options.NIterations
    [aIdx, bIdx, aArray, bArray] = iterate(aIdx, bIdx, aArray, bArray);  % iterate once

    % Print Progess
    if (ismember(i, progressBarNums))
        fprintf('Iteration %i/%i\n',i,options.NIterations);
    end

    % Data Recording
    if (mod(i, options.ExportInterval) == 0)  % record every ExportInterval iterations
        iterRecord(i/options.ExportInterval+1) = i;
        aIdxRecord(i/options.ExportInterval+1,:) = aIdx;  % record new positions
        bIdxRecord(i/options.ExportInterval+1,:) = bIdx;  % record new positions
        tempRecord(i/options.ExportInterval+1) = options.T;
    end

    % Temperature Control
    if (strcmp(options.TempProfileType, 'jump'))
        % Change temp after set number of iterations
        if (i > options.TempEquibIterations)  
            options.T = options.Temp2;  % change T
        end
    elseif (strcmp(options.TempProfileType, 'ramp'))
        % Ramp temp over time after set number of iterations
        if (i > options.TempEquibIterations)  
            options.T = options.T + options.TempRampRate;  % incrementT
        end
    elseif (strcmp(options.TempProfileType, 'anneal'))
        % Immediately start ramping temp down from Temp2 to set temp, then
        % equilibriate for TempEquibIterations
        options.TempRampRate = (T1 - options.Temp2)/(options.NIterations - options.TempEquibIterations);
        if (i == 1)
            options.T = options.Temp2;
        elseif (i <= options.NIterations - options.TempEquibIterations) 
            options.T = options.T + options.TempRampRate;  % incrementT
        end
    end

end

toc

%% Plotting

% Video Plot

if (options.MakeVideo)

    figVideo = figure();
    axis("equal");
    set(figVideo,'Color','w');
    aXRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);
    aYRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nA);
    bXRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nB);
    bYRecord = zeros(floor(options.NIterations/options.ExportInterval)+1, nB);
    
    for i = 1:nA  % for each particle, covert ID record to coords
        aXRecord(:, i) = coordArray(aIdxRecord(:,i));
        aYRecord(:, i) = coordArray(aIdxRecord(:,i) + options.SimSize^2);
    end

    for i = 1:nB  % for each particle, covert ID record to coords
        bXRecord(:, i) = coordArray(bIdxRecord(:,i));
        bYRecord(:, i) = coordArray(bIdxRecord(:,i) + options.SimSize^2);
    end
    
    % Plot for video, updated each frame
    hold on
    videoPlotA = scatter(aXRecord(1,:), aYRecord(1,:), 10,  'r', "filled", 'MarkerFaceAlpha',.5);
    videoPlotB = scatter(bXRecord(1,:), bYRecord(1,:), 10,  'b', "filled", 'MarkerFaceAlpha',.5);
    if (options.VideoTempLabel)
        videoTempText = text(0.65*options.CellLength*options.SimSize/2, 0.8*options.CellLength*options.SimSize/2, ['T: ' num2str(tempRecord(1))] );
    end
    hold off
    xlim([min(coordArray(:,:,1),[],'all'), max(coordArray(:,:,1),[],'all')])
    ylim([min(coordArray(:,:,2),[],'all'), max(coordArray(:,:,2),[],'all')])
    xlabel('X [nm]');
    ylabel('Y [nm]');
    
    writerObj = VideoWriter(['MCBothSide_HAA' num2str(deltaHAA) '_HAB' num2str(deltaHAB) '_' options.EndPlacementType],'MPEG-4');  % videoWrite object
    writerObj.FrameRate = 60;
    open(writerObj)  % start writing video
    for i = 1:size(aXRecord, 1)
        set(videoPlotA, 'XData', aXRecord(i,:), 'YData', aYRecord(i,:));  % update plot data
        set(videoPlotB, 'XData', bXRecord(i,:), 'YData', bYRecord(i,:));  % update plot data
        if (options.VideoTempLabel)
            set(videoTempText, 'String', ['T: ' num2str(tempRecord(i))]);
        end
        F = getframe(figVideo);  % capture frame
        writeVideo(writerObj,F);  % record frame
    end
    close(writerObj);  % stop writing video

end

%% Export Data

if (options.makeCSV)
    writematrix(aIdxRecord, ['AIdxRecord_' num2str(nA) '_' num2str(abs(deltaHAA)) '_' num2str(abs(deltaHAB)) '_' num2str(options.T) '_' num2str(N) '_' num2str(b) '_' num2str(tetherGridSpacing) '_' num2str(options.InitBundleSize) '_' num2str(options.ExportInterval) '.csv']);  % Export record of A idxs
    writematrix(bIdxRecord, ['BIdxRecord_' num2str(nA) '_' num2str(abs(deltaHAA)) '_' num2str(abs(deltaHAB)) '_' num2str(options.T) '_' num2str(N) '_' num2str(b) '_' num2str(tetherGridSpacing) '_' num2str(options.InitBundleSize) '_' num2str(options.ExportInterval) '.csv']);  % Export record of A idxs
end

%% Functions

    % Iteration Main Funtion
    function [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterate(aIdxIn, bIdxIn, aArrayIn, bArrayIn)

        % Record input variables
        aIdxOut = aIdxIn;
        bIdxOut = bIdxIn;
        aArrayOut = aArrayIn;
        bArrayOut = bArrayIn;

        if strcmp(options.IterationMethod, 'sequential')

            [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateSequential(aIdxOut, bIdxOut, aArrayOut, bArrayOut);

        elseif strcmp(options.IterationMethod, 'random')

            randType = (rand(1) < 0.5);  % if true, pick A randomly, else pick B
            randAID = randi(nA);
            randBID = randi(nB);

            if (randType)  % iterate randA
                [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateA(randAID, aIdxOut, bIdxOut, aArrayOut, bArrayOut);
            else 
                [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateB(randBID, aIdxOut, bIdxOut, aArrayOut, bArrayOut);
            end

        end

    end

    % Iteration all A's randomly, then all B's randomly. Only one A-B bond
    % per A/B
    function [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateSequential(aIdxIn, bIdxIn, aArrayIn, bArrayIn)
    
        % Record input variables
        aIdxOut = aIdxIn;
        bIdxOut = bIdxIn;
        aArrayOut = aArrayIn;
        bArrayOut = bArrayIn;

        % Arrays for ends with available A-B pairing
        aArrayAvail = aArrayIn;
        bArrayAvail = bArrayIn;
    
        % Shuffled list of As and Bs
        randAs = randperm(length(aIdxOut));
        randBs = randperm(length(bIdxOut));
        
        % For each A & B in random order, update state
        for j = 1:nA

            %% Calculating Streching Contribution
                
            %% Updating an A

            aID = randAs(j);  % current A (idx from OG list)

            % Coords of A shifted in all directions
            aShift = repmat(aIdxOut(aID), [9 1]) + shiftArray;
            aShift(or(aShift > options.SimSize.^2, aShift < 1)) = NaN;  % replace out of bounds with NaN
        
            % Convert 3D index to 1D indeces
            aShiftXIdxs = aShift;
            aShiftYIdxs = aShiftXIdxs + options.SimSize^2;
            aShiftZIdxs = aShiftYIdxs + options.SimSize^2;
            aShiftOOB = isnan(aShiftXIdxs);  % out of bounds A shifts
        
            % Replace NaN with 1 to avoid errors
            aShiftXIdxs(aShiftOOB) = 1; aShiftYIdxs(aShiftOOB) = 1; aShiftZIdxs(aShiftOOB) = 1;
        
            % Get real space coords of each A shift
            aShiftLocs = [coordArray(aShiftXIdxs), coordArray(aShiftYIdxs), coordArray(aShiftZIdxs)];
            aShiftLocs(repmat(aShiftOOB, [1, 3])) = NaN;  % set OOB values to NaN
        
            % End-end dist from A shift to tether loc
            aShiftDelta = aShiftLocs(:,1:2) - repmat(aTetherLocs(aID,1:2),[9, 1]);
            % Gibbs free energy of stretching (over kT) for each shift
            deltaGAStretch = 3/2*sum(aShiftDelta.^2,2)/N/(b^2);

            %% Calculating Collision

            % Evaluate A array at shift locations to check for collision
            aCollision = aArrayOut(aShiftXIdxs);
            aCollision = [1; ~aCollision(2:end)];  % 0 if collision, else 1
    
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
            aShiftPairIdxs = repmat(aIdxOut(aID), [64 1]) + pairArray;
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

            %% A-B Interactions

            % Set value to 1 if there's an available B present (B is
            % not already paired with an A)
            abShiftPairVals = double(bArrayAvail(aShiftPairIdxs) > 0);
            abShiftPairVals(aShiftPairOOB) = 0;  % set OOB values to 0

            % Group values by shift direction
            abShiftPairSum = [
                sum(abShiftPairVals(1:8)); sum(abShiftPairVals(9:15)); sum(abShiftPairVals(16:22));
                sum(abShiftPairVals(23:29)); sum(abShiftPairVals(30:36)); sum(abShiftPairVals(37:43));
                sum(abShiftPairVals(44:50)); sum(abShiftPairVals(51:57)); sum(abShiftPairVals(58:64))];
            abShiftPairSum = double(abShiftPairSum > 0);  % set value to 1 if not equal to 0
            % Calculate Delta G contribution of A-B pairing based on well depth
            deltaGAB = deltaHAB/options.T * abShiftPairSum;

            % For each potential shift with 1 or more A-B pair, save a random
            % B index from one of the pairs
            abShiftPairBIdxs = {
                find(abShiftPairVals(1:8)); 8 + find(abShiftPairVals(9:15)); 15 + find(abShiftPairVals(16:22));
                22 + find(abShiftPairVals(23:29)); 29 + find(abShiftPairVals(30:36)); 36 + find(abShiftPairVals(37:43));
                43 + find(abShiftPairVals(44:50)); 50 + find(abShiftPairVals(51:57)); 57 + find(abShiftPairVals(58:64))};
            abShiftPairBIdxs = cellfun(@randValue, abShiftPairBIdxs);  % CAN PROB OPTIMIZE
            abShiftPairBIdxsNaN = isnan(abShiftPairBIdxs);
            abShiftPairBIdxs(abShiftPairBIdxsNaN) = 1;
            abShiftPairBIdxs = aShiftPairIdxs(abShiftPairBIdxs);
            abShiftPairBIdxs(abShiftPairBIdxsNaN) = NaN;

            %% Boltzmann Factor
        
            % Calculate Boltzmann weight of each potential state change
            % Multiplies by collision matrix to eliminate any states that have collisions
            boltzmannWeights = aCollision.*exp(-(deltaGAStretch + deltaGAA + deltaGAB));
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
        
            % Update B avail array 
            abPairState = deltaGAB(stateChosen) ~= 0;  % true if an A-B pair was made in the chosen state
            if (abPairState)  % if a pair was made with a B
                % Remove B paired with from bArrayAvail
                bArrayAvail(abShiftPairBIdxs(stateChosen)) = 0;
            end 

        end

        for j = 1:nB

            bID = randBs(j); % current B (idx from OG list)

            % Coords of A shifted in all directions
            bShift = repmat(bIdxOut(bID), [9 1]) + shiftArray;
            bShift(or(bShift > options.SimSize.^2, bShift < 1)) = NaN;  % replace out of bounds with NaN
        
            % Convert 3D index to 1D indeces
            bShiftXIdxs = bShift;
            bShiftYIdxs = bShiftXIdxs + options.SimSize^2;
            bShiftZIdxs = bShiftYIdxs + options.SimSize^2;
            bShiftOOB = isnan(bShiftXIdxs);  % out of bounds B shifts
        
            % Replace NaN with 1 to avoid errors
            bShiftXIdxs(bShiftOOB) = 1; bShiftYIdxs(bShiftOOB) = 1; bShiftZIdxs(bShiftOOB) = 1;
        
            % Get real space coords of each B shift
            bShiftLocs = [coordArray(bShiftXIdxs), coordArray(bShiftYIdxs), coordArray(bShiftZIdxs)];
            bShiftLocs(repmat(bShiftOOB, [1, 3])) = NaN;  % set OOB values to NaN
        
            % End-end dist from B shift to tether loc
            bShiftDelta = bShiftLocs(:,1:2) - repmat(bTetherLocs(bID,1:2),[9, 1]);
            % Gibbs free energy of stretching (over kT) for each shift
            deltaGBStretch = 3/2*sum(bShiftDelta.^2,2)/N/(b^2);

            %% Calculating Collision

            % Evaluate B array at shift locations to check for collision
            bCollision = bArrayOut(bShiftXIdxs);
            bCollision = [1; ~bCollision(2:end)];  % 0 if collision, else 1
    
            % Default select lowest G move (w/out collision) if GStretch is
            % significantly large (to avoid floating point errors)
            if sum(deltaGBStretch > 300) == 9  % if all are very large
              
                [~, stateChosen] = min(deltaGBStretch - 1000*bCollision);
                % Calculate new idx of B given selected move
                bNewIdx = bIdxOut(bID) + shiftArray(stateChosen);
            
                % Update A array
                bArrayOut(bIdxOut(bID)) = 0;  % set old position to 0
                bArrayOut(bNewIdx) = bID;  % set new position to B idx
                bIdxOut(bID) = bNewIdx;  % set new B idx into master list
    
                return;
    
            end

            %% B-B Interactions

            % Look for other A connections with each potential B shift
            bShiftPairIdxs = repmat(bIdxOut(bID), [64 1]) + pairArray;
            bShiftPairIdxs(or(bShiftPairIdxs < 1, bShiftPairIdxs > options.SimSize^2)) = NaN;  % replace out of bounds with NaN
    
            bShiftPairOOB = isnan(bShiftPairIdxs);  % out of bounds B shifts
            % Replace NaN with 1 to avoid errors
            bShiftPairIdxs(bShiftPairOOB) = 1;
            % Set value to 1 if there's another A present
            bShiftPairVals = double(bArrayOut(bShiftPairIdxs) > 0);
            bShiftPairVals(bShiftPairOOB) = 0;  % set OOB values to 0
    
            % Group values by shift direction
            bShiftPairSum = [
                sum(bShiftPairVals(1:8)); sum(bShiftPairVals(9:15)); sum(bShiftPairVals(16:22));
                sum(bShiftPairVals(23:29)); sum(bShiftPairVals(30:36)); sum(bShiftPairVals(37:43));
                sum(bShiftPairVals(44:50)); sum(bShiftPairVals(51:57)); sum(bShiftPairVals(58:64))];
            % Calculate Delta G contribution of B-B pairing based on well depth
            deltaGBB = deltaHBB/options.T * bShiftPairSum;

            %% B-A Interactions

            % Set value to 1 if there's an available A present (A is
            % not already paired with a B)
            abShiftPairVals = double(aArrayAvail(bShiftPairIdxs) > 0);
            abShiftPairVals(bShiftPairOOB) = 0;  % set OOB values to 0

            % Group values by shift direction
            abShiftPairSum = [
                sum(abShiftPairVals(1:8)); sum(abShiftPairVals(9:15)); sum(abShiftPairVals(16:22));
                sum(abShiftPairVals(23:29)); sum(abShiftPairVals(30:36)); sum(abShiftPairVals(37:43));
                sum(abShiftPairVals(44:50)); sum(abShiftPairVals(51:57)); sum(abShiftPairVals(58:64))];
            abShiftPairSum = double(abShiftPairSum > 0);  % set value to 1 if not equal to 0
            % Calculate Delta G contribution of A-B pairing based on well depth
            deltaGAB = deltaHAB/options.T * abShiftPairSum;

            % For each potential shift with 1 or more A-B pair, save a random
            % B index from one of the pairs
            abShiftPairAIdxs = {
                find(abShiftPairVals(1:8)); 8 + find(abShiftPairVals(9:15)); 15 + find(abShiftPairVals(16:22));
                22 + find(abShiftPairVals(23:29)); 29 + find(abShiftPairVals(30:36)); 36 + find(abShiftPairVals(37:43));
                43 + find(abShiftPairVals(44:50)); 50 + find(abShiftPairVals(51:57)); 57 + find(abShiftPairVals(58:64))};
            abShiftPairAIdxs = cellfun(@randValue, abShiftPairAIdxs);  % NEED TO CONFIRM THIS WORKS RIGHT, ALSO CAN PROB OPTIMIZE
            abShiftPairAIdxsNaN = isnan(abShiftPairAIdxs);
            abShiftPairAIdxs(abShiftPairAIdxsNaN) = 1;
            abShiftPairAIdxs = bShiftPairIdxs(abShiftPairAIdxs);
            abShiftPairAIdxs(abShiftPairAIdxsNaN) = NaN;

            %% Boltzmann Factor
        
            % Calculate Boltzmann weight of each potential state change
            % Multiplies by collision matrix to eliminate any states that have collisions
            boltzmannWeights = bCollision.*exp(-(deltaGBStretch + deltaGBB + deltaGAB));
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

            % Calculate new idx of B given selected move
            bNewIdx = bIdxOut(bID) + shiftArray(stateChosen);
        
            % Update B array
            bArrayOut(bIdxOut(bID)) = 0;  % set old position to 0
            bArrayOut(bNewIdx) = bID;  % set new position to B idx
            bIdxOut(bID) = bNewIdx;  % set new B idx into master list

            % Update B avail array 
            abPairState = deltaGAB(stateChosen) ~= 0;  % true if an A-B pair was made in the chosen state
            if (abPairState)  % if a pair was made with a B
                % Remove B paired with from bArrayAvail
                aArrayAvail(abShiftPairAIdxs(stateChosen)) = 0;
            end 
        
        end
        
    end

    % Iterate a given A
    function [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateA(aID, aIdxIn, bIdxIn, aArrayIn, bArrayIn)

        % Record input variables
        aIdxOut = aIdxIn;
        bIdxOut = bIdxIn;
        aArrayOut = aArrayIn;
        bArrayOut = bArrayIn;

        %% Calculating Streching Contribution
                
        % Coords of A shifted in all directions
        aShift = repmat(aIdxOut(aID), [9 1]) + shiftArray;
        aShift(or(aShift > options.SimSize.^2, aShift < 1)) = NaN;  % replace out of bounds with NaN
    
        % Convert 3D index to 1D indeces
        aShiftXIdxs = aShift;
        aShiftYIdxs = aShiftXIdxs + options.SimSize^2;
        aShiftZIdxs = aShiftYIdxs + options.SimSize^2;
        aShiftOOB = isnan(aShiftXIdxs);  % out of bounds A shifts
    
        % Replace NaN with 1 to avoid errors
        aShiftXIdxs(aShiftOOB) = 1; aShiftYIdxs(aShiftOOB) = 1; aShiftZIdxs(aShiftOOB) = 1;
    
        % Get real space coords of each A shift
        aShiftLocs = [coordArray(aShiftXIdxs), coordArray(aShiftYIdxs), coordArray(aShiftZIdxs)];
        aShiftLocs(repmat(aShiftOOB, [1, 3])) = NaN;  % set OOB values to NaN
    
        % End-end dist from A shift to tether loc
        aShiftDelta = aShiftLocs(:,1:2) - repmat(aTetherLocs(aID,1:2),[9, 1]);
        % Gibbs free energy of stretching (over kT) for each shift
        deltaGAStretch = 3/2*sum(aShiftDelta.^2,2)/N/(b^2);

        %% Calculating Collision

        % Evaluate A array at shift locations to check for collision
        aCollision = aArrayOut(aShiftXIdxs);
        aCollision = [1; ~aCollision(2:end)];  % 0 if collision, else 1

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
        aShiftPairIdxs = repmat(aIdxOut(aID), [64 1]) + pairArray;
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

        %% A-B Interactions

        % Look for B connections with each potential A shift
        abShiftPairIdxs = repmat(aIdxOut(aID), [81 1]) + pairArrayAB;
        abShiftPairIdxs(or(abShiftPairIdxs < 1, abShiftPairIdxs > options.SimSize^2)) = NaN;  % replace out of bounds with NaN

        abShiftPairOOB = isnan(abShiftPairIdxs);  % out of bounds A shifts
        % Replace NaN with 1 to avoid errors
        abShiftPairIdxs(abShiftPairOOB) = 1;
        % Set value to 1 if there's a B in the B array over that location
        abShiftPairVals = double(bArrayOut(abShiftPairIdxs) > 0);
        abShiftPairVals(abShiftPairOOB) = 0;  % set OOB values to 0

        % Group values by shift direction
        abShiftPairSum = [
            sum(abShiftPairVals(1:9)); sum(abShiftPairVals(10:18)); sum(abShiftPairVals(19:27));
            sum(abShiftPairVals(28:36)); sum(abShiftPairVals(37:45)); sum(abShiftPairVals(46:54));
            sum(abShiftPairVals(55:63)); sum(abShiftPairVals(64:72)); sum(abShiftPairVals(73:81))];
        % Calculate Delta G contribution of A-A pairing based on well depth
        deltaGAB = deltaHAB/options.T * abShiftPairSum;

        %% Boltzmann Factor
    
        % Calculate Boltzmann weight of each potential state change
        % Multiplies by collision matrix to eliminate any states that have collisions
        boltzmannWeights = aCollision.*exp(-(deltaGAStretch + deltaGAA + deltaGAB));
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

    % Iterate a given B
    function [aIdxOut, bIdxOut, aArrayOut, bArrayOut] = iterateB(bID, aIdxIn, bIdxIn, aArrayIn, bArrayIn)

        % Record input variables
        aIdxOut = aIdxIn;
        bIdxOut = bIdxIn;
        aArrayOut = aArrayIn;
        bArrayOut = bArrayIn;

        %% Calculating Streching Contribution
                
        % Coords of A shifted in all directions
        bShift = repmat(bIdxOut(bID), [9 1]) + shiftArray;
        bShift(or(bShift > options.SimSize.^2, bShift < 1)) = NaN;  % replace out of bounds with NaN
    
        % Convert 3D index to 1D indeces
        bShiftXIdxs = bShift;
        bShiftYIdxs = bShiftXIdxs + options.SimSize^2;
        bShiftZIdxs = bShiftYIdxs + options.SimSize^2;
        bShiftOOB = isnan(bShiftXIdxs);  % out of bounds A shifts
    
        % Replace NaN with 1 to avoid errors
        bShiftXIdxs(bShiftOOB) = 1; bShiftYIdxs(bShiftOOB) = 1; bShiftZIdxs(bShiftOOB) = 1;
    
        % Get real space coords of each B shift
        bShiftLocs = [coordArray(bShiftXIdxs), coordArray(bShiftYIdxs), coordArray(bShiftZIdxs)];
        bShiftLocs(repmat(bShiftOOB, [1, 3])) = NaN;  % set OOB values to NaN
    
        % End-end dist from V shift to tether loc
        bShiftDelta = bShiftLocs(:,1:2) - repmat(bTetherLocs(bID,1:2),[9, 1]);
        % Gibbs free energy of stretching (over kT) for each shift
        deltaGBStretch = 3/2*sum(bShiftDelta.^2,2)/N/(b^2);

        %% Calculating Collision

        % Evaluate B array at shift locations to check for collision
        bCollision = bArrayOut(bShiftXIdxs);
        bCollision = [1; ~bCollision(2:end)];  % 0 if collision, else 1

        if sum(deltaGBStretch > 300) == 9  % if all are very large

            [~, stateChosen] = min(deltaGBStretch - 1000*bCollision);
            % Calculate new idx of B given selected move
            bNewIdx = bIdxOut(bID) + shiftArray(stateChosen);
        
            % Update B array
            bArrayOut(bIdxOut(bID)) = 0;  % set old position to 0
            bArrayOut(bNewIdx) = bID;  % set new position to A idx
            bIdxOut(bID) = bNewIdx;  % set new A idx into master list

            return;

        end

        %% B-B Interactions

        % Look for other B connections with each potential B shift
        bShiftPairIdxs = repmat(bIdxOut(bID), [64 1]) + pairArray;
        bShiftPairIdxs(or(bShiftPairIdxs < 1, bShiftPairIdxs > options.SimSize^2)) = NaN;  % replace out of bounds with NaN

        bShiftPairOOB = isnan(bShiftPairIdxs);  % out of bounds B shifts
        % Replace NaN with 1 to avoid errors
        bShiftPairIdxs(bShiftPairOOB) = 1;
        % Set value to 1 if there's another B present
        bShiftPairVals = double(bArrayOut(bShiftPairIdxs) > 0);
        bShiftPairVals(bShiftPairOOB) = 0;  % set OOB values to 0

        % Group values by shift direction
        bShiftPairSum = [
            sum(bShiftPairVals(1:8)); sum(bShiftPairVals(9:15)); sum(bShiftPairVals(16:22));
            sum(bShiftPairVals(23:29)); sum(bShiftPairVals(30:36)); sum(bShiftPairVals(37:43));
            sum(bShiftPairVals(44:50)); sum(bShiftPairVals(51:57)); sum(bShiftPairVals(58:64))];
        % Calculate Delta G contribution of B-B pairing based on well depth
        deltaGBB = deltaHBB/options.T * bShiftPairSum;

        %% A-B Interactions

        % Look for A connections with each potential A shift
        abShiftPairIdxs = repmat(bIdxOut(bID), [81 1]) + pairArrayAB;
        abShiftPairIdxs(or(abShiftPairIdxs < 1, abShiftPairIdxs > options.SimSize^2)) = NaN;  % replace out of bounds with NaN

        abShiftPairOOB = isnan(abShiftPairIdxs);  % out of bounds A shifts
        % Replace NaN with 1 to avoid errors
        abShiftPairIdxs(abShiftPairOOB) = 1;
        % Set value to 1 if there's a B in the B array over that location
        abShiftPairVals = double(aArrayOut(abShiftPairIdxs) > 0);
        abShiftPairVals(abShiftPairOOB) = 0;  % set OOB values to 0

        % Group values by shift direction
        abShiftPairSum = [
            sum(abShiftPairVals(1:9)); sum(abShiftPairVals(10:18)); sum(abShiftPairVals(19:27));
            sum(abShiftPairVals(28:36)); sum(abShiftPairVals(37:45)); sum(abShiftPairVals(46:54));
            sum(abShiftPairVals(55:63)); sum(abShiftPairVals(64:72)); sum(abShiftPairVals(73:81))];
        % Calculate Delta G contribution of A-A pairing based on well depth
        deltaGAB = deltaHAB/options.T * abShiftPairSum;

        %% Boltzmann Factor
    
        % Calculate Boltzmann weight of each potential state change
        % Multiplies by collision matrix to eliminate any states that have collisions
        boltzmannWeights = bCollision.*exp(-(deltaGBStretch + deltaGBB + deltaGAB));
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

        % Calculate new idx of B given selected move
        bNewIdx = bIdxOut(bID) + shiftArray(stateChosen);
    
        % Update B array
        bArrayOut(bIdxOut(bID)) = 0;  % set old position to 0
        bArrayOut(bNewIdx) = bID;  % set new position to B idx
        bIdxOut(bID) = bNewIdx;  % set new B idx into master list

    end

%% Functions

    % Helper function for abShiftPairIDs
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

