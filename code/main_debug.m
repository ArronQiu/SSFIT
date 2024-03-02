clear; clear global; 
% clc

%% Read data
global imgCoarse
global imgFine
global imgMaskValid imgMaskInvalid

cd '.' % cd to the code folder
dirIn = '..\test_data\';

imgCoarse = load([dirIn 'coarse_series.mat']);
imgCoarse= imgCoarse.coarse_series;
imgFineRaw = load([dirIn, 'fine_series.mat']);
imgFineRaw = imgFineRaw.fine_series;
imgMaskInvalidRaw = load([dirIn, 'mask_fine_series.mat']);
imgMaskInvalidRaw = imgMaskInvalidRaw.mask_fine_series;

ind_fine = csvread([dirIn, 'index_fine_series.txt']);

[nlc,nsc,ntc] = size(imgCoarse);
[nl,ns,ntf] = size(imgFineRaw);
imgFine = zeros(nl,ns,ntc,class(imgFineRaw));
imgFine(:,:,ind_fine) = imgFineRaw;
imgMaskInvalid = true(nl,ns,ntc);
imgMaskInvalid(:,:,ind_fine) = imgMaskInvalidRaw;
imgMaskValid = ~imgMaskInvalid;

global IMG_IS_BACKGROUND
IMG_IS_BACKGROUND = sum(imgMaskValid,3) == 0;

clear imgFineRaw imgMaskInvalidRaw

%% Displaying input data
fprintf('Displaying input data...\n')
ca = [0 1e4];
figure; montage(imgCoarse); caxis(ca)
figure; montage(imgFine); caxis(ca)
figure; montage(imgMaskInvalid);

%% Data parameters
global CLASS_FINE
CLASS_FINE = class(imgFine);

global SIZE_COARSE SIZE_FINE
global N_PIX_COARSE N_PIX_FINE N_TIME
global SCALE_FACTOR

SIZE_COARSE = size(imgCoarse);
SIZE_FINE = size(imgFine);
N_TIME = SIZE_COARSE(3);
N_PIX_FINE = SIZE_FINE(1) * SIZE_FINE(2);
N_PIX_COARSE = SIZE_COARSE(1) * SIZE_COARSE(2);
SCALE_FACTOR = SIZE_FINE(1) / SIZE_COARSE(1);

%% DEBUG
imgMaskValid


%% Set hyper-parameters
global MAX_VALUE
MAX_VALUE = 1e4;        % the maximum value of NDVI

% customized hyper-parameters
global WINDOW_SIZE_HALF % half of the window size in base extraction (e.g., if 1, window size = 1*2+1 = 3)
global THRES_P_VAR_CUM  % threshold of accumulated variance in base extraction
WINDOW_SIZE_HALF = 3;
THRES_P_VAR_CUM = 0.8;

global WINDOW_SIZE_HALF_CR % half of the window size in searching reference (e.g., if 1, window size = 1*2+1 = 3)
global THRES_R_CR          % threshold of r in searching reference 
WINDOW_SIZE_HALF_CR = 3;
THRES_R_CR = 0.8;

%% Run
% for acceleration, parallel computation (parfor) is used
if isempty(gcp('nocreate')), parpool(min(4, feature('numcores') - 1)); end

tic
imgCoarse = single(imgCoarse);
imgFine = single(imgFine);
imgFineFusion = run_ssfit(imgCoarse, imgFine, imgMaskValid);
toc

figure; montage(imgFineFusion); caxis(ca);

%% FN
function [imgFineFusion] = run_ssfit(imgCoarse, imgFine, imgMaskValid)

basePixWiseFine = computeBaseLocal(imgCoarse);

[fineFill, weight] = computeFillSimple(imgCoarse, imgFine, imgMaskValid);

imgFineFusion = computeFitting(fineFill, basePixWiseFine, weight);

imgFineFusion = computePostProcess(imgCoarse, imgFine, imgMaskValid, imgFineFusion);

end

