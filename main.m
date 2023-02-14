% -----------------------------------------------------------------
% Author: Alvaro Flores-Romero
% Task: Show Principal Component Reconstruction on Spectral Signal
% -----------------------------------------------------------------

clc
clear all
close all
currentDir = pwd;
checker24 = load(strcat(currentDir,"\ColorCheckers\reflec_checker24.mat"));
checker24 = checker24.checker;

% The matrix is 61x24 which means that each column corresponds to each one
% of the 24 patches and the 61 rows are the spectral information
% To put multispectral data the matrix needs to be nxp
% where p is the spectra, so we have to transpose 


% Calculate PCA
%  Each column of eV is a vector of the components
% This matrix has wavelengths over the rows
[eV,eD,eW] = eig(checker24*checker24');
[coeff,score,latent] = pca(checker24');
[U,S,V] = svd(checker24*checker24');

% Plot The First 6 Principal Components
figure
plot(eV(:,56:end));
title('First 6 Principal Components of Color Checker')
legend('PC 6','PC 5 ','PC 4','PC 3','PC 2','PC 1');
xlabel('Wavelength Index');

% Compute the CAF for the first ten eigenvectors

% To get the VAF we need to sum the eigen values from the S matrix
% We extract values different from 0 (which are in the diagonal
eigen = flip(diag(eD));
partialEigen = [];
for i = 1:numel(eigen)
    partialEigen(i) = sum(eigen(1:i))/sum(eigen);
end

% Plot the Variance Accounted For (VAF)
figure;
plot(partialEigen);
title('VAF of Principal Components 1-10 (Rest Ignored on Plot)');
xlim([0,10]);
xlabel('Wavelength Index');
ylabel('VAF [0-1]')

%% Reconstruction
% Reconstruct the 24 spectral reflectances using the first four
% eigenvectors
reconstruct = {};
% We do the reconstructioon based on the projection for each eigen vector
x = linspace(numel(partialEigen),1,numel(partialEigen));
numVectores = 8;
% Reconstruiremos cada Patch
for patch = 1:width(checker24)
    tmpVector = [];
    % Usando los 4 vectores principales
    for i = 1:numVectores
        if isempty(tmpVector)
            tmpVector = (dot(checker24(:,patch),eV(:,x(i))).*eV(:,x(i)));
        else
            tmp = (dot(checker24(:,patch),eV(:,x(i))).*eV(:,x(i)));
            tmpVector = tmpVector + tmp;
        end
    end
    reconstruct{patch} = tmpVector;
end

%% RCGFC Calculation
GFC = [];

% Iteramos para cada patch (24 iterations)
for i=1:numel(checker24(1,:))
    % Iteramos para cada Lambda (61 iterations)
    Num = 0;
    den1 = 0;
    den2 = 0;
    for L=1:numel(checker24(:,1))
        Num = Num + checker24(L,i)*reconstruct{1,i}(L);
        den1 = den1 + checker24(L,i)^2;
        den2 = den2 + reconstruct{1,i}(L)^2;
    end
    GFC(i) = Num/(sqrt(den1)*sqrt(den2));
    
end
CGFC = 1- GFC;
RCGFC = sqrt(CGFC);



%% Interactive Plot
% create figure and components
fig = uifigure("Name",'GUI Spectral Reconstruction of 24-Color Checker'); 

ax = uiaxes('Parent', fig, 'Position', [10 10 400 400]);

% Create a plot
x = linspace(400,700,61); % Because we have 61 values on measure refl
y.A = checker24(:,1);
y.B = reconstruct{1};
p = plot(ax, x, y.A,x,y.B);
legend(ax, 'Original','Reconstruction');
title(ax,'Original and Reconstruction');
xlabel(ax, 'Wavelength (nm)');
ylabel(ax, 'Reflectance [0-1]')
h = text(ax,450,min(p(1).YData) + 0.5*(max(p(1).YData)-min(p(1).YData)),strcat('RCGFC: ',num2str(RCGFC(1))));

% Create a dropdown component 
dd = uidropdown(fig, 'Position', [430 210 100 22],...
    'Items', {'Patch 1', 'Patch 2','Patch 3','Patch 4','Patch 5', ...
    'Patch 6','Patch 7','Patch 8','Patch 9','Patch 10','Patch 11', ...
    'Patch 12','Patch 13','Patch 14','Patch 15','Patch 16','Patch 17', ...
    'Patch 18','Patch 19','Patch 20','Patch 21','Patch 22','Patch 23', ...
    'Patch 24'},...
    'Value', 'Patch 1',...
    'ValueChangedFcn', @(dd, event) selection(dd,p,checker24,reconstruct,h,RCGFC));
ddVectors = uidropdown(fig, 'Position', [430 290 100 22],...
    'Items', {'2','4','6','8','10','12','18','26','44'},...
    'Value', '4',...
    'ValueChangedFcn', @(ddVectors, event) calcReconst(ddVectors,dd,p,checker24,h,partialEigen,eV));
lbl_eig = uilabel(fig, "Text", {'Select the';'Principal Components'; 'to Reconstruct'},'Position',[430 340 300 42]);
lbl_patch = uilabel(fig, "Text", 'Select Patch','Position',[430 240 100 22]);
% Create ValueChangedFcn callback
function selection(dd, p,checker24,reconstruct,h,RCGFC)
    split_dd_Value = split(dd.Value,' ');
    numeric_dd_value = str2double(split_dd_Value{2});
    val_original = checker24(:,numeric_dd_value);
    val_recons = reconstruct{numeric_dd_value};
    txt_rcgfc = strcat('RCGFC: ',num2str(RCGFC(numeric_dd_value)));

    %Original Data Update
    p(1).YData = val_original;
    % Reconstructed Data Update
    p(2).YData = val_recons;
    h.String = txt_rcgfc;
    h.Position = [450,min(p(1).YData) + 0.5*(max(p(1).YData)-min(p(1).YData))];
end

function calcReconst(ddVectors,dd,p,checker24,h,partialEigen,eV)   
    numVectores = str2double(ddVectors.Value);
    
    reconstruct = {};
    % We do the reconstructioon based on the projection for each eigen vector
    x = linspace(numel(partialEigen),1,numel(partialEigen));
    
    % Reconstruiremos cada Patch
    for patch = 1:width(checker24)
        tmpVector = [];
        % Usando los 4 vectores principales
        for i = 1:numVectores
            if isempty(tmpVector)
                tmpVector = (dot(checker24(:,patch),eV(:,x(i))).*eV(:,x(i)));
            else
                tmp = (dot(checker24(:,patch),eV(:,x(i))).*eV(:,x(i)));
                tmpVector = tmpVector + tmp;
            end
        end
        reconstruct{patch} = tmpVector;
    end
    
    GFC = [];

    % Iteramos para cada patch (24 iterations)
    for i=1:numel(checker24(1,:))
        % Iteramos para cada Lambda (61 iterations)
        Num = 0;
        den1 = 0;
        den2 = 0;
        for L=1:numel(checker24(:,1))
            Num = Num + checker24(L,i)*reconstruct{1,i}(L);
            den1 = den1 + checker24(L,i)^2;
            den2 = den2 + reconstruct{1,i}(L)^2;
        end
        GFC(i) = Num/(sqrt(den1)*sqrt(den2));

    end
    CGFC = 1- GFC;
    RCGFC = abs(sqrt(CGFC));
    split_dd_Value = split(dd.Value,' ');
    numeric_dd_value = str2double(split_dd_Value{2});

    val_original = checker24(:,numeric_dd_value);
    val_recons = reconstruct{numeric_dd_value};
    txt_rcgfc = strcat('RCGFC: ',num2str(RCGFC(numeric_dd_value)));

    %Original Data Update
    p(1).YData = val_original;
    % Reconstructed Data Update
    p(2).YData = val_recons;
    h.String = txt_rcgfc;
    h.Position = [450,min(p(1).YData) + 0.5*(max(p(1).YData)-min(p(1).YData))];
    
end
