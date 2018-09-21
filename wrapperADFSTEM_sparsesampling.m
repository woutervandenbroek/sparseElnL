% Copyright 2018 Wouter Van den Broek
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


% Wrapper for ADFSTEM

% Input: See setDefaultParameters

% Output
% xR: Raw reconstructed object
% xS: Smoothed reconstructed object
% z: The sparse help variable
% yM: The simulated measurements. 
% La: error function (Augmented Lagrangian)
% L1: Absolute value of z.
% Emeas: Measurement constraint, similarity between model and experiments
% Ezx: Similarity between z and x

clear

% Set noise model to 2, i.e. Poisson and read-out, even though the read-out
% is reall small.
noiseModel = 2;
% Flag for regularization. 1 for TV, 0 for Laplacian
flagTV = 0;
% Problem specific parameters, i.e. of the measurements and the object
n1 = 256;
n2 = 256;

% Optimization parameters for ADF denoising. See setDefaultParameters for
% an explanation
% These parameters are also given in the file parameters_ADFSTEM_sM4_nM2
nLineSearch0 = 10;
nLineSearchIt1 = 1000;
N0 = 20;
N1 = 200;
muLo = 10;
muHi = 150;
beta = 1e-5 * muHi;
c2 = 0.01;
alphaInit = 1e-12;
% Not important in this example as the sensing matrix A is made "by hand" 
% below:
sensingModel = 3;
% Impinging electron dose
I0 = 107.61;
% Load the noisy recording
load('yADF-DN.mat')
y = y(:);

% Set up for inpainting with 10%, 20%, ..., 100% of pixels
M0 = round( ( 0.1:0.1:1.0 )' * n1 * n2 );
nM = numel( M0 );
% Read-out noise variance
c = 0.18160;

if noiseModel == 1
    y( y < ( 2 * sqrt( c ) ) ) = 0;
end
y = y / I0;
c = c / I0;

% Array for all reconstructions
xR0 = zeros( n1, n2, nM );
y0 = y;
for j = 1:nM
    disp( ' ' )
    disp( ['>>  >>  >>  >>  ITERATION NO. ', int2str( j ) '  <<  <<  <<  <<'] )
    disp( ' ' )
    
    % Choose M0(j) pixels at random from y0
    M = M0( j );
    tmp = randperm( n1 * n2 );
    msk = tmp( 1:M );
    clear tmp
    msk = sort( msk );
    A = sparse( 1:M, msk, ones( M, 1 ), M, n1*n2 );
    y = y0( msk );
    
    % Save the measurements
    save y y

    % Starting guess
    xR = A' * y(:);
    xR = xR + 1e-2 * std( xR(:) ) * randn( n1 * n2, 1 );
    xR = xR / ( mean( A(:) )^2 * n1 * n2 * M );
    
    parameters = setParameters( I0, y, A, xR, sensingModel, noiseModel, c, ... 
        nLineSearch0, nLineSearchIt1, N0, N1, n1, n2, ...
        muLo, muHi, beta, c2, alphaInit, flagTV  );
    
    [ xR, z, yM, La, L1, Emeas, Ezx ] = sparseElnL( parameters );
    
    % Uncomment to save your results
    % save(['sparseADFSTEM_', int2str(j), '.mat'], ... 
    %         'A', 'y0', 'yM', 'y', 'xR');
    
    % Amass all reconstructions
    xR0( :, :, j ) = xR;
end