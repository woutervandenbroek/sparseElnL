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

% Input: See setDefaultParameters

% Output
% xR: Raw reconstructed object
% z: The sparse help variable
% yM: The simulated measurements. 
% La: error function (Augmented Lagrangian)
% L1: Absolute value of z.
% Emeas: Measurement constraint, similarity between model and experiments
% Ezx: Similarity between z and x

clear

% See parameter-files parameters_SINPIX_sMx_nMy for examples:
% sMx:  sensing matrix x, with x = 1, 2, 3 or 4; see below
% nMy:  noise model y, with y = 1, 2 or 3; see below

% Number of Monte Carlo iterations
nMC = 10;

% Flag for TV (1) or Laplacian (0) regularization
flagTV = 0;

% Size of the n x n test object
n = 100;
N = n * n;
% Number of measurements
M = 2000;
% Set the fraction of on-pixels
% p = linspace( 0.01, 0.9, 10 ); 
p  = exp( linspace( log( 0.01 ), log(0.50), 3 ) );

% What sensing matrix to use
% 1: Single-pixel with 0/1 pixels
% 2: Single-pixel with fractional pixels
% 3: ADF-STEM CS
% 4: ADF-STEM, scanning
sensingModel = 1;
% Noise model
% 1: Pure Poisson noise
% 2: Poisson + read-out noise
% 3: only read out noise
noiseModel = 2;
% Electron/photon dose y = I0 A x
I0 = 7.460e5;
% Relative read-out noise
gam = 625.0; %Set to 0 for no read-out

% Make the phantom
if flagTV == 1  % Piecewise constant
    x = phantom( n );
    x( x < 0 ) = 0; % Shouldn't be necessary...
    % Add small background to mimic real experiments better
    x = x * 0.9 + 0.1; % 0.1 < x < 1
end
% Or load one if you need piecewise linear
if flagTV == 0  % Piecewise linear
    load('rampDiscs100x100.mat')
    n = 100;
    N = n * n;
    M = 2000;
end

% Optimization parameters. See setDefaultParameters for an explanation.
% Number of line searches.
nLineSearch0 = 10;  
% Number of line searches for Step 1. of the first iteration. This number
% is not supposed to be reached, convergence should've set in before the 
% end. Increase if not!
nLineSearchIt1 = 1000;
% The number of iterations needed to bring mu down from muHi to muLo
N0 = 20;
% The total number of iterations. Set to 1 for the first trial to get a
% estimates of good muHi and muLo values
N1 = 50;
muLo = 10;
muHi = 100;
beta = 1e5 * muHi;
% Convergence criterion of the line search. Strong Wolfe condition. Best
% set to a small value; 0.01 is recommended.
c2 = 0.01;
% Initial step size for the line search. If the algorithm fails, this is
% the first parameter to change
alphaInit = 1e-12;

% Mean squared error, initialize as zeros
mse = zeros( numel(p), nMC );
% Make a mask so that the 2-pixel wide edge is not used for the MSE calcs
msk = zeros( n, n );
msk( 3:(n-2), 3:(n-2) ) = 1;
msk = find( msk );

for jMC = 1:nMC
    for jP = 1:numel( p )
        disp(' ')
        disp(['>>  >>  >>  >>  MONTE CARLO ITERATION ', int2str(jMC), ' and ' int2str(jP), '  <<  <<  <<  <<'])
        disp(' ')
        
        % Build the sensing matrix based on the sensing model
        A = sensingMatrix( x, M, sensingModel, p( jP ) );
        
        % Convert the relative variance gam to an absolute variance c
        [c, gam] = gam2C( x, I0, gam, sensingModel, noiseModel );
        
        % Build the observations y:
        y = I0 * A * x(:);
        
        % add noise
        y = addNoise( y, noiseModel, c );
        
        % Rescale the problem for better optimization:
        y = y / I0;
        c = c / I0;
        
        % Save the measurements
        save y y
        
        % Problem specific parameters, i.e. of the measurements and the object
        n1 = n;
        n2 = n;
        
        % Starting guess
        xR = A' * y(:);
        xR = xR + 1e-2 * std( xR(:) ) * randn( n1 * n2, 1 );
        xR = xR / ( mean( A(:) )^2 * N * M );
        
        parameters = setParameters(  I0, y, A, xR, sensingModel, noiseModel, c, ... 
            nLineSearch0, nLineSearchIt1, N0, N1, n1, n2, ...
            muLo, muHi, beta, c2, alphaInit, flagTV  );
        
        [ xR, z, yM, La, L1, Emeas, Ezx ] = sparseElnL( parameters );
        
        mse( jP, jMC ) = mean( ( xR( msk ) - x( msk ) ).^2 );
        
        % Uncomment if you want to save your results
        % save(['myResults_date_', int2str(jP), '-', int2str(jMC), '.mat'], ... 
        %    'mse', 'p', 'parameters', 'xR');
        
        % Display the MSE
        disp( ['MSE: ', num2str(mse( jP, jMC ))] )
    end
end
