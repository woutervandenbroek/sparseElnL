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

function parameters = setDefaultParameters

% Problem specific parameters, i.e. of the measurements and the object
parameters.n1 = 100;             % First dimension of the solution
parameters.n2 = 100;             % Second dimension of the solution
parameters.intensity = 5e5;      % Number of photons/electrons: y = I0 * A * x
parameters.measurements = 0;     % The measurements y: y = I0 * A * x
parameters.sensingmatrix = 0;    % The sensing matrix A: y = I0 * A * x
parameters.startingguess = 0;    % The starting guess for the unkown object x: y = I0 * A * x
parameters.sensingmodel = 1;     % The sensing model 1: Single-pixel with 0/1 pixels
                                 %                   2: Single-pixel with fractional pixels
                                 %                   3: ADF-STEM CS, p on-pixels
                                 %                   4: ADF-STEM, scanning
parameters.noisemodel = 2;       % Noise model 1: only Poisson noise
                                 %             2: Poisson and read-out noise
                                 %             3: only read-out noise
parameters.readoutvariance = 1;  % Variance of the read-out noise
parameters.flagTV = 1;           % 1 for TV, 0 for Laplacian regularization

% Optimization parameters
parameters.nLineSearch0 = 20;       % No. of iterations of the line searches (except the first one)
parameters.nLineSearchIt1 = 10000;  % No. of iterations for the FIRST line search
parameters.N0 = 50;                 % Number of iterations to drive mu from muHi to muLo
parameters.N1 = 200;                % Total number of iterations
parameters.muLo = 1;                % mu is driven from muHi to muLo. 
                                    % mu controls agreement between x and z
parameters.muHi = 10;               % Typically > 10 * muLo
parameters.beta = 1;                % Controls agreement with measurements y. 
                                    % Typically equal to muLo or muHi
parameters.c2 = 0.1;                % Wolfe condition on the derivative.
                                    % Lower seems better.
parameters.alphaInit = 1e-6;        % Initial step for line search.