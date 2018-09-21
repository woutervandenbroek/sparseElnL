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

function parameters = setParameters( I0, y, A, xR, sensingModel, ...
    noiseModel, c,nLineSearch0, nLineSearchIt1, N0, N1, n1, n2, ...
    muLo, muHi, beta, c2, alphaInit, flagTV )

% Problem specific parameters, i.e. of the measurements and the object
parameters.n1 = n1;
parameters.n2 = n2;
parameters.intensity = I0;
parameters.measurements = y(:);
parameters.sensingmatrix = A;
parameters.startingguess = xR(:);
parameters.sensingmodel = sensingModel;
parameters.noisemodel = noiseModel;
parameters.readoutvariance = c;
parameters.flagTV = flagTV;

% Optimization parameters
parameters.nLineSearch0 = nLineSearch0;
parameters.nLineSearchIt1 = nLineSearchIt1;
parameters.N0 = N0;
parameters.N1 = N1;
parameters.muLo = muLo;
parameters.muHi = muHi;
parameters.beta = beta;
parameters.c2 = c2;
parameters.alphaInit = alphaInit;