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

function [ c, gam ] = gam2C( x, I0, gam, sensingModel, noiseModel )

% sensingModel 1: Single-pixel with 0/1 pixels
%              2: Single-pixel with fractional pixels
%              3: ADF-STEM CS, p on-pixels
%              4: ADF-STEM, scanning

switch sensingModel
    case { 1, 2, 5 }
        N = numel( x );
        c = gam * I0 * mean( x(:) ) / N;
    case { 3, 4 }
        c = gam * I0 * mean( x(:) );
end

if noiseModel == 1
    c = 0;
    gam = 0;
end