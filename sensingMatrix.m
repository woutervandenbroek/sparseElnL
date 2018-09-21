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

function A = sensingMatrix( x, M, sensingModel, p )

% x: Object
% M: Number of measurements
% sensingModel 1: Single-pixel with 0/1 pixels
%              2: Single-pixel with fractional pixels
%              3: ADF-STEM CS, p on-pixels
%              4: ADF-STEM, scanning
%              5: Single-pixel scanning with 0/1 pixels
% p: Percentage of on-pixels

if sensingModel == 1
    if nargin < 4
        p = 0.5;
    end
end

N = numel( x );

switch sensingModel
    case { 1, 2, 3 }
        P = round( p * N );
        if ( P < 1 )
            disp('    !!! WARNING in sensingMatrix: p * N < 1 !!! ')
        end
        if ( sensingModel == 1 ) || ( sensingModel == 2 )
            onVal = 1 / M;
        elseif sensingModel == 3
            onVal = 1 / P;
        end
        A = zeros( M, N );
        for j = 1:M
            s = randperm( N );
            if ( sensingModel == 1 ) || ( sensingModel == 3 )
                A( j, s(1:P) ) = onVal;
            end
            if ( sensingModel == 2 )
                A( j, s(1:P) ) = rand( 1, P ) * onVal;
            end
        end
    case 4
        A = eye( N ) * ( M / N );
        A = sparse( A );
    case 5
        A = eye( N ) * ( 1 / N );
        A = sparse( A );
end

% Sparse A can speed up calculations enormously
if p < 0.34
    A = sparse( A );
end