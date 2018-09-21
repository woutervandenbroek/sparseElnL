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

function getMSE_ADFSTEM( fCS )

% Computes the MSE wrt to the ground truth, after registering and rescaling
% to the ground truth, GroundTruth_ADFSTEM.mat
% fCS is a reconstruction from sensing matrix 3

load('GroundTruth_ADFSTEM.mat')
f0 = f0m;
f0 = f0 / mean( f0(:) );

f1 = fCS;
clear fCS;

[n1, n2] = size( f0 );

f1 = reshape( f1, n1, n2 );


ab1 = [ n1*n2 sum( f1(:) ); sum( f1(:) ) sum( f1(:).^2 ) ] ...
    \ [ sum( f0(:) ); sum( f0(:).*f1(:) ) ];
f1 = ab1(1) + ab1(2) * f1;

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 500;
optimizer.MinimumStepLength = 1e-6;
f1_reg = imregister( f1, f0, 'rigid', optimizer, metric );

figure( 1 );
imshowpair( f1_reg, f0, 'diff' );

mask1 = false( n1, n2 );
mask1( 11:(n1-10), 11:(n2-10) ) = true;
ind1 = find( mask1 );
ab1 = [ numel(ind1) sum( f1_reg(ind1) ); sum( f1_reg(ind1) ) sum( f1_reg(ind1).^2 ) ] ...
    \ [ sum( f0(ind1) ); sum( f0(ind1).*f1_reg(ind1) ) ];

f1_reg = ab1(1) + ab1(2) * f1_reg;

mse01 = sum( ( f1_reg( ind1 ) - f0( ind1 ) ).^2 ) / numel( ind1 );

disp( ['Mean squared error, 0 vs 1: ', num2str( mse01 )] )