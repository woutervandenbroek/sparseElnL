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

function visualizeResults

close all

figure(1);
load('xR.mat')
n = round( sqrt( numel(xR) ) );
imagesc( reshape( xR, n, n ))  ; colormap gray, axis equal, axis tight, colorbar

figure(2);
load('La.mat')
load('L1Sparse.mat')
load('Emeas.mat')
load('Ezx.mat')
subplot(2,2,1); plot( La )
subplot(2,2,2); plot( 1:numel(L1Sparse),L1Sparse,'r' )
title( 'L1. Blue: Atomic. Red: Sparse' )
subplot(2,2,3); plot( Emeas )
title( 'Emeas' )
subplot(2,2,4); plot( Ezx )
title( 'Ezdx' )

figure(3);
load('y.mat')
load('yM.mat')
plot(1:numel(y),y,'b.', 1:numel(yM),yM,'r')
title('Blue: Original, Red: Modeled')