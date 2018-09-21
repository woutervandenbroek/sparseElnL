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

function histogramWidths( f )

N = numel( f );
f = abs( f(:) );

f = sort( f );

iTr = [ 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.98, 0.99 ]';
iTr = max( round( N * iTr ), 1 );

mu = 1 ./ f( iTr );

iTr = [ 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.90, 0.95, 0.98, 0.99 ]';
disp('>>')
disp( '>> Mu values to enclose these fractions of the histogram:' )

for j = 1:numel(iTr)
    disp( ['>> Fraction: ', num2str(iTr( j ), '%3.2f' ), ', Mu-value: ', num2str( mu(j) )] )
end

disp('>>')