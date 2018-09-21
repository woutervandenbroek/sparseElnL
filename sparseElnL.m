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


function [ xR, z, yM, La, L1Sparse, Emeas, Ezx ] = sparseElnL( params )

I0 = params.intensity;
y = params.measurements;
A = params.sensingmatrix;
xR = params.startingguess;
noiseModel = params.noisemodel;
c = params.readoutvariance;
nLineSearch0 = params.nLineSearch0;
nLineSearchIt1 = params.nLineSearchIt1;
N0 = params.N0;
N = params.N1;
n1 = params.n1;
n2 = params.n2;
muLo = params.muLo;
muHi = params.muHi;
beta = params.beta;
c2 = params.c2;
alphaInit = params.alphaInit;
flagTV = params.flagTV;

global nM_kjgahfauzer
nM_kjgahfauzer = noiseModel;

if N < N0
    N0 = N - 1;
end

y = y(:);

La = zeros(N,1);
Emeas = zeros(N,1);
Ezx = zeros(N,1);
L1Sparse = zeros(N,1);

e = c;

filterFctr = 1;
if flagTV == 1
    filterFctr = 2;
end

if flagTV == 1
    z  = zeros( n1 * n2 * filterFctr, 1 ); % z-variable, for shrinkage
    lambda = zeros( n1 * n2 * filterFctr, 1 ); % Langrangian multipliers for agreement between Dx and z
end
if flagTV == 0  % So that there are no wraparound artifacts for Laplacian
    z  = zeros( ( n1 - 2 ) * ( n2 - 2 ) * filterFctr, 1 ); % z-variable, for shrinkage
    lambda = zeros( ( n1 - 2 ) * ( n2 - 2 ) * filterFctr, 1 ); % Langrangian multipliers for agreement between Dx and z
end
nu = 0; % multiplier for agreement with experiments
alpha = alphaInit;
alphaFlag = 0;
alphaArray = alphaInit + zeros( 9, 1 );
mu = 0;
flagLS = 1;

gammaJ = [ relaxParameter( N0/2, N0 ), ones(1,N-N0) ];

DE1 = derivativeLa( xR, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV ); % Derivative of the Augmented Lagrangian

sD = -DE1;

for j = 1:N
    
    if j == 1
        nLineSearch = nLineSearchIt1;
        % mu = 0; Should've already been set above to compute DE1 correctly
    else
        nLineSearch = nLineSearch0;
        mu = muHi * ( 1 - gammaJ(j) ) + muLo * gammaJ(j);
    end
    
    if j > 2 % NOT "j > 1", this is not a typo, the shrinkage for 2 has already been done while j == 1
        % Shrinkage
        z = max( abs( totalVar( xR, n1, n2, flagTV ) - lambda / mu ) - 1 / mu, 0 ) .* sign( totalVar( xR, n1, n2, flagTV ) - lambda / mu );
    end
    
    for k = 1:nLineSearch
        
        [sD, DE1, alpha] = searchDirection( sD, DE1, alpha, xR, lambda, nu, z, A, mu, beta, y, e, n1, n2, flagLS, I0, flagTV );

        if ( k == 1 )
            DE1Norm  = sqrt( sum( DE1.^2 ) );
        else
            DE1Norm = max( DE1Norm, sqrt( sum( DE1.^2 ) ) );
        end
     
        [xR, alpha, E1, DE2, flagLS] = lineSearchRobust( sD, DE1, xR, alpha, c2, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
        disp(['Step size: ', num2str( alpha ), ', LINE SEARCH FLAG: ', int2str(flagLS) ])
        DE2Norm = sqrt( sum( DE2.^2 ) );
                
        if ( flagLS == 1 ) || ( flagLS == 2 ) % Fill up an array with typical step sizes
            alphaArray( 2:9 ) = alphaArray( 1:8 );
            alphaArray( 1 ) = alpha;
        end
        if ( flagLS == 3 ) || ( flagLS == 0 )
            alpha = median( alphaArray ) * 1e-1 ;
            alphaFlag = alphaFlag + 1;
        else
            alphaFlag = max( 0, alphaFlag - 1 );
        end
        
        if ( j == 1 ) && ( ( DE2Norm / ( n1*n2 ) ) < 1e-9 )
            disp(' ')
            disp( ['>>>> Stopping line search in first epoch at iteration ', int2str( k )] )
            break
        end
                
    end
    disp(' ')
    disp(['    Relative norm of the derivative: ', num2str( DE2Norm / DE1Norm )])
    disp(['    Absolute norm of the derivative: ', num2str( DE2Norm / (n1*n2) )])
    
    if j < 2
        % Set mu to muHi for j == 1. The other js are taken care of above
        mu = muHi;
        % Shrinkage
        z = max( abs( totalVar( xR, n1, n2, flagTV ) - lambda / mu ) - 1 / mu, 0 ) ...
            .* sign( totalVar( xR, n1, n2, flagTV ) - lambda / mu );
    end
    
    % Update of the multipliers
    lambda = lambda - mu * ( totalVar( xR, n1, n2, flagTV ) - z );
    nu = nu - beta * MeasurementConstraint( A, xR, y, e, I0 );
    disp( [ '    Multiplier Nu: ', num2str(nu) ] )
    disp(' ')
    
    % Record the augmented Lagrangian
    La( j ) = E1;
    L1Sparse( j ) = mean( abs( z ) );
    Emeas( j ) = beta * MeasurementConstraint( A, xR, y, e, I0 )^2;
    Ezx( j ) = mu * sum( ( totalVar( xR, n1, n2, flagTV ) - z ).^2 );
    
    % Calculate the simulated projections
    yM = A * xR;
    
    % Save the results
    saveDisplayResults( j, xR, z, yM, La, Emeas, Ezx, L1Sparse );
    
    if alphaFlag > 11
        disp('Step size was below 5e-16 too often, exiting program. ')
        break
    end
end

% N = 1 means a test run to get realistic values for muLo and muHi
if N == 1
    disp('>> It is recommended to chose muHi at around the 0.1 or 0.2 value')
    disp('>> and muLo at around the 0.8 or 0.9 value.')
    disp('>> ')
    tmp = totalVar( xR, n1, n2, flagTV );
    histogramWidths( tmp );
    clear tmp
end

xR = reshape( xR, n1, n2 );
if flagTV == 0
    z = reshape( z, n1 - 2, n2 - 2 );
    z = [ zeros( 1, n2 - 2 ); z; zeros( 1, n2 - 2 ) ];
    z = [ zeros( n1, 1 ), z, zeros( n1, 1 ) ];
end
z = reshape( z, filterFctr * n1, n2 );


function [x, alpha1, E1, DE, flag] = lineSearchRobust( d, DE, x, alpha, c2, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV )

% Line search like in FDES

c1 = 1e-4;

alpha = abs( alpha );
alphaLo = 0.3162 * alpha;
flag = 0;

dO = d' * DE;

alpha0 = 0;
E0 = partialError( x, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
EO = E0;
d0 = dO;
alpha1 = alpha;

% First test if strong Wolfe condition is fullfilled immediately
E1 = partialError( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
DE = derivativeLa( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
d1 = d' * DE;
if ( ( E1 < ( EO + c1 * alpha1 * dO ) ) && ( abs( d1 ) < c2 * abs( dO ) ) )
    flag = 1;
    x = x + alpha1 * d;
    return;
end

% Find an alpha1 with positive derivative, then start line search
for i = 1:10
    if d1 > 0
        if i ~= 1
            E1 = partialError( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
        end
        break
    end
    alpha1 = alpha1 + i * alpha;
    DE = derivativeLa( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
    d1 = d' * DE;
    if i == 10
        E1 = partialError( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
    end
end

alphaHi = 2 * alpha1;
alpha2 = alpha1;
E2 = E1;
d2 = d1;

for i = 1:20
    % First ten times we try the strong wolfe conditions
    % disp(['        Internal iteration: ', int2str(i)])
    if i < 11
        % Test for strong Wolfe conditions
        if ( ( E2 < ( EO + c1 * alpha2 * dO ) ) && ( abs( d2 ) < c2 * abs( dO ) ) )
            flag = 1;
            E1 = E2;
            alpha1 = alpha2;
            break;
        end
        
        % Find the position of the minimum of the cubic function fitted to E0, d0 and E1, d1:
        alpha2 = minimumThirdOrderPoly( alpha0, alpha1, E0, E1, d0, d1 );
        
        % Apply safeguards:
        if ( alpha2 > alphaHi )
            alpha2 = alphaHi;
            % Increase alphaHi Linearly
            alphaHi = alphaHi + i * alpha;
        end
        if ( alpha2 < alphaLo )
            alpha2 = alphaLo;            
            alphaLo = 0.3162 * alphaLo;
        end
        
        % get value of error at x + alpha2*d
        E2 = partialError( x + alpha2*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
        % derivative along d
        DE = derivativeLa( x + alpha2*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
        d2 = d' * DE;
        
        % Keep a positive and a negative derivative, i.e. keep the minimum
        % bracketed
        if ( d2 < 0 )
            E0 = E2;
            d0 = d2;
            alpha0 = alpha2;
        else
            E1 = E2;
            d1 = d2;
            alpha1 = alpha2;
        end
        
        %Stop if alpha2 is too small
        if ( alpha2 < 5e-16 )
            flag = 3;
            alpha1 = 0;
            E1 = EO;
            DE = derivativeLa( x, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
            break; 
        end
         
    % Now try the weak wolfe conditions, continue with alpha1
    else
        c1 = 0.20;
        E1 = partialError( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
        if ( E1 < ( EO + c1 * alpha1 * dO )  )
            flag = 2;
            DE = derivativeLa( x + alpha1*d, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
            break;
        end
        % Decrease step size with factor of sqrt( 0.1 ) until weak wolfe condition is met. 
        alpha1 = alpha1 * 0.3162;
        %Stop if alpha1 is too small
        if ( alpha1 < 5e-16 )
            flag = 3;
            alpha1 = 0;
            E1 = EO;
            DE = derivativeLa( x, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
            break; 
        end
    end
end

if ( flag == 0 ) % Means neither Weak nor Strong Wolfe conditions were met. The whole optimization is stopped, because something must be wrong.
    alpha1 = 0;
    E1 = EO;
    DE = derivativeLa( x, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );
end

x = x + alpha1 * d;


function xMin = minimumThirdOrderPoly( x0, x1, f0, f1, dF0, dF1 )

% Find the coefficients of the polynomial a * x^3 + b * x^2 + c * x + d = f
% but d is not needed
c      = dF0 * ( x1 - x0 );
dF1x01 = dF1 * ( x1 - x0 );
a =  2.0 * f0 - 2.0 * f1 +       c + dF1x01;
b = -3.0 * f0 + 3.0 * f1 - 2.0 * c - dF1x01;

% Now determine the position of the minimum
Det = b * b - 3.0 * a * c;
if ( Det > 0.0 )
    Det = sqrt( Det );
    xMin = ( -b + Det ) / ( 3.0 * a );
    if( ( 6.0 * a * xMin + 2.0 * b ) < 0.0 )
        xMin = ( -b - Det ) / ( 3.0 * a );
    end
    xMin = xMin * ( x1 - x0 );
    xMin = xMin + x0;
elseif ( ( ( f1 - f0 ) * ( x1 - x0 ) ) > 0.0 ) % Should be a division "/", but only the sign is needed, so a multiplication gives that as well, but without the risk of dividing by zero
    xMin = -5e16;
else
    xMin =  5e16;
end

% NOTE: This is the expression for a quadratic function:
% x1 = x1 - x0;
% x1 =  x0 + 0.5f * dF0 * x1 * x1 / ( f0 - f1 + x1 * dF0 + 1e-8f );
% This does not use dF1, only the value and the derivative in x0 and the value in x1


function [sD, DE2, alpha] = searchDirection( sD, DE1, alpha0, x, lambda, nu, z, A, mu, beta, y, e, n1, n2, flag, I0, flagTV )
% Polak-Ribiere search direction for conjugate gradients

DE2 = derivativeLa( x, lambda, nu, z, A, mu, beta, y, e, n1, n2, I0, flagTV );

if ( flag == 3 ) || ( flag == 0 )
    betaPR = 0;
else
   % Polak Ribiere: Has automatic reset and a soft reset
    betaPR = ( DE2' * ( DE2 - DE1 ) ) / ( DE1' * DE1 );
    betaPR = max( 0, betaPR );
end

alpha = alpha0 * abs( sD' * DE1 );
sD = -DE2 + betaPR * sD;
alpha = alpha / abs( sD' * DE2 ); % See p. 58 of Nocedal & Wright

if -( sD' * DE2 ) < 0.0  % Restart for PR-CG, sD is no descent direction 
    sD = -DE2;
    alpha = alpha0;
end


function saveDisplayResults( j, xR, z, yM, La, Emeas, Ezx, L1Sparse )

k = j;
if j == 1
    k = 0;
end

if ( mod( k, 1 ) == 0 ) || ( isPowOfHalf( k ) ) % play around with condition if you only want to save at powers of two.
    save xR xR
    save z z
    save yM yM
    save La La
    save Emeas Emeas
    save Ezx Ezx
    save L1Sparse L1Sparse
end


function La = partialError( x, lambda, nu, z, A, mu, beta, y, c, n1, n2, I0, flagTV )

temp = MeasurementConstraint( A, x, y, c, I0 );

La = sum( abs( z ) );

La = La - (lambda') * ( totalVar( x, n1 ,n2, flagTV ) - z ) ... 
    + 0.5 * mu * sum( ( totalVar( x, n1, n2, flagTV ) - z ).^2 ) ...
    - nu * temp + 0.5 * beta * temp.^2;


function dLa = derivativeLa( x, lambda, nu, z, A, mu, beta, y, c, n1, n2, I0, flagTV )

global nM_kjgahfauzer
nM = nM_kjgahfauzer;

temp = A * x;
if ( nM == 1 )
    e = max( 1e-1 / I0, 1e-9 );
    temp = A' * ( ( y + 0.5 / I0 ) .* dsoftlog( temp, e ) - 1 );
end
if ( nM == 2 )
    e = max( c*1e-1, 1e-6 );
    temp = ( temp - y ) ./ ( softplus( temp, e ) + c );
    temp = A' * ( temp .* ( 2 - temp .* dsoftplus( A * x, e ) ) );
end
if ( nM == 3 )
    temp = ( temp - y ) ./ ( c );
    temp = A' * ( 2 * temp );
end
dLa = temp * ( -nu + beta * MeasurementConstraint( A, x, y, c, I0 ) );
temp = lambda + mu * ( z - totalVar( x, n1, n2, flagTV ) );
temp = totalVarTranspose( temp, n1, n2, flagTV );
dLa = dLa - temp;


function  E = MeasurementConstraint( A, x, y, c, I0 ) 

global nM_kjgahfauzer
nM = nM_kjgahfauzer;

if ( nM == 1 )
    e = max( 1e-1 / I0, 1e-9 );
    tmp = A * x;
    E = sum( ( y + 0.5 / I0 ) .* ( 1 - softlog( y, e ) + softlog( tmp, e ) ) - tmp );
end
if ( nM == 2 )
    e = max( c*1e-1, 1e-6 );
    E = sum( ( A * x - y ).^2 ./ ( softplus( A * x, e ) + c ) - 1 / I0 );
end
if ( nM == 3 )
    E = sum( ( A * x - y ).^2 ./ ( c ) - 1 / I0 );
end


function iHalf = isPowOfHalf( j )

iHalf = false;

if isPowOfTwo( j )
    iHalf = true;
else
    jTest = round( 2^( floor( log2( j ) ) + 0.5 ) );
    if ( j == jTest )
        iHalf = true;
    end
end


function iHalf = isPowOfTwo( j )

iHalf = false;

jTest = round( 2^( floor( log2( j ) ) ) );
if ( j == jTest )
    iHalf = true;
end


function z = totalVar( x, n1, n2, flag )

if flag == 1
    z = zeros( 2 * n1, n2 );
    
    ker = zeros( n1, n2 );
    ker( [1 2], 1 ) = [ -1; 1];  % MUST match TotalVarTranspose
    z( 1:n1, : ) = real( ifft2( fft2( reshape( x, n1, n2 ) ) .* fft2( ker ) ) );
    
    ker = zeros( n1, n2 );
    ker( 1, [1 2] ) = [ -1, 1];  % MUST match TotalVarTranspose
    z( (n1 + 1):(2*n1), : ) = real( ifft2( fft2( reshape( x, n1, n2 ) ) .* fft2( ker ) ) );
    
    z = z(:);
else
    ker = zeros( n1, n2 );
    ker( 1:2, 1:2 ) = [ -4 1; 1 0 ] * 0.2236;  % MUST match TotalVarTranspose
    ker( 1, n2  ) = 1 * 0.2236;
    ker( n1, 1 ) = 1 * 0.2236;
    z = real( ifft2( fft2( reshape( x, n1, n2 ) ) .* fft2( ker ) ) );
    z = z( 2:( n1-1 ), 2:( n2-1 ) );
    
    z = z(:);
end
 

function x = totalVarTranspose( z, n1, n2, flag )

if flag == 1
    z = reshape( z, 2 * n1, n2 );
    
    ker = zeros( n1, n2 );
    ker( [1 2], 1 ) = [ -1; 1];  % MUST match TotalVar
    ker = rotConvKernel( ker );
    x = real( ifft2( fft2( z( 1:n1, : ) ) .* fft2( ker ) ) );
    
    ker = zeros( n1, n2 );
    ker( 1, [1 2] ) = [ -1, 1];  % MUST match TotalVar
    ker = rotConvKernel( ker );
    x = x + real( ifft2( fft2( z( (n1+1):(2*n1), : ) ) .* fft2( ker ) ) );
    
    x = x(:);
else
    z = reshape( z, n1 - 2, n2 - 2 );
    z = [ zeros( 1, n2 - 2 ); z; zeros( 1, n2 - 2 ) ];
    z = [ zeros( n1, 1 ), z, zeros( n1, 1 ) ];
    
    ker = zeros( n1, n2 );
    ker( 1:2, 1:2 ) = [ -4 1; 1 0 ] * 0.2236;  % MUST match TotalVar
    ker( 1, n2  ) = 1 * 0.2236;
    ker( n1, 1 ) = 1 * 0.2236;
    ker = rotConvKernel( ker );
    
    x = real( ifft2( fft2( z ) .* fft2( ker ) ) );
    x = x(:);
end


function ker = rotConvKernel( ker )

[ n1, n2 ] = size( ker );

ker = rot90( fftshift( ker ), 2 );
if mod( n1, 2 ) == 0
    ker = circshift( ker, 1, 1 );
end
if mod( n2, 2 ) == 0
    ker = circshift( ker, 1, 2 );
end
ker = ifftshift( ker );


function y = softplus( x, e )

y = x;

t = ( x / e ) < 4.5;

y( t ) = e * log( 1 + exp( x( t ) / e ) );


function dy = dsoftplus( x, e )

dy = ones( size( x ) );

t = ( x / e ) < 4.5;

dy( t ) = ( 1 + exp( -x( t ) / e ) ).^(-1);


function y = softlog( x, e )

x = x / e;
y = zeros( size( x ) );

tr = 5.5;

t = x < ( -tr + 1e-2 );
y( t ) = x( t ) - 0.5 * exp( x( t ) ) + log( e );

t = x > ( tr - 1e-2 );
y( t ) = log( x( t ) + exp( -x( t ) ) ) + log( e );

t = ( ( x > -tr ) & ( x < tr ) );
y( t ) = log( e * log( 1 + exp( x( t ) ) ) );


function dy = dsoftlog( x, e )

x = x / e;
dy = zeros( size( x ) );

tr = 5.5;

t = x < ( -tr + 1e-2 );
dy( t ) = 1/e - ( 0.5 / e ) * exp( x( t ) );

t = x > ( tr - 1e-2 );
dy( t ) = ( 1/e - exp( -x( t ) ) / e ) ./ ( x( t ) + exp( -x( t ) ) );

t = ( ( x > -tr ) & ( x < tr ) );
dy( t ) = exp( x( t ) ) ./ ( e * log( 1 + exp( x( t ) ) ) .* ( 1 + exp( x( t ) ) ) );



function gamma = relaxParameter( n, N )

% Polynomial with derivative k in the origin and first and second
% derivative zero in 1
k = N / n;
Ainv = [ 6 -3; -8 5; 3 -2 ];
a234 = Ainv * ([1-k -k]');
gamma = (0:(N-1)) / (N-1);
gamma = k*gamma + a234(1) * gamma.^2 + a234(2)*gamma.^3 + a234(3)*gamma.^4;
gamma = gamma * ( 1 - 1e-6 ) + 1e-6;