%% FunctionDefined.
function UDF = FunctionDefined
UDF.ColorCoordinates1931_SD = @ColorCoordinates1931_SD;
UDF.ColorCoordinates1931_xy = @ColorCoordinates1931_xy;
UDF.ColorCoordinates1960_uv = @ColorCoordinates1960_uv;
UDF.FinalColorCoordinate = @FinalColorCoordinate;
UDF.ColourTemperature = @ColourTemperature;
UDF.ConstructSpectrum = @ConstructSpectrum;
UDF.CameraResponsivity = @CameraResponsivity;
UDF.f = @f;
end

%% ColorCoordinates1931_SD: Calculating the 3 coefficients of monochromatic light.
function [ X, Y, Z, K ] = ColorCoordinates1931_SD( SpectrumData, TristimulusValues, bata )

delta_lambda = 5;  % Interval of wavelength.
num = (780-380)/delta_lambda+1;  % Number of wavelength partitions.

% Calculating Normalization Coefficient K.
Y_K = 0;
for ii=1:1:num
    Y_K = Y_K + SpectrumData( ( ii - 1 ) * delta_lambda + 1, 2 ) * TristimulusValues(ii,2) * delta_lambda;
end
K = 100 / Y_K;

X = 0;
for ii=1:1:num
    X = X + SpectrumData( ( ii - 1 ) * delta_lambda + 1, 2 ) * bata(ii) * TristimulusValues(ii,1) * delta_lambda;
end
X = X * K;

Y = 0;
for ii=1:1:num
    Y = Y + SpectrumData( ( ii - 1 ) * delta_lambda + 1, 2 ) * bata(ii) * TristimulusValues(ii,2) * delta_lambda;
end
Y = Y * K;

Z = 0;
for ii=1:1:num
    Z = Z + SpectrumData( ( ii - 1 ) * delta_lambda + 1, 2 ) * bata(ii) * TristimulusValues(ii,3) * delta_lambda;
end
Z = Z * K;
end

%% ColorCoordinates1931_xy: Calculate color coordinates defined in 1931.
function [ x, y, z ] = ColorCoordinates1931_xy( X, Y, Z )
x = X / ( X + Y + Z );
y = Y / ( X + Y + Z );
z = Z / ( X + Y + Z );
end

%% ColorCoordinates1960_uv: Calculate color coordinates defined in 1964.
function [ u, v ] = ColorCoordinates1960_uv( x, y )
u = 4 * x / ( -2 * x + 12 * y + 3 );
v = 6 * y / ( -2 * x + 12 * y + 3 );
end

%% FinalColorCoordinate: Final ColorCoordinate.
function ColorCoordinate =  FinalColorCoordinate( x, y, u, v )
ColorCoordinate( 1, 1 ) = x; ColorCoordinate( 1, 2 ) = y;
ColorCoordinate( 2, 1 ) = u; ColorCoordinate( 2, 2 ) = v;
end

%% Solve ColourTemperature.
function CCT = ColourTemperature( BlackbodyColorCoordinate_uv, x, y )
CCT = 0;
error = 0.00001;
[ u, v ] = ColorCoordinates1960_uv( x, y );

[ row, ~ ] = size( BlackbodyColorCoordinate_uv );
d_cache = sqrt( ( BlackbodyColorCoordinate_uv( 1, 2 ) )^2 + ( v - BlackbodyColorCoordinate_uv( 1, 3 ) )^2 ); % Is it on track?
if abs( d_cache ) <= error
    CCT = BlackbodyColorCoordinate_uv( 1, 1 );
else
    d_left = d_cache;
end

for i_r = 2 : 1 : row
    d_cache = sqrt( ( u - BlackbodyColorCoordinate_uv( i_r, 2 ) )^2 + ( v - BlackbodyColorCoordinate_uv( i_r, 3 ) )^2 ); % Is it on track?
    if abs( d_cache ) <= error
        CCT = BlackbodyColorCoordinate_uv( i_r, 1 );
        break;
    else
        d_cache = ( BlackbodyColorCoordinate_uv( i_r, 4 ) * u + (-1) * v + ( BlackbodyColorCoordinate_uv( i_r, 3 ) - BlackbodyColorCoordinate_uv( i_r, 4 ) * BlackbodyColorCoordinate_uv( i_r, 2 ) ) ) / ( sqrt( BlackbodyColorCoordinate_uv( i_r, 4 )^2 + 1 ) ); % Is it out of track?
        if d_cache * d_left < 0
            CCT = ( 1 / BlackbodyColorCoordinate_uv( i_r, 1 ) + ( d_cache / ( d_cache - d_left ) ) * ( 1 / BlackbodyColorCoordinate_uv( i_r - 1, 1 ) - 1 / BlackbodyColorCoordinate_uv( i_r, 1 ) ) )^(-1);
            break;
        else
            d_left = d_cache;
        end
    end
end
end

%% Construct Standard Spectrum.
function SpectrumDataSTD = ConstructSpectrum( S_lambda, CCT )
delta_lambda = 1; % Interval of wavelength
num = (780-380)/delta_lambda+1; % Number of wavelength partitions
if CCT <= 3400
    % Spectral energy distribution of Planckian radiation
    for ii=1:1:num
        SpectrumDataSTD(ii,1) = 380+(ii-1)*delta_lambda;
        SpectrumDataSTD(ii,2) = 100*(560/SpectrumDataSTD(ii,1))^5*(exp((1.435e7)/(560*CCT))-1)/(exp((1.435e7)/(SpectrumDataSTD(ii,1)*CCT))-1);
    end
elseif CCT > 3400 && CCT <= 5000
    for ii=1:1:num
        SpectrumDataSTD(ii,1) = 380+(ii-1)*delta_lambda;
    end
    for ii=1:1:num
        P3400(ii) = 100*(560/SpectrumDataSTD(ii,1))^5*(exp((1.435e7)/(560*3400))-1)/(exp((1.435e7)/(SpectrumDataSTD(ii,1)*3400))-1);
    end
    xD = -4.6070*1e9/(5000^3)+2.9678*1e6/(5000^2)+0.09911*1e3/5000+0.244063;
    yD = -3.000*xD^2+2.870*xD-0.275;
    M1 = (-1.34674-1.77861*xD+5.90757*yD)/(0.02387+0.25539*xD-0.73217*yD);
    M2 = (0.03638-31.44464*xD+30.06400*yD)/(0.02387+0.25539*xD-0.73217*yD);
    for ii=1:1:num
        D5000(ii)=S_lambda(ii,1)+M1*S_lambda(ii,2)+M2*S_lambda(ii,3);
    end
    for ii=1:1:num
        SpectrumDataSTD(ii,2) = (D5000(ii)*(CCT-3400)+P3400(ii)*(5000-CCT))/(5000-3400);
    end
else
    if CCT <=7000
        xD = -4.6070*1e9/(CCT^3)+2.9678*1e6/(CCT^2)+0.09911*1e3/CCT+0.244063;
    else
        xD = -2.0064*1e9/(CCT^3)+1.9018*1e6/(CCT^2)+0.24748*1e3/CCT+0.237040;
    end
    yD = -3.000*xD^2+2.870*xD-0.275;
    
    M1 = (-1.34674-1.77861*xD+5.90757*yD)/(0.02387+0.25539*xD-0.73217*yD);
    M2 = (0.03638-31.44464*xD+30.06400*yD)/(0.02387+0.25539*xD-0.73217*yD);
    
    for ii=1:1:num
        SpectrumDataSTD(ii,1) = 380+(ii-1)*delta_lambda;
        SpectrumDataSTD(ii,2) = S_lambda(ii,1)+M1*S_lambda(ii,2)+M2*S_lambda(ii,3);
    end
end
end

%% Camera responsivity in TLCI.
function [ X, Y, Z ] = CameraResponsivity( Rc,Gc,Bc,M,MS,MD )
% Linear matrix.
temp = M*[Rc;Gc;Bc];RM=temp(1);GM=temp(2);BM=temp(3);
temp = MS*[RM;GM;BM];RB=temp(1);GB=temp(2);BB=temp(3);
% Opto-electronic transfer characteristic (gamma correction).
if RB<0.018
    RC_prime = 4.5*RB;
else
    RC_prime = 1.099*RB^(0.45)-0.099;
end
if GB<0.018
    GC_prime = 4.5*GB;
else
    GC_prime = 1.099*GB^(0.45)-0.099;
end
if BB<0.018
    BC_prime = 4.5*BB;
else
    BC_prime = 1.099*BB^(0.45)-0.099;
end
% Electro-optic transfer characteristic (gamma)
% RD = RC_prime^(2.4);GD = GC_prime^(2.4);BD = BC_prime^(2.4);
RD = (RC_prime^12)^(1/5);GD = (GC_prime^12)^(1/5);BD = (BC_prime^12)^(1/5);
% Display primaries
temp=MD*[RD;GD;BD];X=temp(1);Y=temp(2);Z=temp(3);
end

%% f(var).
function fvar = f(var)
if var<(24/116)^3
    fvar = 1/3*(116/24)^2*var+16/116;
else
    fvar = var^(1/3);
end
end
