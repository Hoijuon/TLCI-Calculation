%% Read.

% Read the spectrum.
path = '.\Spectrum\';
FileName = 'SpectrumDataTest';
file = dir( fullfile( path, strcat( FileName, '.xlsx' )));
FullFileName = file.name;
[ ~, sheet, ~] = xlsfinfo( strcat( path, FullFileName ));
SpectrumDataTest = xlsread( strcat( path, FullFileName ), sheet{1} );

% Read tristimulus values of CIE in 1931.
[ data1, data2, data3 ] = textread( 'TristimulusValues1931.txt', '%n%n%n' );
TristimulusValues1931 = [ data1 data2 data3 ];

% Read EBU standard camera curves in 2012.
[ data1, data2, data3 ] = textread( 'EBUStandardCameraCurves2012.txt', '%n%n%n' );
EBUStandardCameraCurves2012 = [ data1 data2 data3 ];

% Read blackbody trajectory in xy and isothermSlope.
[ data1, data2, data3, data4 ] = textread( 'BlackbodyTrajectoryInxyAndIsothermSlope.txt', '%n%n%n%n' );
BlackbodyColorCoordinate_xy  = [ data1 data2 data3 data4 ];

% Read blackbody trajectory in uv and isothermSlope.
[ data1, data2, data3, data4 ] = textread( 'BlackbodyTrajectoryInuvAndIsothermSlope.txt', '%n%n%n%n' );
BlackbodyColorCoordinate_uv = [ data1 data2 data3 data4 ];

% Read Spectral Brightness Coefficient.
[ data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15, data16, data17, data18, data19, data20, data21, data22, data23, data24 ] = textread( 'SpectralBrightnessCoefficient.txt', '%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n%n' );
bata_lambda = [ data1, data2, data3, data4, data5, data6, data7, data8, data9, data10, data11, data12, data13, data14, data15, data16, data17, data18, data19, data20, data21, data22, data23, data24 ];

% Read Standard Spectrum.
[ data1, data2, data3 ] = textread( 'StandardSpectrumLambdaFit.txt', '%n%n%n' );
S_lambda = [ data1 data2 data3 ];