%% main program to calculate TLCI.
clear;clc;
format long;

% Variable.
% bata = ones(81,1);
S = 90; % Saturation level S.
Kn = 1;

% Constants.
FunctionConstants;

% Read.
FunctionRead;

% FunctionCalculate.
UDF = FunctionDefined;

% Calculating the 3 coefficients (XYZ) of monochromatic light.
[Xt,Yt,Zt,Kt] = UDF.ColorCoordinates1931_SD( SpectrumDataTest, TristimulusValues1931, ones(81,1) );
% Calculate xy color coordinates defined in 1931.
[xt,yt,zt] = UDF.ColorCoordinates1931_xy( Xt, Yt, Zt );
% Calculate color coordinates defined in 1964.
[ut,vt] = UDF.ColorCoordinates1960_uv( xt, yt );
% Final ColorCoordinate.
ColorCoordinate =  UDF.FinalColorCoordinate( xt, yt, ut, vt );

% ColourTemperature.
CCT = round( UDF.ColourTemperature( BlackbodyColorCoordinate_uv, xt, yt ) );
if CCT < 2000 || CCT >25000
    disp( 'Warning: Colour Temperature!' );
end

% Construct Standard Spectrum.
SpectrumDataSTD = UDF.ConstructSpectrum( S_lambda, CCT );
% Construct white point D65.
SpectrumDataW = UDF.ConstructSpectrum( S_lambda, 6500 );

% The tristimulus values of the display white point.
% [XW,YW,ZW,~] = UDF.ColorCoordinates1931_SD( SpectrumDataW, EBUStandardCameraCurves2012, ones(81,1) );
[XW,YW,ZW,~] = UDF.ColorCoordinates1931_SD( SpectrumDataW, TristimulusValues1931, ones(81,1) );

for Ri=1:1:24 % R1 - R24
    bata = bata_lambda(:,Ri);
    % Camera responsivity in TLCI.
    [Rc,Gc,Bc,~] = UDF.ColorCoordinates1931_SD( SpectrumDataTest, EBUStandardCameraCurves2012, bata ); % Responsivity curves
    % colour balancing.
    [Rb,Gb,Bb,~] = UDF.ColorCoordinates1931_SD( SpectrumDataTest, EBUStandardCameraCurves2012, ones(81,1) );
    temp=Gb;
    Rb=Rb*Kn/temp;Gb=Gb*Kn/temp;Bb=Bb*Kn/temp;
    Rcb=Rc*1/Rb;Gcb=Gc*1/Gb;Bcb=Bc*1/Bb;
    [Xti,Yti,Zti] = UDF.CameraResponsivity( Rcb,Gcb,Bcb,M,MS,MD );
    % Camera responsivity in TLCI.
    [Rc,Gc,Bc,~] = UDF.ColorCoordinates1931_SD( SpectrumDataSTD, EBUStandardCameraCurves2012, bata );
    % colour balancing.
    [Rb,Gb,Bb,~] = UDF.ColorCoordinates1931_SD( SpectrumDataSTD, EBUStandardCameraCurves2012, ones(81,1) );
    temp=Gb;
    Rb=Rb*Kn/temp;Gb=Gb*Kn/temp;Bb=Bb*Kn/temp;
    Rcb=Rc*1/Rb;Gcb=Gc*1/Gb;Bcb=Bc*1/Bb;
    [Xri,Yri,Zri] = UDF.CameraResponsivity( Rcb,Gcb,Bcb,M,MS,MD );
    
    % calculate CIELAB values for the Test- and Reference-illuminated colour samples.
    Lti_asterisk = 116*UDF.f(Yti/YW)-16; Lri_asterisk = 116*UDF.f(Yri/YW)-16;
    ati_asterisk = 500*(UDF.f(Xti/XW)-UDF.f(Yti/YW)); ari_asterisk = 500*(UDF.f(Xri/XW)-UDF.f(Yri/YW));
    bti_asterisk = 200*(UDF.f(Yti/YW)-UDF.f(Zti/ZW)); bri_asterisk = 200*(UDF.f(Yri/YW)-UDF.f(Zri/ZW));
    
    % calculate intermediate values for the manipulation of the CIELAB results.
    Cti_asterisk = sqrt(ati_asterisk^2+bti_asterisk^2); Cri_asterisk = sqrt(ari_asterisk^2+bri_asterisk^2);
    Ci_asterisk_ave = (Cti_asterisk+Cri_asterisk)/2;
    gi = 0.5*(1-sqrt(Ci_asterisk_ave^7/(Ci_asterisk_ave^7+25^7)));
    ati_prime = (1+gi)*ati_asterisk; ari_prime = (1+gi)*ari_asterisk;
    Cti_prime = sqrt(ati_prime^2+bti_asterisk^2); Cri_prime = sqrt(ari_prime^2+bri_asterisk^2);
    % The hue, h, is measured in degrees, and uses a four-quadrant arctangent.
    if bti_asterisk==0 && ati_prime==0
        hti = 0;
    elseif ati_prime==0 && bti_asterisk<0
        hti = -90;
    elseif ati_prime==0 && bti_asterisk>0
        hti = 90;
    elseif ati_prime<0 && bti_asterisk<0
        hti = (atan(bti_asterisk/ati_prime)-pi)/pi*180;
    elseif ati_prime<0 && bti_asterisk>=0
        hti = (atan(bti_asterisk/ati_prime)+pi)/pi*180;
    else
        hti = atan(bti_asterisk/ati_prime)/pi*180;
    end
    if hti>360
        hti = hti - floor(hti/360)*360;
    elseif hti<0
        hti = hti + abs(floor(hti/360))*360;
    end
    
    if bri_asterisk==0 && ari_prime==0
        hri = 0;
    elseif ari_prime==0 && bri_asterisk<0
        hri = -90;
    elseif ari_prime==0 && bri_asterisk>0
        hri = 90;
    elseif ari_prime<0 && bri_asterisk<0
        hri = (atan(bri_asterisk/ari_prime)-pi)/pi*180;
    elseif ari_prime<0 && bri_asterisk>=0
        hri = (atan(bri_asterisk/ari_prime)+pi)/pi*180;
    else
        hri = atan(bri_asterisk/ari_prime)/pi*180;
    end
    if hri>360
        hri = hri - floor(hri/360)*360;
    elseif hri<0
        hri = hri + abs(floor(hri/360))*360;
    end
    
    Li_prime_ave = (Lti_asterisk+Lri_asterisk)/2;
    Ci_prime_ave = (Cti_prime+Cri_prime)/2;
    if Cti_prime*Cri_prime==0
        hi_ave = hti+hri;
    else
        if abs(hti-hri)<=180
            hi_ave = (hti+hri)/2;
        else
            if hti+hri<360
                hi_ave = (hti+hri+360)/2;
            else
                hi_ave = (hti+hri-360)/2;
            end
        end
    end
    
    if Cti_prime*Cri_prime==0
        delta_hi = hti-hri;
    else
        if hri-hti<-180
            delta_hi = hri-hti+360;
        elseif hri-hti>180
            delta_hi = hri-hti-360;
        else
            delta_hi = hri-hti;
        end
    end
    
    % Calculate further intermediate values for the colour under test.
    S_Li = 1+(0.015*(Li_prime_ave-50)^2)/(sqrt(20+(Li_prime_ave-50)^2));
    S_Ci = 1+0.045*Ci_prime_ave;
    Ti = 1-0.17*cos((hi_ave-30)/180*pi)+0.24*cos((2*hi_ave)/180*pi)+0.32*cos((3*hi_ave+6)/180*pi)-0.2*cos((4*hi_ave-63)/180*pi);
    S_Hi = 1+0.015*Ci_prime_ave*Ti;
    Rci = 2*sqrt(Ci_prime_ave^7/(Ci_prime_ave^7+25^7));
    delat_theta_i = 30*exp(-((hi_ave-275)/25)^2);
    Rti = -Rci*sin(2*delat_theta_i/180*pi);
    
    % calculate the resulting CIEDE2000 difference values.
    Delta_Li = (Lri_asterisk-Lti_asterisk)/1/S_Li;
    Delta_Ci = (Cri_prime-Cti_prime)/1/S_Ci;
    Delta_Hi = (2*sin(delta_hi/2/180*pi)*sqrt(Cti_prime*Cri_prime))/1/S_Hi;
    
    % Coiour-Difference Formula:
    Dealta_E_asterisk(Ri) = sqrt(Delta_Li^2+Delta_Ci^2+Delta_Hi^2+Rti*Delta_Ci*Delta_Hi);

end

Dealta_Ea_asterisk = 1/24*sum(Dealta_E_asterisk(1:24));
Ra = 100-4.6*Dealta_Ea_asterisk;
% fprintf('Ra = %.2f\n',Ra);

Dealta_Ea_asterisk = (1/18*sum(Dealta_E_asterisk(1:18).^4))^(1/4);
Qa = 100/(1+(Dealta_Ea_asterisk/k)^p);
fprintf('Qa = %.2f\n',Qa);
