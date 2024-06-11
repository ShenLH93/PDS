%% Polarization-Driven Solution for Mitigating Scattering 
%% and Uneven Illumination in Underwater Imagery
%% Linghao Shen;Mohamed Reda;Xun Zhang;Yongqiang Zhao;Seong G. Kong
%% Accepted by IEEE Transactions on Geoscience and Remote Sensing, 2024


clear all;
close all;


testimg_dir = '.\inputs\'; 
save_dir='.\outputs\';

files = dir(testimg_dir);
size0 = size(files);
length_img = size0(1);

for i = 3:1:size0(1)
        % read image
        fileName = strcat(testimg_dir,files(i).name);   
        Image = im2double(imread(fileName));

        % Calculate I0-I35,This code takes the following sample as 
        % an example and can be replaced by the demosaicing algorithm.
        Image =Image(:,:,1);
        [m,n,~] = size(Image);
        I0   = Image(1:2:m,1:2:n,1);
        I45  = Image(1:2:m,2:2:n,1);
        I90  = Image(2:2:m,2:2:n,1);
        I135 = Image(2:2:m,1:2:n,1);
        % [I0,I45,I90,I135] = Newton_Polynomial_Interpolation(I);

        % Calculate Stokes Vector
        I = 0.5*(I0+I45+I90+I135);
        Q = I0-I90;
        U = I45-I135;

        % Calculate Imax and Imin
        Imax = 0.5*(I +sqrt(Q.*Q+U.*U));
        Imin = 0.5*(I -sqrt(Q.*Q+U.*U));
        
        % Calculate Illumination I, reflectance R, 
        % transmission map t and degree of liner polarization DoLP.
        [I, R,t ,p] = monoPDS(Imax,Imin);

        % Adjust Illumination
%         Iy = I.^(0.5);
%         Result = Iy.*R;

        I_max = max(max(max(I)),1);
        I_min = min(min(I));
        I_mean = mean(mean(I));
        PX = [I_min I_mean I_max];
        PY = [I_min+2/3*(I_mean-I_min) I_mean (I_mean+1/3*(I_max-I_mean))];
        Ik = polyfit(PX,PY,2);
        I2 = polyval(Ik,I);
        Result = I2.*R;

        maxS = max(max(Result));
        minS = min(min(Result));
        Result = uint8(255.*(Result-minS)./(maxS-minS));

        % Save results
        imwrite(Result, strcat(save_dir,files(i).name));
end



