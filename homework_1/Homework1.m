%% Homework 1

clear all, clc, close all

load('Testdata.mat')

%n = (length(Undata))^(1/3);
%n = uint8(n);
n = 64;
L = 15;
x2 = linspace(-L,L,n+1);
x = x2(1:n);
y = x;
z = y;
[X,Y,Z]=meshgrid(x,y,z);

k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1];
ks=fftshift(k);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
Unt = 0;

for k = 1 : 20
    
    Un(:, :, :) = reshape(Undata(k,:), n, n, n);
    Unt = Unt + fft(Un(:, :, :));
    Unt = Unt / k;
    close all,
    isosurface(Kx,Ky,Kz,abs(Un),1)
    axis([-10 10 -10 10 -10 10]),
    grid on,
    drawnow
    pause(5)
    Unt = Unt * k;
    
end

%close all

%Unt = Unt / k;

%Untf = ifft(Unt);

%isosurface(X,Y,Z,abs(Untf), 0.4)


%Unf = 0;

%for k = 1
    
%    Un(:, :, :) = reshape(Undata(k,:), n, n, n);
%    Unf = Unf + Un(:, :, :);
    
%end

%Unf/k;

%isosurface(X,Y,Z,abs(Un), 0.4)



%for k = 1:2
%    
%    Un(:, :, :) = reshape(Undata(k,:), n, n, n);
%    t(1, k) = k-1;
%    close all
%    isosurface(X,Y,Z,abs(Un), 0.4)
%    
%end


%% Tips

clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
Uns = 0;
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
Uns = Uns + fft(Un);
Uns = Uns / j;
close all, isosurface(X,Y,Z,abs(Uns),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
pause(5)
Uns = Uns * j;
end

%%
clear; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);
Unf = zeros(64, 64, 64)
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n);
Unf = Unf + fft(Un);
Unf = Unf / j;
close all, isosurface(X,Y,Z,abs(Unf),1)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
pause(5)
Unf = Unf * j;
end