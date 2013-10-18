Img=imread('gourd.bmp');
% Img=imread('E:\geo\testbmps\ls2.bmp');
Img=double(Img(:,:,1));
%% parameter setting
timestep=1;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=5;
iter_outer=20;
lambda=5; % coefficient of the weighted length term L(phi)
alfa=-3;  % coefficient of the weighted area term A(phi)
epsilon=1.5; % papramater that specifies the width of the DiracDelta function

% sigma=.8;    % scale parameter in Gaussian kernel
% G=fspecial('gaussian',15,sigma); % Caussian kernel
% Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
% Img_smooth=medfilt2(Img,[3 3]);
Img_smooth=Img;
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.
c0=2;
initialLSF = c0*ones(size(Img));
% generate the initial region R0 as two rectangles
%initialLSF(25:35,20:25)=-c0; 
initialLSF(25:35,40:50)=-c0;
% initialLSF(150:200,150:200)=-c0; 
% initialLSF(25:35,40:50)=-c0;
phi=initialLSF;
potential=2;  
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end  
for n=1:200

f=phi;
g0 = f;
[nrow,ncol] = size(f);
[vx, vy]=gradient(g);
g0([1 nrow],[1 ncol]) = g0([3 nrow-2],[3 ncol-2]);  
g0([1 nrow],2:end-1) = g0([3 nrow-2],2:end-1);          
g0(2:end-1,[1 ncol]) = g0(2:end-1,[3 ncol-2]);  
phi=g0;

    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); % add a small positive number to avoid division by zero
Ny=phi_y./(s+smallNumber);
[nxx,junk]=gradient(Nx);  
[junk,nyy]=gradient(Ny);
f=nxx+nyy;

    curvature=f;
          distRegTerm = 4*del2(phi)-curvature;  % compute distance 
     
%diracPhi=Dirac(phi,epsilon);
x=phi;
sigma=1.5;
%function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

diracPhi=f;
    areaTerm=diracPhi.*g; % balloon/pressure force
    edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    phi=phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm);
%imshow(phi);
imwrite(phi,strcat(num2str(n),'.bmp'));
end;