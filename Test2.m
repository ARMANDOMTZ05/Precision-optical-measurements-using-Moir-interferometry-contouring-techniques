%% Método Amalia corrimiento de fase con interferómetro

Bfase1 = load("BCD5_fase0.txt");
Bfase2 = load("BCD5_fase1.txt");
Bfase3 = load("BCD5_fase2.txt");
Bfase4 = load("BCD5_fase3.txt");
OFase1 = load("OBJ5_fase0.txt");
OFase2 = load("OBJ5_fase1.txt");
OFase3 = load("OBJ5_fase2.txt");
OFase4 = load("OBJ5_fase3.txt");


%% 

% Filtramos las frecuencias bajas
FTBfase1 = fftshift(fft2(Bfase1/max(max(Bfase1))));
FTBfase2 = fftshift(fft2(Bfase2/max(max(Bfase2))));
FTBfase3 = fftshift(fft2(Bfase3/max(max(Bfase3))));
FTBfase4 = fftshift(fft2(Bfase4/max(max(Bfase4))));
FTOfase1 = fftshift(fft2(OFase1/max(max(OFase1))));
FTOfase2 = fftshift(fft2(OFase2/max(max(OFase2))));
FTOfase3 = fftshift(fft2(OFase3/max(max(OFase3))));
FTOfase4 = fftshift(fft2(OFase4/max(max(OFase4))));

R = 25;
[row,col] = size(Bfase1);
[xGrid,yGrid] = meshgrid(1:col,1:row);
Bw = sqrt((xGrid - col/2).^2 + (yGrid - row/2).^2) <= R;

Filtered_FFTBfase1 = FTBfase1.*Bw;
Filtered_FFTBfase2 = FTBfase2.*Bw;
Filtered_FFTBfase3 = FTBfase3.*Bw;
Filtered_FFTBfase4 = FTBfase4.*Bw;
Filtered_FFTOfase1 = FTOfase1.*Bw;
Filtered_FFTOfase2 = FTOfase2.*Bw;
Filtered_FFTOfase3 = FTOfase3.*Bw;
Filtered_FFTOfase4 = FTOfase4.*Bw;

NBfase1 = real(ifft2(ifftshift(Filtered_FFTBfase1)));
NBfase2 = real(ifft2(ifftshift(Filtered_FFTBfase2)));
NBfase3 = real(ifft2(ifftshift(Filtered_FFTBfase3)));
NBfase4 = real(ifft2(ifftshift(Filtered_FFTBfase4)));

NOfase1 = real(ifft2(ifftshift(Filtered_FFTOfase1)));
NOfase2 = real(ifft2(ifftshift(Filtered_FFTOfase2)));
NOfase3 = real(ifft2(ifftshift(Filtered_FFTOfase3)));
NOfase4 = real(ifft2(ifftshift(Filtered_FFTOfase4)));


%% 

%Recuperamos la fase envuelta

Bphase = atan2((NBfase4 - NBfase2),(NBfase1 - NBfase3));
Ophase = atan2((NOfase4 - NOfase2),(NOfase1 - NOfase3));

T_phase = Bphase - Ophase;
NFT_phase = NFBphase - NFOphase;

%Desenvolvemos la fase
Unwraped_phaseB = phase_unwrap(Bphase); %Metodo 2

Unwraped_phaseO = phase_unwrap(Ophase); %Método 2
%% 

%Hacemos la diferencia de fase
Tphase = Unwraped_phaseB - Unwraped_phaseO;
Tphase2 = Tphase(39:643,155:759);

%% 


% Reconstrucción de la dimensión z
theta = acos(0.15/0.17);
Z = (Tphase2/(2*pi))*((0.0084)/tan(theta));
Z2 = (Tphase2/(2*pi))*(0.0084/tan(theta));

[rowZ,colZ] = size(Z);

x1 = linspace(0,0.0762,colZ);
y1 = linspace(0,0.0593,rowZ);



Z(Z < 0) = 0;

[rowZ2,colZ2] = size(Z);

CM = 0.0512;

x2 = linspace(-CM/2,CM/2,colZ2);
y2 = linspace(-CM/2,CM/2,rowZ2);

% Filtrado para eliminar el ruido de la reoncstrucción
[X2,Y2] = meshgrid(x2,y2);

FC = Z;
FC(X2.^2 + Y2.^2 >= 0.024.^2) = 0;
FC(X2.^2 + Y2.^2 >= 0.0253.^2 & FC > 0.0035) = 0;
FC(X2.^2 + Y2.^2 >= 0.021.^2 & FC > 0.01) = 0;


SNZ = smoothdata(FC,2);


%% Parte de las gráficas

close all

% figure(1)
% tiledlayout(2,2)
% 
% nexttile
% surf(Unwraped_phase3,EdgeColor="none")
% title('Unwraped phase (Método 1) sin filtro')
% colormap gray
% 
% nexttile
% surf(Unwraped_phase4,EdgeColor="none")
% title('Unwraped phase (Método 2) sin fitro')
% 
% nexttile
% surf(Unwraped_phase, EdgeColor="none")
% title('Unwraped phase (Método 1) con filtro')
% colormap gray
% 
% nexttile
% surf(Unwraped_phase2,EdgeColor="none")
% title('Unwraped phase (Método 2) con fitro')
% colormap gray
% 
figure(1)
surf(Unwraped_phaseO,EdgeColor="none")
title('Unwraped phase (Método 2) con fitro')
colormap gray


figure(2)
surf(x2,y2,Z,EdgeColor="none")
colormap gray

figure(3)
surf(x2,y2,FC,EdgeColor="none")
xlabel('X(m)', Interpreter='latex')
ylabel('Y(m)', Interpreter='latex')
zlabel('Z(m)', Interpreter='latex')
grid off
colormap gray
% xlim([-CM/2 CM/2])
% ylim([-CM/2 CM/2])
% zlim([0 CM])
axis equal

% figure(4)
% tiledlayout(2,2)
% nexttile
% imagesc(NBfase1);colormap gray
% title('$\delta = 0$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NBfase2);colormap gray
% title('$\delta = \pi/2$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NBfase3);colormap gray
% title('$\delta = \pi$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NBfase4);colormap gray
% title('$\delta = 3\pi/4$',Interpreter='latex')
% axis off

figure(4)
imagesc(NBfase1);colormap gray
title('$\delta = 0$',Interpreter='latex')
axis off

figure(5)
imagesc(NBfase2);colormap gray
title('$\delta = \pi/2$',Interpreter='latex')
axis off

figure(6)
imagesc(NBfase3);colormap gray
title('$\delta = \pi$',Interpreter='latex')
axis off

figure(7)
imagesc(NBfase4);colormap gray
title('$\delta = 3\pi/2$',Interpreter='latex')
axis off

% figure(5)
% tiledlayout(2,2)
% nexttile
% imagesc(NOfase1);colormap gray
% title('$\delta = 0$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NOfase2);colormap gray
% title('$\delta = \pi/2$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NOfase3);colormap gray
% title('$\delta = \pi$',Interpreter='latex')
% axis off
% nexttile
% imagesc(NOfase4);colormap gray
% title('$\delta = 3\pi/4$',Interpreter='latex')
% axis off

figure(8)
imagesc(NOfase1);colormap gray
title('$\delta = 0$',Interpreter='latex')
axis off

figure(9)
imagesc(NOfase2);colormap gray
title('$\delta = \pi/2$',Interpreter='latex')
axis off

figure(10)
imagesc(NOfase3);colormap gray
title('$\delta = \pi$',Interpreter='latex')
axis off

figure(11)
imagesc(NOfase4);colormap gray
title('$\delta = 3\pi/2$',Interpreter='latex')
axis off


% figure(12)
% tiledlayout(1,2)
% nexttile
% imagesc(Tphase2)
% colormap gray
% nexttile
% % surf(Diff_phase,EdgeColor="none");

figure(12)
imagesc(Bphase)
axis equal
axis off
colormap gray

figure(13)
imagesc(Ophase)
colormap gray
axis equal
axis off






%% 


function phi = phase_unwrap(psi, weight)
    if (nargin < 2) % unweighted phase unwrap
        % get the wrapped differences of the wrapped values
        dx = [zeros([size(psi,1),1]), wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])];
        dy = [zeros([1,size(psi,2)]); wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])];
        rho = diff(dx, 1, 2) + diff(dy, 1, 1);
        
        % get the result by solving the poisson equation
        phi = solvePoisson(rho);
        
    else % weighted phase unwrap
        % check if the weight has the same size as psi
        if (~all(size(weight) == size(psi)))
            error('Argument error: Size of the weight must be the same as size of the wrapped phase');
        end
        
        % vector b in the paper (eq 15) is dx and dy
        dx = [wrapToPi(diff(psi, 1, 2)), zeros([size(psi,1),1])];
        dy = [wrapToPi(diff(psi, 1, 1)); zeros([1,size(psi,2)])];
        
        % multiply the vector b by weight square (W^T * W)
        WW = weight .* weight;
        WWdx = WW .* dx;
        WWdy = WW .* dy;
        
        % applying A^T to WWdx and WWdy is like obtaining rho in the unweighted case
        WWdx2 = [zeros([size(psi,1),1]), WWdx];
        WWdy2 = [zeros([1,size(psi,2)]); WWdy];
        rk = diff(WWdx2, 1, 2) + diff(WWdy2, 1, 1);
        normR0 = norm(rk(:));
        
        % start the iteration
        eps = 1e-8;
        k = 0;
        phi = zeros(size(psi));
        while (~all(rk == 0))
            zk = solvePoisson(rk);
            k = k + 1;
            
            if (k == 1) pk = zk;
            else 
                betak = sum(sum(rk .* zk)) / sum(sum(rkprev .* zkprev));
                pk = zk + betak * pk;
            end
            
            % save the current value as the previous values
            rkprev = rk;
            zkprev = zk;
            
            % perform one scalar and two vectors update
            Qpk = applyQ(pk, WW);
            alphak = sum(sum(rk .* zk)) / sum(sum(pk .* Qpk));
            phi = phi + alphak * pk;
            rk = rk - alphak * Qpk;
            
            % check the stopping conditions
            if ((k >= numel(psi)) || (norm(rk(:)) < eps * normR0)) break; end;
        end
    end
end
function phi = solvePoisson(rho)
    % solve the poisson equation using dct
    dctRho = dct2(rho);
    [N, M] = size(rho);
    [I, J] = meshgrid([0:M-1], [0:N-1]);
    dctPhi = dctRho ./ 2 ./ (cos(pi*I/M) + cos(pi*J/N) - 2);
    dctPhi(1,1) = 0; % handling the inf/nan value
    
    % now invert to get the result
    phi = idct2(dctPhi);
    
end
% apply the transformation (A^T)(W^T)(W)(A) to 2D matrix
function Qp = applyQ(p, WW)
    % apply (A)
    dx = [diff(p, 1, 2), zeros([size(p,1),1])];
    dy = [diff(p, 1, 1); zeros([1,size(p,2)])];
    
    % apply (W^T)(W)
    WWdx = WW .* dx;
    WWdy = WW .* dy;
    
    % apply (A^T)
    WWdx2 = [zeros([size(p,1),1]), WWdx];
    WWdy2 = [zeros([1,size(p,2)]); WWdy];
    Qp = diff(WWdx2,1,2) + diff(WWdy2,1,1);
end
