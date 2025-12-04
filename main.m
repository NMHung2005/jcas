close all;
clear;
clc;
%% Variable initialization
theta = (-90:0.1:90-0.1)*pi/180; % rad
lambda = 1; % Wavelength
% It is assumed that the arrays dedicated to communication and sensing have
% the same number of elements.
% Number of array elements
M = 12;
% Sensing & communication steering vectors
A = generateSteeringVector(theta, M, lambda);
desDirs_c = 0.0;
% Container to save the reference BF vector, which is a basis for the
% multibeam modelling
W_ref = zeros(M, size(desDirs_c, 2));
%% Array response for the equivalent directions
% Q*phi is half of the number of values in the equivalent directions grid
Q = 160; % With Q=40, 12 ULA curve is exactly as in the paper
% Quanization step
phi = 1;
% Equivalent scanning directions (value of sin)
eqDir = -1:phi/Q:1-phi/Q;
% Array response associated with the equivalent directions
Aq = generateQuantizedArrResponse(M, eqDir);
%% Reference beam
% Generate desired and reference patterns as well as intial BF vector
[PdM, P_refGen, W0] = generateDesPattern(eqDir, sin(desDirs_c), Aq);
% [PdM, P_refGen, W0] = generateDesPattern(theta, desDirs_c, A);
P_init = ones(size(eqDir));
PM = P_init;
% Indices of the eq. directions that are to be approximated
alpha = sort([find(ismember(eqDir, eqDir(1:4:end))), find(PdM)]);
%% Optimize
% --- Thiet lap va chay Differential Evolution (voi ham muc tieu co trong so, khong tiem W0) ---

% 1. Dinh nghia cac tham so cho DE
de_params.pop_size = 100;     % Kich thuoc quan the
de_params.F = 0.8;            % He so dot bien
de_params.CR = 0.9;           % Ty le lai ghep
de_params.max_gens = 1000;    % So the he toi da

% --- Tao vector trong so de "phat" sai so o bup phu ---
sidelobe_penalty = 10; % "Phat" sai so bup phu manh gap 10 lan
weights = ones(1, length(alpha));
% Tim cac diem alpha nam trong vung bup phu (noi PdM = 0)
sidelobe_indices_in_alpha = find(PdM(alpha) == 0);
% Tang trong so o cac diem do
weights(sidelobe_indices_in_alpha) = sidelobe_penalty;

% 2. Dinh nghia ham muc tieu, truyen ca `weights` vao
objective_func = @(w_vec) calculate_DE_cost(w_vec, M, Aq, PdM, alpha, weights);

% 3. Goi ham toi uu DE (phien ban khong co W0)
best_w_vector_real = optimize_DE(objective_func, M, de_params);

% 4. Tai tao lai vector trong so phuc tu ket qua
W_ref_real_part = best_w_vector_real(1:M);
W_ref_imag_part = best_w_vector_real(M+1:end);
w_final_complex = W_ref_real_part + 1j * W_ref_imag_part;

% Chuan hoa ket qua cuoi cung
W_ref(:, 1) = w_final_complex / norm(w_final_complex);
%% Plot optimized beams
plot(eqDir, zeros(size(eqDir)));
hold on
plot(eqDir, 10*log10(PdM/max(PdM)), 'm-*')
hold on
plot(eqDir, 10*log10(P_refGen/max(P_refGen)), '--black')
hold on
plot(eqDir, 10*log10(abs(W_ref'*Aq)/max(abs(W_ref'*Aq))), 'r')
legend('Initial', 'Desired', 'Conventional 12-element ULA', 'Optimized', ...
    'Location', 'northoutside', 'NumColumns', 4)
xlabel("Equivalent directions")
ylabel("|A|, dB")
xlim([-1 1])
ylim([-35, 1])
%matlab2tikz('InitDesConvOpt.tex', 'height', '5cm', 'width', '8cm', ...
%    'showInfo', false)
%% Displace the reference beam
%spacing = 2*asin(1.2/M); % Fractional part may spoil the subbeam location
spacing=0.2;
deltas = -0.8:spacing:0.8; % approx. -53...53 degrees
%deltas = [0, sin(10.8*pi/180)]; % For Mengshuai
W_dd = zeros(M, size(deltas, 2));
for i=1:size(deltas, 2)
    W_dd(:, i) = displacePattern(W_ref, deltas(i), M);
end%for
%% Plot the beams coming from the optimized beamforming vectors
figure;
hold on
for i = 1:size(deltas, 2)
    plot(eqDir, 10*log10(abs(W_dd(:, i)'*Aq)/max(abs(W_dd(:, i)'*Aq))))    
end%for
xlim([-1 1])
ylim([-35, 0])
xlabel("Equivalent directions")
ylabel("|A|, dB")
grid on
%% Combination of the sensing and communicating beams
% Communication-sensing trade-off parameter
ro = 0.5; 
% Method 1
% Normalize beams
W_dd = W_dd./vecnorm(W_dd, 2, 2);
W_t = zeros(M, size(deltas, 2)-1);
% It is assumed that central beam is the one used for communication
comBeamIdx = cast(size(deltas, 2)/2, 'uint32');
%comBeamIdx = 1; % for Mengshuai
j=1;
for i = 1:size(W_dd, 2)
    if i~=comBeamIdx
        W_t(:, j) = sqrt(ro)*W_dd(:, comBeamIdx) + sqrt(1-ro)*W_dd(:, i);
        j = j + 1;
    end%if
end%for
%% Plot a combined beam
plot(eqDir, 10*log10(abs(W_t(:, 1)'*Aq)/max(abs(W_t(:, 1)'*Aq))))
xlabel("\theta, rad")
ylabel("|A(\theta)|, dB")
xlim([-1 1])
grid on
%matlab2tikz('Multibeam.tex', 'height', '5cm', 'width', '8cm', ...
%    'showInfo', false)
%% Plot the combined beams
figure;
hold on;
for i = 1:size(W_t, 2)
    % Equivalent directions
    %plot(eqDir, 10*log10(abs(W_t(:, i)'*Aq)));
    % Angle, rad
    plot(theta, 10*log10(abs(W_t(:, i)'*A)/max(abs(W_t(:, i)'*A))));
end%for\
%title(sprintf("M=%d", M))
xlabel("\theta, rad")
ylabel("|A(\theta)|, dB")
grid on
xlim([-pi/2, pi/2])
ylim([-15, 0])
%matlab2tikz('Multiple_multibeams.tex', 'height', '5cm', 'width', '8cm', ...
    %'showInfo', false)
