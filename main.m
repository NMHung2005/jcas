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
%% Optimize & Compare DE Variants
% --- Thiet lap chung cho viec toi uu ---

% 1. Dinh nghia cac tham so cho DE
de_params.pop_size = 150;     % Kich thuoc quan the
de_params.F = 0.8;            % He so dot bien
de_params.CR = 0.9;           % Ty le lai ghep
de_params.max_gens = 200;    % So the he toi da (cho vong so sanh)

% --- Tao vector trong so de "phat" sai so o bup phu ---
sidelobe_penalty = 100; % "Phat" sai so bup phu manh gap 100 lan
weights = ones(1, length(alpha));
sidelobe_indices_in_alpha = find(PdM(alpha) == 0);
weights(sidelobe_indices_in_alpha) = sidelobe_penalty;

% 2. Dinh nghia ham muc tieu
objective_func = @(w_vec) calculate_DE_cost(w_vec, M, Aq, PdM, alpha, weights);

% --- Chay so sanh cac bien the ---
variants_to_run = {'rand/1/bin', 'best/1/bin', 'rand-to-best/1/bin', 'rand/2/bin', 'best/2/bin'};
all_cost_histories = cell(length(variants_to_run), 1);
final_costs = zeros(length(variants_to_run), 1); % Luu sai so cuoi cung cua moi bien the
colors = lines(length(variants_to_run)); 

figure; % Tao figure moi cho do thi hoi tu
hold on;

fprintf('--- BAT DAU SO SANH CAC BIEN THE DE (200 a generation) ---\n');
for i = 1:length(variants_to_run)
    variant = variants_to_run{i};
    [~, cost_history] = optimize_DE(objective_func, M, de_params, variant);
    all_cost_histories{i} = cost_history;
    final_costs(i) = cost_history(end); % Luu lai sai so cuoi cung
    
    semilogy(1:de_params.max_gens, cost_history, 'DisplayName', variant, 'Color', colors(i,:), 'LineWidth', 1.5);
end
fprintf('--- KET THUC SO SANH ---\n\n');

% --- Hoan thien do thi so sanh ---
title('So sánh sự hội tụ của các biến thể DE');
xlabel('Thế hệ (Generation)');
ylabel('Log(Sai số tốt nhất)');
legend('show', 'Interpreter', 'none', 'Location', 'northeast');
grid on;
hold off;

% --- Tim ra bien the tot nhat va chay lai de lay ket qua cuoi cung ---
[~, best_variant_idx] = min(final_costs);
best_variant_name = variants_to_run{best_variant_idx};
fprintf('--- Bien the tot nhat trong so sanh la: %s ---\n', best_variant_name);

fprintf('--- Chay lai bien the tot nhat (%s) voi 1000 a generation de lay ket qua cuoi cung ---\n', best_variant_name);
de_params.max_gens = 1000; % Su dung lai so the he day du
[best_w_vector_real, ~] = optimize_DE(objective_func, M, de_params, best_variant_name);

% 4. Tai tao lai vector trong so phuc tu tu ket qua
W_ref_real_part = best_w_vector_real(1:M);
W_ref_imag_part = best_w_vector_real(M+1:end);
w_final_complex = W_ref_real_part + 1j * W_ref_imag_part;

% Chuan hoa ket qua cuoi cung
W_ref(:, 1) = w_final_complex / norm(w_final_complex);

%% Run ILS for comparison
fprintf('--- Bat dau toi uu bang Two-Step ILS ---\n');
iter_nr_ils = 10; % So vong lap cho ILS
% Goi ham ILS. W0 va P_init duoc khoi tao o phan %% Reference beam
W_ils = twoStepILS(iter_nr_ils, alpha, Aq, W0, P_init, PdM);
% Chuan hoa ket qua ILS de so sanh cong bang
W_ils = W_ils / norm(W_ils);
fprintf('--- ILS hoan thanh ---\n');

%% Plot optimized beams (DE vs ILS)
figure; % Tao figure moi de khong ve de len do thi so sanh hoi tu
hold on
grid on

% Ve cac duong tham chieu
plot(eqDir, 10*log10(PdM/max(PdM)), 'm-*', 'DisplayName', 'Desired')
plot(eqDir, 10*log10(P_refGen/max(P_refGen)), '--', 'Color', [0.3 0.3 0.3], 'DisplayName', 'Conventional')

% Ve ket qua cua DE
plot(eqDir, 10*log10(abs(W_ref'*Aq)/max(abs(W_ref'*Aq))), 'r', 'LineWidth', 1.5, 'DisplayName', 'DE Optimized')

% Ve ket qua cua ILS
plot(eqDir, 10*log10(abs(W_ils'*Aq)/max(abs(W_ils'*Aq))), 'b--', 'LineWidth', 1.5, 'DisplayName', 'ILS Optimized')

% Cau hinh do thi
title('So sánh búp sóng tối ưu: Desired, Conventional, DE, ILS');
xlabel("Equivalent directions")
ylabel("|A|, dB")
xlim([-1 1])
ylim([-40, 5])
legend('Location', 'northoutside', 'NumColumns', 4)
hold off

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
