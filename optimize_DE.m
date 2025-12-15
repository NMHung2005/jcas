function [best_w_vec, best_cost_history] = optimize_DE(objective_func, M, de_params, variant)

pop_size = de_params.pop_size; % Kích thước quần thể
F = de_params.F;               % Hệ số đột biến
CR = de_params.CR;             % Tỷ lệ lai ghép
max_gens = de_params.max_gens; % Số thế hệ tối đa
dim = 2 * M;                   % Số chiều (phần thực + phần ảo)

% --- Mảng lưu lịch sử hội tụ ---
best_cost_history = zeros(max_gens, 1);

% --- Định nghĩa không gian tìm kiếm ---
lower_bound = -1 * ones(1, dim);
upper_bound =  1 * ones(1, dim);

% 1. KHỞI TẠO QUẦN THỂ
% Tạo ngẫu nhiên các vector trọng số trong không gian tìm kiếm
population = lower_bound + rand(pop_size, dim) .* (upper_bound - lower_bound);

% Đánh giá độ lỗi ban đầu cho từng cá thể
costs = zeros(pop_size, 1);
for i = 1:pop_size
    costs(i) = objective_func(population(i, :));
end

% 2. VÒNG LẶP TỐI ƯU
fprintf('Bat dau toi uu bang DE (bien the: %s)...\n', variant);
for gen = 1:max_gens
    % Tìm cá thể tốt nhất của thế hệ hiện tại (cần cho các biến thể 'best' và 'rand-to-best')
    [~, best_idx_pop] = min(costs);
    best_vector = population(best_idx_pop, :);

    for i = 1:pop_size
        % --- Đột biến (Mutation) ---
        % Chọn các cá thể ngẫu nhiên r1, r2, ... (khác i)
        indices = 1:pop_size;
        indices(i) = []; % Loại bỏ chỉ số của cá thể hiện tại
        
        switch variant
            case 'rand/1/bin'
                r = indices(randperm(length(indices), 3));
                mutant_vector = population(r(1), :) + F * (population(r(2), :) - population(r(3), :));
            
            case 'best/1/bin'
                r = indices(randperm(length(indices), 2));
                mutant_vector = best_vector + F * (population(r(1), :) - population(r(2), :));

            case 'rand-to-best/1/bin'
                r = indices(randperm(length(indices), 2));
                mutant_vector = population(i, :) + F * (best_vector - population(i, :)) + F * (population(r(1), :) - population(r(2), :));

            case 'rand/2/bin'
                r = indices(randperm(length(indices), 5));
                mutant_vector = population(r(1), :) + F * (population(r(2), :) - population(r(3), :)) + F * (population(r(4), :) - population(r(5), :));

            case 'best/2/bin'
                r = indices(randperm(length(indices), 4));
                mutant_vector = best_vector + F * (population(r(1), :) - population(r(2), :)) + F * (population(r(3), :) - population(r(4), :));
            
            otherwise
                error('Bien the DE "%s" khong duoc ho tro.', variant);
        end
        
        % Xử lý nếu vượt ngoài không gian tìm kiếm
        mutant_vector = max(mutant_vector, lower_bound);
        mutant_vector = min(mutant_vector, upper_bound);
        
        % --- Lai ghép (Crossover) ---
        trial_vector = population(i, :);
        j_rand = randi(dim);
        for j = 1:dim
            if rand() < CR || j == j_rand
                trial_vector(j) = mutant_vector(j);
            else
                trial_vector(j) = population(i, j);
            end
        end
        
        % --- Lựa chọn (Selection) ---
        trial_cost = objective_func(trial_vector);
        if trial_cost < costs(i)
            population(i, :) = trial_vector;
            costs(i) = trial_cost;
        end
    end
    
    % Ghi lại và hiển thị tiến trình
    best_cost_history(gen) = min(costs);
    if mod(gen, 10) == 0
        fprintf('The he %d, Sai so tot nhat: %f\n', gen, best_cost_history(gen));
    end
end

% 3. TRẢ VỀ KẾT QUẢ TỐT NHẤT
[~, best_idx] = min(costs);
best_w_vec = population(best_idx, :);
fprintf('DE (%s) hoan thanh! Sai so cuoi cung: %f\n', variant, min(costs));

end
