function cost = calculate_DE_cost(w_real_imag, M, Aq, PdM, alpha, weights)
% HÀM MỤC TIÊU CHO THUẬT TOÁN DE
% w_real_imag: Vector chứa phần thực và phần ảo của các trọng số w
% M: Số phần tử ăng-ten
% Aq: Ma trận đáp ứng mảng
% PdM: Búp sóng mong muốn
% alpha: Các chỉ số hướng cần tối ưu
% weights: Trọng số cho các điểm sai số

% 1. Tái tạo lại vector trọng số phức w từ vector thực và chuyển vị thành vector cột
w_complex = (w_real_imag(1:M) + 1j * w_real_imag(M+1:end)).';

% 2. Chuẩn hóa vector w để đảm bảo ràng buộc công suất (tùy chọn nhưng nên có)
w_complex = w_complex / norm(w_complex);

% 3. Tính búp sóng được tạo ra tại các hướng alpha
pattern_generated = abs(w_complex' * Aq(:, alpha));

% 4. Lấy búp sóng mong muốn tại các hướng alpha
pattern_desired = PdM(alpha);

% 5. Tính sai số có trọng số (weighted sum of squared errors)
cost = sum(weights .* (pattern_generated - pattern_desired).^2);

end
