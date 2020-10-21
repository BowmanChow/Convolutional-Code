function  y = KalmanFilter(z, A, Q, H, R, init)

% X(k) = A * X(k-1) + W(k) ,  Cov(W) = Q
% Z(k) = H * X(k) + V(k) , Cov(V) = R

if size(z, 1) == 1
    z = z.';
end

X_pre = zeros(size(z));
P_pre = X_pre;
X = X_pre;
Kg = X_pre;
P = X_pre;

X_pre(1) = init;
P_pre(1) = Q;
Kg(1) = P_pre(1) * H / (H^2 * P_pre(1) + R);
X(1) = X_pre(1) + Kg(1) * (z(1) - H * X_pre(1));
P(1) = (1 - Kg(1) * H) * P_pre(1);

for i = 2:length(z)
    X_pre(i) = A * X(i - 1);
    P_pre(i) = A^2 * P(i-1) + Q;
    Kg(i) = P_pre(i) * H / (H^2 * P_pre(i) + R);
    X(i) = X_pre(i) + Kg(i) * (z(i) - H * X_pre(i));
    P(i) = (1 - Kg(i) * H) * P_pre(i);
end

y = X;


