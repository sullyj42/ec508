clear all; close all; 
b = [0 1 1 1]'; 
A = [1 0; 0 1; 1 1]; 
for i = 1:2
   c1(:,i)    =  A*b( (1:2).*i );
   c((1:3)*i) =  A*b( (1:2).*i );
end
c = mod(c, 2);
disp(c); 

x = (-1).^c; 
disp(x); 

T = 5; % Arbitrary? 
t = linspace(-10, 100); 
xb = zeros(size(t)); 
for i = 1:length(c)
    cond = p0(t, i, T); 
    temp = x(i).*cond;
    xb = xb+temp; 
end
figure(1); 
subplot(3,2,1:2); plot(t, xb); title('BPSK'); 
figure(2); hold on; plot(x, 0.*x, 'o');
clear xb x; xb = zeros(size(t)); 

for i = 1:length(c)/2
    cond = p0(t, i, T); 
    if     isequal(c((2*i-1):(2*i)),  [0, 0])
        x(i) = 3; 
    elseif isequal(c((2*i-1):(2*i)),  [0, 1])
        x(i) = 1; 
    elseif isequal(c((2*i-1):(2*i)),  [1, 1])
        x(i) = -1; 
    elseif isequal(c((2*i-1):(2*i)),  [1, 0])
        x(i) = -3; 
    else
        warning('bb'); 
    end

    
    temp = x(i).*cond;
    xb = xb+temp; 
end
figure(1); 
subplot(3,2,3:4); plot(t, xb); % Phase amplitude modulation
figure(2); hold on; plot(x, 0.*x, '+'); title('Constellation'); 


for i = 1:length(c)/2
    cond = p0(t, i, T); 
    if     isequal(c((2*i-1):(2*i)),  [0, 0])
        x(i) = 1+sqrt(-1); 
    elseif isequal(c((2*i-1):(2*i)),  [0, 1])
        x(i) = 1-sqrt(-1); 
    elseif isequal(c((2*i-1):(2*i)),  [1, 1])
        x(i) =  -1 - sqrt(-1); 
    elseif isequal(c((2*i-1):(2*i)),  [1, 0])
        x(i) = -1 + sqrt(-1); 
    else
        warning('bb'); 
    end
    
    temp = x(i).*cond;
    xb = xb+temp; 
end
figure(1); subplot(3,2,5); plot(t, real(xb)); title('QPSK'); xlabel('Real');  % Phase amplitude modulation
subplot(3,2,6); plot(t, imag(xb)); title('QPSK'); xlabel('Imaginary'); 

figure(2); hold on; plot(x, '*'); title('Constellation'); xlabel('Re'); ylabel('Im'); 
legend('BPSK', 'PAM', 'QPSK'); axis equal; 
xlim([-2 2]);
ylim([-2 2]); 

function y = p0(t, m, T)
        cond = ((t-m*T) >= 0).*((t-m*T) <= m*T); 
        y = cond; 
    
end


