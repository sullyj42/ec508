% b[m] -->     (c[m])    -->        (x[m])      --> (xb(t)
% bits --> channel code  --> constellation map  --> D/A

b      = [0 1 1 1]; 
chCode = [1 0; 0 1; 1 1]; 


if mod(length(b), 2) % Lazy way to do this
    disp('Logic not implimented to code odd-bits');
end
n = 1; 
for i = 1:floor(length(b)/2)
    
    C(n:n+2) = chCode*b(i:(i+1))';
    n = n + 3; % Hardcode this better later
end
C = mod(C, 2); 

