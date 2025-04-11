    
fprintf('Test: 000%%')
for i = 1:20
    fprintf('\b\b\b\b\b')
    fprintf(' %03d%%',floor(i))
end
fprintf('\n')