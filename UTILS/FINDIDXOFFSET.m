function [offset_idx, correlation] = FINDIDXOFFSET(A, B, windows)
%FINDTIMEOFFSET The function calculates the time offset between two
%timetables A and B that contain CO2 data. Be aware that columns 'T' of
%both timetables should denote the time column.



corr_coeff = zeros(1,1);

range = min(height(A), height(B));


for start = 1:range-windows-1
    
    comp_a = A(1:windows+1);
    comp_b = B(start:windows+start);
    corr_coeff(start) = corr(comp_a, comp_b);
    
end

[offset_idx, correlation] = max(corr_coeff);

end

