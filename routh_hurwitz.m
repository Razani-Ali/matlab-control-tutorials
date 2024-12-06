function rh_table = routh_hurwitz(coeffs)
    % Input: coeffs - a symbolic vector of coefficients of the polynomial
    % Output: rh_table - Routh-Hurwitz table as a symbolic matrix

    % Degree of the polynomial
    n = length(coeffs) - 1;

    % Dimensions of the Routh-Hurwitz table
    rows = n + 1;
    cols = ceil(n / 2) + 1;

    % Initialize Routh-Hurwitz table with zeros
    rh_table = sym(zeros(rows, cols + 1));
    
    k = 0;
    % Fill the first two rows
    for i = 1:length(coeffs)
        k = k + mod(i, 2);
        if mod(i, 2) == 1
            rh_table(1, k) = coeffs(i);
        else
            rh_table(2, k) = coeffs(i);
        end
    end

    % Compute the rest of the table
    for i = 3:rows
        for j = 1:cols - 1
            % Formula:
            % (two up one right) - (two up first column) * (one up one right) / (one up first column)
            denominator = rh_table(i - 1, 1);
            if denominator == 0
                % Handle divide-by-zero with a small epsilon perturbation
                denominator = sym('eps'); 
            end
            rh_table(i, j) = rh_table(i - 2, j + 1) - rh_table(i - 2, 1) * rh_table(i - 1, j + 1) / denominator;
        end

        % Check if the current row is all zeros
        if all(rh_table(i, :) == 0)
            ind = n - i + 1; % Current polynomial order
            mul = ind + 1:-2:0;
            rem = cols - length(mul) + 1; % Remaining columns to fill
            mul = [mul, zeros(1, rem)]; % Multiplication factors
            fprintf('Row %d is all zeros. Using polynomial derivatives:\n', i);
            fprintf('Coefficients of pervious row:')
            disp(rh_table(i - 1, :)); % Convert to numeric for display
            fprintf('Coefficients * s raised to a power of: %s\n', mat2str(mul));

            % Replace current row with derived polynomial row
            rh_table(i, :) = mul .* rh_table(i - 1, :);
        end
    end

    % Trim extra columns with all zeros
    if mod(n, 2) == 1
        rh_table = rh_table(:, 1:end - 2);
    else
        rh_table = rh_table(:, 1:end - 1);
    end
end
