%% Ressitor calculations


% green = 180Ω, 1/4W
% UV = 6.8Ω, 3W

Vsupply = 12; 
Vf_values = [9 10.5 11.9];
I_values = [0.35]; 

% Create a matrix to store resistor values
Resistor_Values = zeros(length(Vf_values), length(I_values));

% Calculate resistor values
for i = 1:length(Vf_values)
    for j = 1:length(I_values)
        Vf = Vf_values(i);
        I = I_values(j);
        R = (Vsupply - Vf) / I; % Ohm's law
        Resistor_Values(i, j) = R;
    end
end

% Display results in a table
Resistor_Table = array2table(Resistor_Values, ...
    'VariableNames', strcat("I_", string(I_values * 1000), "mA"), ...
    'RowNames', strcat("Vf_", string(Vf_values), "V"));

disp('Calculated Resistor Values (Ohms):');
disp(Resistor_Table);


%

R = 4.7*1.05
I = (Vsupply - Vf_values) ./ R% Current in Amperes

P = I.^2 .* R % Power in Watts



% green = 180ohm, 1/4W