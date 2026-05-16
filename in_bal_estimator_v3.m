function [tauArray, In_balArray,AD4SpecificTau] = in_bal_estimator_v3(R1,R2,T1,T2,N1,N2,Nb,I1,A20_gain,S_ab,DeltaT_ab,FlagSpecificTau,SpecificTau)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT PARAMETERS (IMPORTED FROM THE GUI)
        % R1 and R2 - resistors on the master and slave branches, respectively (in ohms)
            % NOTE: R2 = 12906.4037297 ohms corresponds to QHR, i = 2
        % T1 and T2 - temperatures of R1 and R2, respectively (given in °C) 
        % N1, N2, Nb = number of turns of the master, slave, and balance windings, respectively (integer)
        % I1 - nominal value of the current applied to the master branch (in mA)
        % A20_gain = null detector analog gain
        % S_ab = Value of the Seebeck coefficient - it depends on the materials. Given in uV/K
        % DeltaT_ab = Value of the Seebeck coefficient - it depends on the materials. Given in uV/K
        %%%%%       NOTE -- IT IS ASSUMED THAT S_ab AND DeltaT_ab ARE THE       %%%%%
        %%%%%        SAME ACROSS ALL 16 CONNECTORS (4-WIRE MEASUREMENTS)        %%%%%
    
    % The following parameters (imported from the GUI) need to be converted to numeric
    R1 = str2double(R1);
    R2 = str2double(R2);
    N1 = str2double(N1);
    N2 = str2double(N2);
    Nb = str2double(Nb);
    A20_gain = 10^str2double(A20_gain);
    
    % The original SpecificTau value (imported from the GUI)
    SpecificTauOrig = SpecificTau;      % ... is stored in an auxiliary variable (just for display)
    SpecificTau = round(SpecificTau);   % ... needs to be rounded
            
    % The following parameter (imported from the GUI) needs to be converted to boolean
    FlagSpecificTau_x = false;
    if FlagSpecificTau == "Sim"
        FlagSpecificTau_x = true;       % True if the user wants to evaluate the Allan deviation for an specific tau
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %%%%% CONSTANTS %%%%%
    kB = 1.380649*10^-23;               % Boltzmann constant (J/K)
    AD_A20_voltage = 0.14;              % Null detector Allan deviation (voltage noise) for tau = 10 s (nV)
    AD_A20_current = 0.20;              % Null detector Allan deviation (current noise) for tau = 10 s (pA)
    AD_SQUID_a = 2.7*10^-11;            % SQUID Allan deviation given in the slope-intercept form (y = _ax + _b)
    AD_SQUID_b = 3.4*10^-10;            % SQUID Allan deviation given in the slope-intercept form (y = _ax + _b)
    N_ref = 32;                         % Number of winding turns to which the SQUID noise is scaled to
    Thermal_EMF_baseline = 0.2*10^-9;   % Least thermal EMF of all. It corresponds to a connection with Seebeck
                                        % coefficient of 0.2 uV/K and temperature difference of 1 mK. Estimate = 0.2 (nV)
    AD_floor = 0.1*10^-9;               % 1/f noise floor of connections in terms of Allan deviation. Estimate = 0.1 (nV)
    F_res = 1.25;                       % Resolution factor, i.e., DAC bit-to-ADC bit ratio
    Rb = 10^6;                          % Nominal value (in ohms) of the balance branch resistor (Rb = 1 Mohm)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% SETTING UP THINGS %%%%%
    tauArray = 10:1:100;                                        % create x-axis -- time (tau), given in s
        % Import I1 Allan deviations into a matrix:
    I1_ADMatrix = readmatrix('I_1_Allan_deviations.xlsx', 'Sheet', I1, 'Range', 'A2:B95');
    
    if R2 == 12906.4037297
        FlagQHR = true;
    else
        FlagQHR = false;
    end
    
    AD4SpecificTau = "";
    if ((SpecificTau < 10 || SpecificTau > 100) && FlagSpecificTau_x)
        AD4SpecificTau = "ATENÇÃO: o tempo Tau especificado está fora do intervalo válido !!!";
    end

    %%%%% NOISE INTENSITY COEFFICIENTS %%%%%
        % Noise intensity coefficients h0 of the resistors
    h0_R1 = 4*kB*R1*(T1+273.15);                                % (V^2/Hz)
    if FlagQHR == false
        h0_R2 = 4*kB*R2*(T2+273.15);                            % (V^2/Hz)
    else
        h0_R2 = 0;
    end

        % Noise intensity coefficients of the null detector
    h0_A20_voltage = 2*10*(AD_A20_voltage*10^-9)^2;             % (V^2/Hz)
    h0_A20_current = 2*10*(AD_A20_current*10^-12)^2;            % (A^2/Hz)
    h1_A20_voltage = h0_A20_voltage/(2*40*2*log(2));            % h_(-1), in V^2
    h1_A20_current = h0_A20_current/(2*40*2*log(2));            % h_(-1), in A^2
    
        % Thermal EMF noise intensity coefficients
    E_Seebeck = S_ab*DeltaT_ab/10^6;                            % Thermal EMF (in V)
    h0_EMF = E_Seebeck^2;                                       % (V^2/Hz)
    E_Seebeck_scaling = E_Seebeck/Thermal_EMF_baseline;         % Scaling factor
    AD_floor_scaling = E_Seebeck_scaling*AD_floor;              % (nV)
    h1_EMF = AD_floor_scaling^2/(2*log(2));                     % h_(-1), in V^2
    
    %%%%% NOISE (GIVEN IN VOLTS) :: SQUID, A20, AND RESISTORS %%%%%
        % Predefine output array size before loop - for optimal code execution
    In_balArray = double.empty(0,length(tauArray));
    
        % Noise estimation
    for index = 1:length(tauArray)
    
        %%%%% CALCULATE THE EQUIVALENT I1 CURRENT NOISE (IN VOLTS) %%%%%
        I1_AD = I1_ADMatrix(index,2);                           % From 'I_1_Allan_deviations.xlsx' (in A)
        I1_eq_noise = I1_AD*R1;                                 % (in V)

        %%%%% CALCULATE THE NOISE FROM THERMAL EMF (IN VOLTS) %%%%%
        AD_thermal_EMF = sqrt(4*(h0_EMF/(2*tauArray(index))));
        thermal_EMF_noise = sqrt(4*AD_thermal_EMF^2 + sqrt(h1_EMF/(1/tauArray(index)))^2);  % 4x4 connections

        %%%%% CALCULATE SQUID NOISE (IN VOLTS) %%%%%
        h2_SQUID = (tauArray(index)*AD_SQUID_a + AD_SQUID_b)^2*1.5/(tauArray(index)*pi^2);  % h_(-2), in V^2*Hz
        SQUID_noise = sqrt((sqrt(pi^2*tauArray(index)*h2_SQUID/1.5)*(N_ref/N1))^2+...
                        (sqrt(pi^2*tauArray(index)*h2_SQUID/1.5)*(N_ref/(N2+Nb)))^2);
        
        %%%%% CALCULATE A20 NOISE VOLTAGE AND EQUIVALENT CURRENT NOISE (IN VOLTS) %%%%%
        if tauArray(index) <= 40
            A20_noise_voltage = sqrt(h0_A20_voltage/(2*tauArray(index)));
            A20_noise_current = sqrt((sqrt(h0_A20_current/(2*tauArray(index)))*R1)^2+...
                                (sqrt(h0_A20_current/(2*tauArray(index)))*R2)^2);
        else
            A20_noise_voltage = sqrt(h1_A20_voltage*2*log(2));
            A20_noise_current = sqrt((sqrt(h1_A20_current*2*log(2))*R1)^2+...
                                (sqrt(h1_A20_current*2*log(2))*R2)^2);
        end
        
        %%%%% CALCULATE THERMAL NOISE VOLTAGE FROM RESISTORS (IN VOLTS) %%%%%
        thermal_noise_R1 = sqrt(h0_R1/(2*tauArray(index)));
        thermal_noise_R2 = sqrt(h0_R2/(2*tauArray(index)));
        thermal_noise_Resistors = (sqrt(thermal_noise_R1^2 + thermal_noise_R2^2));
        
        %%%%% CALCULATE TOTAL NOISE VOLTAGE AND BALANCE CURRENT NOISE %%%%%
        total_noise_voltage = sqrt(thermal_noise_Resistors^2 + A20_noise_voltage^2 + ...
                                A20_noise_current^2 + SQUID_noise^2 + ...
                                thermal_EMF_noise^2 + I1_eq_noise^2);        
        In_bal = total_noise_voltage*A20_gain*10*F_res/Rb;
            % Evaluate the Allan deviation for an specific tau
            if tauArray(index) == SpecificTau
                FlagMatchSpecificTau = true;                % True if the evaluated time within the loop is equal to
            else                                            % ... the specific tau set by the user through the interface
                FlagMatchSpecificTau = false;
            end

            if FlagSpecificTau_x && FlagMatchSpecificTau == 1
                AD4SpecificTau_x = round(In_bal*10^9, 2);
                AD4SpecificTau = ['O desvio de Allan para o tempo Tau especificado (', num2str(SpecificTauOrig), ...
                    ' s) é de ', num2str(AD4SpecificTau_x), ' nA'];
            end
            % Store balance current noise values in an array:
        In_balArray(index) = In_bal*10^9;                   % Given in nA
    end
    
    %%%%% UNCOMMENT IF YOU WANT TO PLOT %%%%%
    % plot(tauArray, In_balArray, 'Color', [0.4 0 0.4], 'LineWidth', 2)
    % xlabel('Intervalo de medição \tau (s)')
    % xscale('log')
    % xlim([10 100])
    % ylabel('Desvio de Allan (nA)')
    % yscale('log')
    % ylim([1 100])
    % grid on
end