% If you use this code for your research work please cite :
% [1] Korosh Mahmoodi, Scott E. Kerick, Piotr J. Franaszczuk, Paolo Grigolini, and Bruce J. West, Complexity Control, 
 

% This code creates two time series (Pi1 and Pi2), the drive and the driven, using the Manneville map and the equations 3 to 5 in [1]
% The complexity of the drive system changees over time as equation 6 in [1]

% The second part of the code slices the time series Pi1 and Pi2 to pieces and for each piece evaluates
% the slope (mu -1) of the survival probability of the taus (time distances between the two consequtive crossings of the origin (y=0) of the time series)


tic
clc;
clear all ;
close all ;

PLOT = 0 ; % If = 1 plots the survival probabilities in log-log plot

ST = 10 ; % Start of the linear fit
EN = 100 ; % End of the linear fit

TimeStep = 2e7 ; % Lenght of the time series

r0 = 0.05 ; % Strength of the interaction

ThreeDCrossCorrelation = zeros(21, 20, 1) ;

mu1 = zeros(21, 1) ; % Scaling index of the drive system
mu2 = zeros(20, 1) ; % Scaling index of the driven system

Survivea1 = zeros(21, 20, 1) ; % Slope of the Survival probability of the drive system in a log-log plot
Survivea2 = zeros(21, 20, 1) ;

PerfectMatch = zeros(21, 20, 1) ; % The plane of perfect matching 

TT = 1 ; % Parameter of the Manneville map

        mu0Drive = 2.25 ;

        mu02 =  1.5 ; % Scaling parameter of the driven system
          gg2 = 1/(mu02 -1) ;



                                        % Drive systtem

        Pi1 = zeros( TimeStep, 1) ;
        tav = zeros( TimeStep, 1) ;

                                    MuT = zeros(TimeStep, 1) ;

        Start = 2 ;

        AA = 0 ;
        Sign = 1 ;


        while Start <  TimeStep
            AA = AA + 1 ;

                                        mu01 =  mu0Drive + 0.5* cos(Start/3e5) ;  
                                                       MuT(Start) = mu01 ;
                                                gg1 = 1/(mu01 -1) ; 

            r = rand ;
            S = TT * (-1 + 1/(r.^gg1)) ;
            s = round(S) ;

            while s <= 2  ||  s >= 1e4
                r = rand;
                S = TT * (-1 + 1/(r.^gg1)) ;
                s = round(S) ;
            end
            tav(AA) = s ;

            tav(AA) = s ;

            remain =  tav(AA) - 2 ;     

            Part1 = 1 ;
            Part2 = remain ;
            Part3 = 1 ;

            for j = Start  + 1  : Start + Part1

                if j >= TimeStep
                    break
                end
                Pi1(j) =  Pi1(j-1) +  Sign * (1 / Part1 ) ;

            end
            Start = Start + Part1 ;

            for j = Start + 1  : Start + Part2
                if j >= TimeStep
                    break
                end
                Pi1(j) = Pi1(j-1) +  Sign * 0 ;

            end
            Start = Start + Part2 ;

            for j = Start + 1 : Start + Part3

                if j >= TimeStep
                    break
                end
                Pi1(j) =  Pi1(j-1) + Sign * ( - ( 1/ Part3 )  ) ;

            end
            Start = Start + Part3 ;

            if j >= TimeStep
                break
            end

            r = rand ;
            if r < 0.5
                Sign = 1 ;
            else
                Sign = -1 ;
            end

        end

                                        % Driven systtem

        Pi2 = zeros( TimeStep, 1) ;
        tav2 = zeros( TimeStep, 1) ;

        Start = 2 ;
        Sign = 1 ;
        AA = 0 ;
        while Start <  TimeStep
            AA = AA + 1 ;

            r = rand;
            S = TT * (-1 + 1/(r.^gg2)) ;
            s = round(S) ;

            while s <=  2  ||  s >= 1e4
                r = rand ;
                S = TT * (-1 + 1/(r.^gg2)) ;
                s = round(S) ;
            end

            tav2(AA) = s ;

            remain = tav2(AA) - 2 ; 

            Part1 = 1 ;
            Part2 = remain ;
            Part3 = 1 ;

            for j = Start  + 1  : Start + Part1

                if j >= TimeStep
                    break
                end
                Pi2(j) =  Pi2(j-1) + Sign * (1 / Part1 ) +  r0 * ( ( 1 + Pi1(j-1)) * ( 1 - Pi2(j-1)) - ( 1 - Pi1(j-1)) * ( 1 + Pi2(j-1)) ) ;

                if Pi2(j) > 1
                   Pi2(j) = 1 ;
                end

                if Pi2(j) < -1
                   Pi2(j) = -1 ;
                end
            end
            Start = Start + Part1 ;

            for j = Start + 1  : Start + Part2

                if j >= TimeStep
                    break
                end
                Pi2(j) = Pi2(j-1) +  Sign * 0  +  r0 * ( ( 1 + Pi1(j-1)) * ( 1 - Pi2(j-1))  - ( 1 - Pi1(j-1)) * ( 1 + Pi2(j-1)) )  ;

                if Pi2(j) > 1
                   Pi2(j) = 1 ;
                end

                if Pi2(j) < -1
                   Pi2(j) = -1 ;
                end
            end
            Start = Start + Part2 ;

            for j = Start + 1 : Start + Part3

                if j >= TimeStep
                    break
                end
                Pi2(j) = Pi2(j-1) + Sign * ( - ( 1/ Part3 ) ) +  r0 * (  ( 1 + Pi1(j-1)) * ( 1 - Pi2(j-1))  - ( 1 - Pi1(j-1)) * ( 1 + Pi2(j-1)) )  ;

                if Pi2(j) > 1
                   Pi2(j) = 1 ;
                end

                if Pi2(j) < -1
                   Pi2(j) = -1 ;
                end
            end
            Start = Start + Part3 ;

            if j >= TimeStep
                break
            end

            r = rand ;
            if r < 0.5
                Sign = 1 ;
            else
                Sign = -1 ;
            end

        end

        Pi1 = round( Pi1(1:TimeStep, 1), 3 ) ;
        Pi2 = round( Pi2(1:TimeStep, 1), 3) ;



% Evaluating the scaling time series of Pi1 and Pi2


CHANNELS = 2 ;  % We have two time series to analyze

data(:, 1) =  Pi1 ; % First time series (drive)
data(:, 2) =  Pi2 ; % Second time series (driven)

Slice = 3e5 ; % Lenght of the slice of the data used to evaluate its \mu 

nn = floor( ( length(data(:, 1)) - 0)/(Slice/3)  ) - 2 ; % Number of slices

ScaleDur = zeros(nn, CHANNELS) ; % Matrix to save the \mu value of each slice

TimeSlice = zeros(nn, 1) ; % Indexing the slices
for hh11 = 1 : nn
    TimeSlice(hh11) = (Slice/2) + (hh11-1)*(Slice/3) ;
end

PLOT00 = 0 ; % Index to limit the number of plots

for hh11 =  1 : CHANNELS

    for gg11 = 1 : nn

        % each section of data has 2/3 of the previous one
        sta = floor((gg11 -1)*(Slice/3) ) ;
        DaTaa = zeros(Slice, 1) ;
        for yytt = 1 : Slice
            DaTaa(yytt) = data( yytt  + sta, hh11 ) ;
        end

        Surv =  mu(DaTaa,ST, EN, PLOT) ;
        ScaleDur(gg11, hh11) = 1 -Surv  ;

        PLOT00 = PLOT00 + 1 ;
        if PLOT00 > 5
            PLOT = 0 ;
        end

    end

end

% Evaluating the cross-correlation between the two scaling time series
XX = ScaleDur(:, 1) ;
XXTWO = ScaleDur(:, 2) ;

LEG = 20 ;
[c, lags] = xcorr(XX-mean(XX),XXTWO-mean(XXTWO), LEG, 'normalized') ; % xcorr(diff(data),LEG, 'normalized');

c = c' ;
Corrr = max(c)

xIndex = find(c == max(c), 1, 'first') ;
% plot(c)


plot(TimeSlice, ScaleDur,'LineWidth',1.5)

xlabel('Time'), ylabel('\mu') ;
ax = gca ;
ax.FontSize = 10 ;
% xlim([1 3])
ylim([1 3.5])
% legend('Drive: \mu = 2.25 + 0.5*cos(t/3e5)','Driven: \mu = 2.5','FontSize',10) ;
legend('\mu_{D}','\mu_{S}') ;
% title(['c_{max} = ', num2str(Corrr)])
title(['\mu_{D} = 2.25 + 0.5*cos(t/3e5) ', ', \mu_{S} = 1.5'],'FontWeight','Normal') % Change the title according to the paramters used


toc