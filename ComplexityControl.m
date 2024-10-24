% If you use this code for your research work please cite :
% [1] Korosh Mahmoodi, Scott E. Kerick, Piotr J. Franaszczuk, Paolo Grigolini, and Bruce J. West, Complexity Control, 
 

% This code creates two time series, the drive and the driven, using the Manneville map and the equations 3 to 5 in [1].
% Then, for different values of the mu (temporal complexity index) for the drive and the driven systems it plots the 3D graphs 
% of maximum cross correlation of the drive an driven systems (C_{max}), and mu of the driven after connection.

tic
clc;
clear all ;
close all ;

PLOT = 0 ; % If = 1 plots the graphs

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

for bbb = 1 :  21
    bbb
    for ccc = 1 :  20

        mu1(bbb)  =  1 +  bbb * 0.1 ;
        mu2(ccc)  =  1 + ccc * 0.1 ;

        mu01 =  mu1(bbb) ;
        mu02 =  mu2(ccc) ;

        gg1 = 1/(mu01 -1) ;
        gg2 = 1/(mu02 -1) ;

                                        % Drive systtem

        Pi1 = zeros( TimeStep, 1) ;
        tav = zeros( TimeStep, 1) ;

        Start = 2 ;

        AA = 0 ;
        Sign = 1 ;

        while Start <  TimeStep
            AA = AA + 1 ;

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

        AveX = 0 ;
        AveY = 0 ;

        XX = Pi1 ;
        XXTWO = Pi2 ;

        LEG=20 ;
        [c,lags] = xcorr(XX-mean(XX),XXTWO-mean(XXTWO), LEG, 'normalized'); 

        c=c' ;
        Corrr = max(c) ;

        ThreeDCrossCorrelation(bbb, ccc) = Corrr ;

        % Surv1 =  mu(Pi1, ST, EN, PLOT ) ;
        Surv2 =  mu(Pi2, ST, EN, PLOT ) ;

        % Survivea1(bbb, ccc) = 1-Surv1  ;
        Survivea2(bbb, ccc) = 1-Surv2 ;

        PerfectMatch(bbb, ccc) =  mu1(bbb) ;

    end
end


surf(mu2, mu1, ThreeDCrossCorrelation)
xlabel('\mu of driven', 'rotation', 15, 'FontSize', 12), ylabel('\mu of drive', 'rotation', -28, 'FontSize', 12) ; zlabel('C_{max}', 'FontSize', 12) ;
ylim([1 3])
zlim([0 1])

figure

surf(mu2, mu1, Survivea2)
hold on ;

VV = ' [0.9350, 0.4080, 0.3840]' ; 

mesh(mu2, mu1, PerfectMatch,'EdgeColor', VV,'FaceAlpha',0.4,'FaceColor', VV) ;

xlabel('\mu of driven', 'rotation', 15), ylabel('\mu of drive', 'rotation', -28) ; zlabel('\mu of driven after connection') ;
xlim([1 3])
ylim([1 3])
zlim([1 5])

alpha 0.5

hold off ;

toc