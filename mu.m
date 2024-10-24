% If you use this code for your research work please cite :
% [1] Korosh Mahmoodi, Scott E. Kerick, Piotr J. Franaszczuk, Paolo Grigolini, and Bruce J. West, Complexity Control, 

% This function evaluates the slope (= mu - 1) of the survival probability of the taus (time distances between the two consequtive crossings of the origin (y=0) of the time series.)

function slopeSurvi =  mu(data, ST, EN, Plot )
% function [slope, slopeSurvi] =  muR(data, ST, EN, Plot ) % use this if you want both the slope of the wating time probability density functuion (pdf) and survival pdf.

% data is the input
% ST is the start of the linear fit, e.g., 10
% EN is the end of the linear fit, e.g., 60
% if Plot = 1 the function plots the distribution 

% slope is the slope of the waiting time pdf of the taus (= mu)
% slopeSurvi is the slope of the survival probability of the waiting time distribution of the taus (= mu - 1)

aveY = 0 ; %  mean(data) ;
L = length(data) ;
Di = zeros(L, 1) ;

% dichotomizing the data
for cc = 1 : L

    if data(cc) > aveY
        Di(cc) = 1 ;
    else
        Di(cc) = -1 ;
    end

end

Tau = zeros(L , 1) ;

i = 1 ;
j = 1 ;
                   % evaluating time intervals between two consecutive crossings of the origin
while i <= L

    if Di(i) == 1
        while Di(i) == 1

            Tau(j) = Tau(j) + 1 ;
            i = i + 1 ;

            if i >= L
                break
            end

        end

        j = j + 1 ;

    else

        while Di(i) == -1

            Tau(j) = Tau(j) + 1 ;
            i = i + 1 ;

            if i >= L
                break
            end

        end
        j = j + 1 ;
    end

end

% Evaluating the histogram of the taus
XF = Tau(Tau~=0) ;
nbins = floor( max( abs(XF) ) /1 ) ;
nbins(nbins==0) = [] ;
[counts] = hist(XF, nbins) ; 

% discount any entry that is zero
% counts = counts(counts ~= 0) ;

PiS = nbins ;
PY2 = counts' ;

su = cumsum(PY2) ;
PY2 = PY2 ./ su(length(su)) ;

%  plot
de3 = zeros(PiS, 1) ;
DE3 = zeros(PiS, 1) ;

for ttt = 1:PiS
    de3(ttt) = log(ttt)/log(10) ;
    DE3(ttt) = log((PY2(ttt)))/log(10) ;
end

DEafterFirstDecade = DE3(ST :EN) ;
deAfterFirstDecade = de3(ST :EN) ;

FitLine = polyfit(deAfterFirstDecade,DEafterFirstDecade,1) ;
slope = FitLine(1) ;    % this parameter is the scaling

if Plot == 1

    plot(de3(3:length(de3)-3),DE3(3:length(DE3)-3),'.',deAfterFirstDecade,FitLine(1)*deAfterFirstDecade+FitLine(2),'r--','LineWidth',2) ;
    xlabel('log(\tau)'), ylabel('log\Psi(\tau)') ;
    legend(['\mu = ' num2str( -1*slope )],'Location','northeast') ;
    title('Duration') ;

end
%  end of histogram

% Evaluating the survival probability of the taus
SS = length(PY2) ;
Survival = zeros(SS, 1);
for j = 1 : SS

    for jj = 1 : j -1
        Survival(j) =  Survival(j) + PY2(jj) ;
    end
    Survival(j) = 1 - Survival(j) ;
end

%  plot
de3 = zeros(SS, 1) ;
DE3 = zeros(SS, 1) ;

for ttt = 1 : SS
    de3(ttt) = log(ttt)/log(10) ;
    DE3(ttt) = log((Survival(ttt)))/log(10) ;
end

DEafterFirstDecade = DE3(ST :EN) ;
deAfterFirstDecade = de3(ST :EN) ;

% performs a polyfit of degree 1 to find the slope, delta
FitLine = polyfit(deAfterFirstDecade,DEafterFirstDecade,1) ;
slopeSurvi = FitLine(1) ;  % slope of the linear part of the graph  

% plots 

if Plot == 1
    figure
    plot(de3(3:length(de3)-3),DE3(3:length(DE3)-3),'.',deAfterFirstDecade,FitLine(1)*deAfterFirstDecade+FitLine(2),'r--','LineWidth',2) ;
    xlabel('log(\tau)'), ylabel('log\Psi(\tau)');
    legend(['\mu - 1 = ' num2str( -1*slopeSurvi )],'Location','northeast') ;
    title('Survival') ;

end

end
