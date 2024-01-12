%%%%%ANNUAL%%%%%%
AnnualFlux = 0;
AnnualSeries = zeros(3,14);
Moment2 = zeros(1,14);
j = 1;
for i=32:5145 %starts on January 1
    if Day(i)==1
        AnnualFlux = 0;
        Firstday = i;
    end
    AnnualFlux = AnnualFlux + FluxDay(i);
    if Day(i)==366
        AnnualSeries(1,j)=AnnualFlux;
        Lastday = i;
        AnnualSeries(2,j)=mean(Q(Firstday:Lastday));
        AnnualSeries(3,j)=std(Q(Firstday:Lastday));
        j=j+1;
    end
end
for i=1:14
    Moment2(i)=AnnualSeries(2,i)^2 + AnnualSeries(3,i)^2;
end
Pyieldperha = AnnualSeries(1,:)/243000;

% %%%%%%MONTHLY%%%%%%
% MonthlyFlux = 0;
% MonthlySeries = zeros(3,168);
% Moment2 = zeros(1,168);
% Firstday = 32;
% j = 1;
% for i=32:5145 %starts on January 2
%     if Month(i)==Month(i+1)
%         MonthlyFlux = MonthlyFlux + FluxDay(i);
%         Lastday = i+1;
%     else
%         MonthlyFlux = MonthlyFlux + FluxDay(i);
%         MonthlySeries(1,j)=MonthlyFlux;
%         MonthlySeries(2,j)=mean(Q(Firstday:Lastday));
%         MonthlySeries(3,j)=std(Q(Firstday:Lastday));
%         j = j+1;
%         MonthlyFlux = 0;
%         Firstday = i+1;
%     end
% end
% for i=1:168
%     Moment2(i)=MonthlySeries(2,i)^2 + MonthlySeries(3,i)^2;
% end
% Pyieldperha = MonthlySeries(1,:)/243000;

%%%%%%%VISUALIZE AND REGRESSION%%%%%%%%
scatter(Moment2,Pyieldperha);
X = transpose(Moment2);
Y = transpose(Pyieldperha);
regrcoeff = X\Y;