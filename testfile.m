clc,clear
sum1 = 0;
for i=1:20
    [best, best_len]=Improved_GAforTSP();
    sum1 = sum1 + best_len;
end
average1 = sum1 / 20

sum2 = 0;
for i=1:20
    [best, best_len]=GAforTSP();
    sum2 = sum2 + best_len;
end
average2 = sum2 / 20