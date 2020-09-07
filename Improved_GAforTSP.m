function [best, best_len]=Improved_GAforTSP()
clc,clear
[N,A] = readfile('att48.tsp'); 
D=zeros(N);
D=get_distance(A);
%种群大小pop, 最大世代数gem，交叉概率pc，变异概率pm
pop=50;
Pc=0.9;
Pm=0.02;
gem=200;
%初始化随机数发生器
rand('state',sum(clock)); 

P0=zeros(pop,N+1);
for i=1:pop
    P0(i,2:N)=randperm(N-1)+1; %生成初始群体P0
end;
for i=1:pop
    PP(i,:)=[1 P0(i,2:N) 1];
end; %构造改良圈算法初始矩阵
for k=1:pop
    flag=1;
    while flag
        flag=0;
        for i=1:N-2
            for j=i+2:N
                if D(PP(k,i),PP(k,j))+D(PP(k,i+1),PP(k,j+1))<D(PP(k,i),PP(k,i+1))+D(PP(k,j),PP(k,j+1))
                    PP(k,(i+1):j)=PP(k,j:-1:(i+1));
                    flag=1;
                end
            end
        end
    end
end %PPP为优化过的父代
P1=PP(:,1:N);

[Fit2,average]=fitness(P1,D);
rrr=sort(Fit2);
best_len=rrr(1);
best_fitness=best_len;
aver = zeros(gem,1);
while gem
    P2=[]; %P2是子代，先存取父代的值
    P3=[];
    ss=randperm(pop);
    nn=0;
    sss=0;
    for i=1:2:pop-1
        B=P1(ss(i),:);
        C=P1(ss(i+1),:);
        if rand<get_pc(ss(i),201-gem,Fit2,average,best_fitness)
            [B,C]=cross(B,C);
            nn=nn+1;
            P2(nn,:)=B;
            P2(nn+1,:)=C;
        end 
         if rand>1-get_pm(ss(i),201-gem,Fit2,average,best_fitness)
            pp=randperm(N);
            [B(pp(1)),B(pp(2))]=exchange(B(pp(1)),B(pp(2)));
            sss=sss+1;
            P3(sss,:)=B;
         end
        if rand>1-get_pm(ss(i+1),201-gem,Fit2,average,best_fitness)
            ppp=randperm(N);
            [C(ppp(1)),C(ppp(2))]=exchange(C(ppp(1)),C(ppp(2)));
            P3(sss+1,:)=C;
        end
    end
    P=[P1;P2;P3];
    [Fit2,average]=fitness(P,D);
    aver(201-gem)=average;
    [rrr,bbb]=sort(Fit2);
    best_len=rrr(1);
    best_fitness=best_len;
    best=P(bbb(1) ,:);
    P4=P(bbb(1:round(pop*0.25)),:);
    P5=P(bbb(round(pop*0.25)+1:pop),:);
    P6=[P4;P4;P5];
     P1=P6(1:pop,:);
     clear P2;
     clear P3;
     clear P4;
     clear P5;
     clear P6;
     gem=gem-1;
end
best_len
figure(1)
scatter(A(:,1),A(:,2),15,'o');
hold on;
plot([A(best(1),1),A(best(N),1)],[A(best(1),2),A(best(N),2)],'k');
hold on;
for i=1:N-1
    x0=A(best(i),1);
    x1=A(best(i+1),1);
    y0=A(best(i),2);
    y1=A(best(i+1),2);
    plot([x0, x1],[y0, y1],'k');
    hold on;
end
for i=1:N
   text(A(i,1),A(i,2),num2str(i));
end   
figure(2)
plot(aver);
 
function [Fit,average]=fitness(P,D)        %计算闭合路径值(矩阵）
    pop=size(P,1);
    Fit=zeros(pop,1);        %数值
    [N,coll]=size(D);
    sum = 0;
    for i=1:pop
        for j=1:N-1
          Fit(i,1)=Fit(i,1)+D(P(i,j),P(i,j+1));
        end
          Fit(i,1)=Fit(i,1)+D(P(i,N),P(i,1));
          sum = sum + Fit(i,1);
    end
    average = sum/pop;
    
function [B,C]=cross(B,C)       %部分匹配交叉法
    L=length(B);
    cp1=randperm(L-2);
    if cp1(1)>cp1(2)
        cpmax=cp1(1)+1;
        cpmin=cp1(2)+1;
    else cpmax=cp1(2)+1;
        cpmin=cp1(1)+1;
    end
    B1=B;
    for i=cpmin:cpmax
        x=find(B==C(i));
        [B(x),B(i)]=exchange(B(x),B(i));
    end
    for i=cpmin:cpmax
        y=find(C==B1(i));
        [C(y),C(i)]=exchange(C(y),C(i));
    end
    
function [x,y]=exchange(x,y)
        temp=x;
        x=y;
        y=temp;


function D= get_distance(A) %计算距离矩阵
    [N,col1]=size(A);
	D=zeros(N);
    for i=1:N
        for j=i:N
            D(i,j)=sqrt((A(i,1)-A(j,1)).^2+(A(i,2)-A(j,2)).^2);
            D(j,i)=D(i,j);
        end
    end
    
function pc= get_pc(i, iteration, Fit2, average, best_fitness)
    pc_max = 0;
    pop = 50;
    if iteration <= pop * 0.25
        pc_max = 0.9;
    end
    if pop * 0.25 < iteration <= pop * 0.25 * 3
        pc_max = 0.8;
    end
    if pop * 0.25 * 3 < iteration <= pop
        pc_max = 0.7;
    end

    if Fit2(i) < average
        pc = pc_max;
    else 
        pc = pc_max - (pc_max - 0.6) * (iteration / (2 * pop) +...
            (Fit2(i) - average) / (2 * (best_fitness - average)));
    end
     
function pm = get_pm(i, iteration, Fit2, average, best_fitness)
    pm_min = 0;
    pop = 50;
    if iteration <= pop * 0.25
        pm_min = 0.05;
     end
    if pop * 0.25 < iteration <= pop * 0.25 * 3
        pm_min = 0.1;
    end
    if pop * 0.25 * 3 < iteration <= pop
        pm_min = 0.15;
    end
    
    if Fit2(i) < average
        pm = pm_min;
    else 
        pm = pm_min + (0.2 - pm_min) * ( iteration / (2 * pop) +...
            (Fit2(i) - average) / (2 * (best_fitness - average)));
    end
    


