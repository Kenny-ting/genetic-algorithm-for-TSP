function [best, best_len]=GAforTSP()
% A是读取的坐标，N*3列，N表示城市数
clc,clear
[N,A] = readfile('att48.tsp'); 
D=zeros(N);
D=get_distance(A);
%种群大小pop=80, 最大世代数gem=100，交叉概率Pc=0.9，变异概率Pm=0.2，C为复制个数
pop=50;
Pc=0.75;
Pm=0.05;
gem=200;
%初始化参数定义部分
P0=zeros(pop,N+1);

P0=zeros(pop,N+1);
for i=1:pop 
    tag=zeros(N,1);
    P0(i,1)=1; P0(i,N+1)=1;
    P0(i,2)=randperm(N-1,1)+1; %随机化P0第二个位置
    for k=2:N-1
        now_pos=P0(i,k);
        tag(now_pos)=1;
        next_pos=0; length=inf;
        for j=2:N %贪婪算法
            if ~tag(j) && D(now_pos,j)<length
                length=D(now_pos,j);
                next_pos = j;
            end
        end
        tag(next_pos)=1;
        P0(i,k+1)=next_pos;
    end
end
P1=P0(:,1:N);
aver = zeros(gem,1);
while gem
    P2=[];                              %P2是子代，先存取父代的值
    P3=[];
    ss=randperm(pop);
    nn=0;
    sss=0;
    for i=1:2:pop-1
        B=P1(ss(i),:);
        C=P1(ss(i+1),:);
        if rand<Pc
            [B,C]=cross(B,C);
            nn=nn+1;
            P2(nn,:)=B;
            P2(nn+1,:)=C;
        end; 
         if rand>1-Pm
            pp=randperm(N);
            [B(pp(1)),B(pp(2))]=exchange(B(pp(1)),B(pp(2)));
            sss=sss+1;
            P3(sss,:)=B;
        end;
        if rand>1-Pm
            ppp=randperm(N);
            [C(ppp(1)),C(ppp(2))]=exchange(C(ppp(1)),C(ppp(2)));
            P3(sss+1,:)=C;
        end;  
    end;
    P=[P1;P2;P3];
    [Fit2,average]=fitness(P,D);
    aver(201-gem)=average;
    [rrr,bbb]=sort(Fit2);
    best_len=rrr(1);
    best=P(bbb(1) ,:);
    P4=P(bbb(1:pop),:);
     P1=P4;
     clear P2;
     clear P3;
     clear P4;
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
end;
for i=1:N
   text(A(i,1)+0.1,A(i,2)+0.1,num2str(i),'fontsize',8);
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
    
function D= get_distance(A) %计算距离矩阵
    [N,col1]=size(A);
	D=zeros(N);
    for i=1:N
        for j=i:N
            D(i,j)=sqrt((A(i,1)-A(j,1)).^2+(A(i,2)-A(j,2)).^2);
            D(j,i)=D(i,j);
        end;
    end;    
    
function [B,C]=cross(B,C)       %部分匹配交叉法
    L=length(B);
    cp1=randperm(L-2);
    if cp1(1)>cp1(2)
        cpmax=cp1(1)+1;
        cpmin=cp1(2)+1;
    else cpmax=cp1(2)+1;
         cpmin=cp1(1)+1;
    end;
    B1=B;
    for i=cpmin:cpmax
        x=find(B==C(i));
        [B(x),B(i)]=exchange(B(x),B(i));
    end
    for i=cpmin:cpmax
        y=find(C==B1(i));
        [C(y),C(i)]=exchange(C(y),C(i));
    end;
    
function [x,y]=exchange(x,y)
        temp=x;
        x=y;
        y=temp;






    