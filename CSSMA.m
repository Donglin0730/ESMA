
function [Best_pos,Best_score,curve]=CSSMA(pop,Max_iter,lb,ub,dim,fobj)
z=0.03;
if(max(size(ub)) == 1)
    ub = ub.*ones(1,dim);
    lb = lb.*ones(1,dim);
end
X0=initialization(pop,dim,ub,lb);
X = X0;
fitness = zeros(1,pop);
for i = 1:pop
    fitness(i)=fobj(X(i,:));
end
[fitness, index]= sort(fitness);
GBestF = fitness(1);

for i = 1:pop
    X(i,:) = X0(index(i),:);
end
GBestX = X(1,:);
curve=zeros(1,Max_iter);
LW = zeros(pop,dim);
for t = 1: Max_iter
    worstFitness = fitness(end);
    bestFitness = fitness(1);
    S=bestFitness-worstFitness+eps;
    for i = 1: pop
        if i<pop/2
            LW(i,:)= 1+Levy(dim).*log10((bestFitness-fitness(i))/(S)+1);
        else
            LW(i,:)= 1- Levy(dim).*log10((bestFitness-fitness(i))/(S)+1);
        end
    end
    a = atanh(-(t/Max_iter)+1);
    b = 1-t/Max_iter;
    for i=1:pop
        Mean = mean(X);   
        if rand<z
            X(i,:) = (ub-lb)*rand+lb;
        else
            p =tanh(abs(fitness(i)-GBestF));
            vb = unifrnd(-a,a,1,dim);
            vc = unifrnd(-b,b,1,dim);
            A = randi([1,pop]);
            B = randi([1,pop]);
            di = sqrt((X(A,:)-GBestX).^2)-sqrt((Mean-GBestX).^2);
            if fobj(Mean)<fobj(X(A,:))
                flag =1;
            else
                flag = -1;
            end
            for j=1:dim
                r = rand();
                if r<p & di>0
                    if flag == 1
                        X(i,j) = GBestX(j)+ vb(j)*(LW(i,j)*Mean(j)-X(B,j));
                    else
                        X(i,j) = GBestX(j)+ vb(j)*(LW(i,j)*X(A,j)-X(B,j));
                    end
                else
                    X(i,:) = vc.* X(i,:);
                end
            end
        end
    end 
    k = (1+(t/Max_iter)^0.5)^10;
    Temp = (max(GBestX)+ min(GBestX))/2 + (max(GBestX)+ min(GBestX))/(2*k) - GBestX/k;
    Flag4Upper_bound=Temp>ub;
    Flag4Lower_bound=Temp<lb;
    Temp=(Temp.*(~(Flag4Upper_bound+Flag4Lower_bound)))+ub.*Flag4Upper_bound+lb.*Flag4Lower_bound;
    if fobj(Temp)<fobj(GBestX)
        GBestX = Temp;
    end
    for j = 1:pop
        for a = 1: dim
            if(X(j,a)>ub(a))
                X(j,a) =ub(a);
            end
            if(X(j,a)<lb(a))
                X(j,a) =lb(a);
            end
        end
    end
    
    for j=1:pop
        fitness(j)=fobj(X(j,:));
    end
    X_newTemp= X;
    RandIndex = randperm(pop);
    for a = 1:pop/2
        r1 =rand;
        r2 = rand;
        j = pop/2 + a;
        ii = RandIndex(a);
        jj = RandIndex(j);
        MSp=r1.*X_newTemp(ii,:) + (1 - r1).*X_newTemp(jj,:) + sin(2*pi*r1).*(X_newTemp(ii,:)-X_newTemp(jj,:));
        MSp1 = r2.*X_newTemp(jj,:) + (1 - r2).*X_newTemp(ii,:) +  cos(2*pi*r1).*(X_newTemp(jj,:)-X_newTemp(ii,:));
        MSp(MSp>ub) = ub(MSp>ub);
        MSp1(MSp1>ub) = ub(MSp1>ub);
        MSp(MSp<lb) = lb(MSp<lb);
        MSp1(MSp1<lb) = lb(MSp1<lb);
        fit = fobj(MSp);
        fit1 =fobj(MSp1);
        if(fit <fitness(ii))
            X(ii,:) = MSp;
        end
        if(fit1<fitness(jj))
            X(jj,:) = MSp1;
        end
    end
    for a= 1:pop
        RandIndexZ = randperm(dim);
        index1 = RandIndexZ(1);
        index2 = RandIndexZ(2);
        r = rand;
        MSv = X_newTemp(a,:);
        MSv(index1) = r.*X_newTemp(a,index1) + (1 - r).*X_newTemp(a,index2);
        MSv(MSv>ub) = ub(MSv>ub);
        MSv(MSv<lb) = lb(MSv<lb);
        fitv =  fobj(MSv);
        if(fitv < fitness(a))
            X(i,:) = MSv;
        end
    end
    for j=1:pop
        fitness(j)=fobj(X(j,:));
    end
    [fitness, index]= sort(fitness);
    for j = 1:pop
        X(j,:) = X(index(j),:);
    end
    if(fitness(1)<GBestF)
        GBestF=fitness(1);
        GBestX = X(1,:);
    end
    
    curve(t) = GBestF;
end
Best_pos = GBestX;
Best_score = curve(end);
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

