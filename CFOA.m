%%_______________________________________________________________________%%
%             Catch Fish Optimization Algorithm code                      %
%                                                                         %
%  Developed in MATLAB R2022b                                             %
%                                                                         %
%  Author and programmer:                                                 %
%    Catch Fish Optimization Algorithm: A New Human Behavior Algorithm    %
%                   for solving clustering problems                       %
%  e-Mail:                                                                %
%                                                                         %
%  DOI:                                                                   %
%                                                                         %
%_________________________________________________________________________%
%% Catch Fish Optimization Algorithm
function [Best_score,Best_pos,cg_curve]=CFOA(SearchAgents_no,Max_EFs,lb,ub,dim,fobj)
%% ---------------------Initialization parameter--------------------------%
Fisher=initialization(SearchAgents_no,dim,ub,lb);
newFisher=Fisher;
EFs=0;
Best_score=inf;
Best_pos=zeros(1,dim);
cg_curve=zeros(EFs,1);
fit=inf.*ones(SearchAgents_no,1);
newfit=fit;
%% -----------------------Start iteration run-----------------------------%

while EFs<Max_EFs
    for i=1:SearchAgents_no
        Flag4ub = newFisher(i,:)>ub;
        Flag4lb = newFisher(i,:)<lb;
        newFisher(i,:)=(newFisher(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        newfit(i)=fobj(newFisher(i,:));
        if newfit(i)<=fit(i)
            fit(i)=newfit(i);
            Fisher(i,:)=newFisher(i,:);
        end
        if newfit(i)<=Best_score
            Best_pos=Fisher(i,:);
            Best_score=fit(i);
        end
        EFs=EFs+1;
        cg_curve(EFs)=Best_score;
        if EFs>=Max_EFs
            break;
        end
    end
   
    if EFs<Max_EFs/2
            alpha=((1-3*EFs/(2*Max_EFs))^(3*EFs/(2*Max_EFs)));
            p=rand;
            pos=randperm(SearchAgents_no);
            i=1;
            while i <= SearchAgents_no
                per=randi([3 4]);                                          %Randomly determine the size of the group
%% ---------------------Independent search (p < α)------------------------%
                if p < alpha || i+per-1>SearchAgents_no
                    r=randi([1 SearchAgents_no]);
                    while r==i
                        r=randi([1 SearchAgents_no]);
                    end
                    Exp=((fit(pos(i))-fit(pos(r)))/(max(fit)-Best_score));%
                    rs=rand(1,dim)*2-1;
                    rs=norm(Fisher(r,:)-Fisher(i,:)).*rand.*(1-EFs/Max_EFs)*rs/(rs*rs')^0.5;
                    newFisher(pos(i),:)=Fisher(pos(i),:)+(Fisher(pos(r),:)-Fisher(pos(i),:)).*Exp+(abs(Exp)^0.5).*rs;
                    i=i+1;
%% ------------------------Group capture (p ≥ α)--------------------------%
                else
                    aim=sum(fit(pos(i:i+per-1))/sum(fit(pos(i:i+per-1))).*Fisher(pos(i:i+per-1),:));
                    newFisher(pos(i:i+per-1),:)=Fisher(pos(i:i+per-1),:)+rand(per,1).*(aim-Fisher(pos(i:i+per-1),:))+(1-2*EFs/Max_EFs).*(rand(per,dim)*2 -1);
                    i=i+per;
                end
            end
        
            
    else
%% -------------------------Collective capture----------------------------%
        sigma=(2*(1-EFs/Max_EFs)/((1-EFs/Max_EFs)^2+1)).^0.5;
        for i=1:SearchAgents_no
            W=abs(Best_pos-mean(Fisher)).*(randi([1 3])/3)*sigma;
            newFisher(i,:)=Best_pos+(normrnd(0,W,1,dim));
        end
    end

end
end

