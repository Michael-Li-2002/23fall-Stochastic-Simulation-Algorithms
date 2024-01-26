%% Project 02 ��� 2100010793
% Exit problem
% ʹ���㷨 1, �޸� option.nsde = @EM_1; ʹ���㷨 2, �޸� option.nsde = @EM_2; 
% ��� L^\infty�����ʱ�䲽�� tau �仯��˫����ͼ
%% Settings
option.checkpoint = checkpoints();                         % ������ֵ��ĵ�
% ʱ�䲽�� tau
option.taulist = [2^(-2);2^(-3);2^(-4);2^(-5);2^(-6);2^(-7);2^(-8);2^(-9);2^(-10);2^(-11);2^(-12)]; 
option.sampletimelist = 1000 * (1:size(option.taulist,1)); % ��������
option.pde = modelpde();                                   % ����pde
option.nsde = @EM_1 ;                                      % ��ֵ��ʽ(Euler-Maruyama ��ȷ�� FirstExitTime �ķ���)

%% Results
NT = size(option.taulist,1); 
NP = size(option.checkpoint,1);
uh = zeros(NT,NP);
Linf_err = zeros(NT,1);

for time = 1:NT
    tau = option.taulist(time);
    sampletime = option.sampletimelist(time);
    for ind = 1:NP
        point = option.checkpoint(ind,:);
        U = 0;
        for i = 1:sampletime
            U = U + option.nsde(option.pde,point,tau);
            if mod(i,200) == 0
                disp(['finish: time = ',num2str(time),' index = ',num2str(ind),' sample ',num2str(i)])
            end
        end
        uh(time,ind) = U/sampletime;
    end
    Linf_err(time) = Linferr(option.pde,option.checkpoint,uh(time,:));
end

%% plot convergence rate
loglog(option.taulist,Linf_err)