% 初始化
close all; % 关闭所有图形窗口
clearvars; % 清除所有变量

% -- 设置默认/自定义参数
disp('使用默认仿真设置和参数...');
% 设置默认仿真参数
par.runId = 0;              % 仿真ID（用于复现结果）
par.U = 16;                 % 单天线用户数量
par.B = 256;                % 基站天线数量（B>>U）
par.T = 10;                 % 时间槽数量
par.C = 16;                 % 天线簇数量
par.mod = '64QAM';          % 调制类型：'BPSK','QPSK','16QAM','64QAM','8PSK'
par.trials = 1e2;           % 蒙特卡洛试验次数（传输次数）
par.NTPdB_list = -10:2:40;  % 归一化发射功率[dB]值列表
par.rho2 = 1;               % rho^2=1（不应影响结果！）
par.precoder = {'RZF','rKA','SwoR-rKA'}; % 预编码器类型
par.channel = 'Imperfect CSI XL-MIMO Channel'; % 信道模型
par.betaest = 'pilot';      % 信道估计方法
par.save = false;           % 是否保存结果（true/false）
par.plot = true;            % 是否绘制结果（true/false）

% 算法依赖参数
% FD_WF（如果更改场景，请调整）
par.FD_WF.stomp = 0.125*par.C; % 确定正则化步长（与C相关）
% DP_legacy（如果更改场景，请调整）
par.DP_legacy.delta = 0.3; % 拉格朗日乘子缩放因子（1.0）
par.DP_legacy.gamma = 1.0; % 拉格朗日步长（1.0）
par.DP_legacy.maxiter = 2; % 保持为2

% -- 初始化
% 使用runId作为随机种子（确保可复现性）
rng(par.runId);
% 仿真名称（用于保存结果）
par.simName = ['ERR_',num2str(par.U),'x',num2str(par.B), '_C', ...
    num2str(par.C), '_', par.betaest, '_', par.mod, '_', num2str(par.trials),'Trials'];

% 设置格雷映射星座符号（部分符号根据IEEE 802.11选择）
switch (par.mod)
    case 'BPSK',
        par.symbols = [ -1 1 ]; % BPSK符号
    case 'QPSK',
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ]; % QPSK符号
    case '16QAM',
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ]; % 16QAM符号
    case '64QAM',
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ]; % 64QAM符号
    case '8PSK',
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ]; % 8PSK符号
    case '16PSK',
        par.symbols = [ ...
            exp(1i*2*pi*0/16)  ... % 0000
            exp(1i*2*pi*1/16)  ... % 0001
            exp(1i*2*pi*3/16)  ... % 0010
            exp(1i*2*pi*2/16)  ... % 0011
            exp(1i*2*pi*7/16)  ... % 0100
            exp(1i*2*pi*6/16)  ... % 0101
            exp(1i*2*pi*4/16)  ... % 0110
            exp(1i*2*pi*5/16)  ... % 0111
            exp(1i*2*pi*15/16) ... % 1000
            exp(1i*2*pi*14/16) ... % 1001
            exp(1i*2*pi*12/16) ... % 1010
            exp(1i*2*pi*13/16) ... % 1011
            exp(1i*2*pi*8/16)  ... % 1100
            exp(1i*2*pi*9/16)  ... % 1101
            exp(1i*2*pi*11/16) ... % 1110
            exp(1i*2*pi*10/16) ];  % 1111
end

% 计算符号能量
par.Es = mean(abs(par.symbols).^2);

% 每个簇的天线数量
par.S = par.B/par.C;

% 预计算比特标签
par.bps = log2(length(par.symbols)); % 每个符号的比特数
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% 跟踪仿真时间
time_elapsed = 0;
j=sqrt(-1); % 虚数单位
rng(1); % 设置随机种子

%% 仿真参数
M=256; % 总天线数量
K=16; % 用户数量
D_vec=2:2:60; % 每个用户的激活天线数量
rho_dB=10; % 信噪比（dB）
rho=10^(rho_dB/10); % 将dB转换为线性比例
P=1; % 总功率
Normalization=2; % 归一化类型（可能的值：1 - 归一化1，2 - 归一化2）

% 独立信道实现的集合平均
ens=1e2; % 独立信道实现次数

%% 主循环（遍历激活天线数量D）
for ii=1:length(D_vec)
    D=D_vec(ii); % 当前激活天线数量
    %% 最佳情况
    setDbest=[ones(D,K);zeros(M-D,K)]; % 初始化最佳情况下的激活天线索引
    for kk=1:K
        setDbest(:,kk)=circshift(setDbest(:,kk),(kk-1)*D); % 循环移位以获得最佳情况下的激活天线索引
    end
    %% 最差情况
    setDworst=[ones(D,K);zeros(M-D,K)]; % 最差情况下的激活天线索引
    %% 集合平均循环   
        %% 静态仿真
        H=1/sqrt(2)*(randn(M,K)+j*randn(M,K)); % 静态情况下的信道矩阵
        if(Normalization==1)
            Hbest=H.*setDbest*sqrt(M/D); % 最佳情况下的信道矩阵
            Hworst=H.*setDworst*sqrt(M/D); % 最差情况下的信道矩阵
        elseif(Normalization==2)
            Hbest=H.*setDbest; % 最佳情况下的信道矩阵
            Hworst=H.*setDworst; % 最差情况下的信道矩阵
        end
end    
% -- 开始仿真

% - 初始化结果数组（预编码器 x 归一化发射功率）
[res.PER, res.SER, res.BER ] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));
[res.TxPower, res.RxPower, res.TIME] = deal(zeros(length(par.precoder),length(par.NTPdB_list)));

% 计算需要考虑的噪声方差：NTP = rho^2/N0
N0_list = par.rho2*10.^(-par.NTPdB_list/10);

% 试验循环
tic
for t=1:par.trials

    % 生成数据
    for qq=1:par.T
        % 生成随机比特流
        B(:,:,qq) = randi([0 1],par.U,par.bps);
        % 生成发射符号
        Idx(:,qq) = bi2de(B(:,:,qq),'left-msb')+1;
        S(:,qq) = par.symbols(Idx(:,qq)).';
    end

    % 生成掩码
    MaskI = true(par.U,par.T);
    MaskT = true(1,par.T);
    % 在这里可以添加其他估计方法
    switch par.betaest
        case {'pilot'}
            S(:,1) = ones(par.U,1)*sqrt(par.Es); % 在第一个时间槽中发送全1
            MaskI(:,1) = false;
            MaskT(1,1) = false;
        otherwise,
            error('par.betaest未指定')
    end

    % 生成独立同分布高斯信道矩阵和噪声矩阵
    N = sqrt(0.5)*(randn(par.U,par.T)+1i*randn(par.U,par.T));

    % 在这里可以添加自己的信道模型
    switch par.channel
        case 'Imperfect CSI XL-MIMO Channel'
            a=0.1;
            tau=0.3; 
            Phi=a.^toeplitz([0:par.B-1]);
            sqrtmPhi = sqrtm(Phi);
            %H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));
            W=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    Hn=sqrtmPhi*Hbest;
    error_channel=(randn(par.B,par.U)+1i*randn(par.B,par.U))/sqrt(2*par.U);
    H = (sqrt(1-tau^2)*Hn + tau*sqrtmPhi*error_channel)';
        otherwise
            % 完美CSI XL-MIMO信道
            H=Hbest';
    end

    % 算法循环
    for d=1:length(par.precoder)

        % 归一化发射功率循环
        for k=1:length(par.NTPdB_list)

            % 设置噪声方差
            N0 = N0_list(k);

            % 记录波束形成器的运行时间
            starttime = toc;

            % 波束形成器
            switch (par.precoder{d})
               
                case 'rKA',
                    [X,beta]=RKA(par,S,H,N0); % rKA波束形成器
                case 'SwoR-rKA',
                    [X,beta]=RKA2(par,S,H,N0); % SwoR-rKA波束形成器
                case 'RZF',
                    [X,beta]=RZF(par,S,H,N0); % RZF波束形成器
                    case 'ZF',
                    [X,beta]=ZF(par,S,H,N0); % ZF波束形成器
               
                otherwise,
                    error('par.precoder未指定')

            end

            % 记录波束形成仿真时间
            res.TIME(d,k) = res.TIME(d,k) + (toc-starttime);

            % 在噪声信道上传输数据
            HX = H*X;
            Y = HX + sqrt(N0)*N;

            % 提取发射和接收功率
            res.TxPower(d,k) = res.TxPower(d,k) + mean(sum(abs(X(:)).^2))/par.T;
            res.RxPower(d,k) = res.RxPower(d,k) + mean(sum(abs(HX(:)).^2))/par.U/par.T;

            % 用户必须估计beta
            switch par.betaest
                case 'genie', % 完美的beta直接来自波束形成器
                    betaest = ones(par.U,1)*beta;
                case 'pilot', % 知道第一个符号用于训练
                    betaest = real(1./Y(:,1)*sqrt(par.Es)); % ML估计，因为我们对beta没有先验知识
            end

            % 执行估计
            Shat = (betaest*ones(1,par.T)).*Y;

            % 用户端硬输出数据检测
            for qq=1:par.T
                [~,Idxhat(:,qq)] = min(abs(Shat(:,qq)*ones(1,length(par.symbols))-ones(par.U,1)*par.symbols).^2,[],2);
                Bhat(:,:,qq) = par.bits(Idxhat(:,qq),:);
            end

            % -- 计算错误和复杂度指标
            err = (Idx(MaskI)~=Idxhat(MaskI));
            res.PER(d,k) = res.PER(d,k) + any(err(:));
            res.SER(d,k) = res.SER(d,k) + sum(err(:))/par.U/par.T;
            tmpBER = B(:,:,MaskT)~=Bhat(:,:,MaskT);
            res.BER(d,k) = res.BER(d,k) + sum(tmpBER(:))/(par.U*par.bps*sum(MaskT));

        end % 归一化发射功率循环

    end % 算法循环

    % 跟踪仿真时间
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('估计剩余仿真时间：%3.0f分钟。\n',...
            time_elapsed*(par.trials/t-1)/60);
        tic
    end

end % 试验循环

% 归一化结果
res.PER = res.PER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.TIME = res.TIME/par.trials;

% 手动（或视觉）检查预编码器的发射功率是否正确
% （这是许多论文中常见的错误...）
res.TxPower

% -- 保存最终结果（par和res结构）

if par.save
    save([ par.simName '_' num2str(par.runId) ],'par','res');
end

% -- 显示结果（生成Matlab图形）

if par.plot

    % - BER结果
    marker_style = {'kx-','bo:','rs--','mv-.','gp-.','bs--','y*--'};
    h = figure(1);
    for d=1:length(par.precoder)
        semilogy(par.NTPdB_list,res.BER(d,:),marker_style{d},'LineWidth',2);
        if (d==1)
            hold on
        end
    end
    hold off
    grid on
    box on
    xlabel('归一化发射功率 [dB]','FontSize',12,'Interpreter','Latex')
    ylabel('未编码比特错误率 (BER)','FontSize',12,'Interpreter','Latex');
    if length(par.NTPdB_list) > 1
        axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-3 1]);
    end
    legend(par.precoder,'FontSize',12,'Interpreter','Latex','location','northeast')
    set(gca,'FontSize',12);
    if par.save
        % 保存eps图形（彩色，合理边界框）
        print(h,'-loose','-depsc',[ par.simName '_' num2str(par.runId) ])
    end
end
save('ber3',"par","res")

% RKA函数
function [X, beta] = RKA(par,S,H,N0)
numIterations = 100;
updateSchedule = ["power","uniform","aa"];
numRealizations = 20;
% 遍历所有界限
for b = 1:2
    % 运行RKA
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations,H',updateSchedule(2));
end
P=(reshape(V_RKA(:,3,:),[par.B par.U]));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);
% 对信号进行平均缩放
beta = 1/betainv;
end

% RKA2函数
function [X, beta] = RKA2(par,S,H,N0)
numIterations = 100;
updateSchedule = ["power","uniform","aa"];
numRealizations = 20;
% 遍历所有界限
for b = 1:2
    % 运行RKA
    V_RKA = functionRKA(par.B,par.U,1,numRealizations,numIterations,H',updateSchedule(1));
end
P=(reshape(V_RKA(:,3,:),[par.B par.U]));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);
% 对信号进行平均缩放
beta = 1/betainv;
end

%% 最大比传输（MRT）波束形成
function [X, beta, P] = MRT(par,S,H,N0)

% 传输信号
P = H';
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% 对信号进行平均缩放
beta = 1/betainv;

%beta = 1.0*(norm(s,2)^2+N0*par.U)/(s'*H*x); % 这是作弊行为

end

% ZF波束形成
function [X, beta] = ZF(par,S,H,N0)

% 传输信号
P = zfinv(par,H);
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% 对信号进行平均缩放
beta = 1/betainv;

end

% RZF波束形成
function [X, beta] = RZF(par,S,H,N0)

% 传输信号
P =  H'*inv(H*H'+ par.U/par.rho2*eye(par.U));
betainv = sqrt(par.rho2)/sqrt(par.Es*trace(P*P'));
X = betainv*(P*S);

% 对信号进行平均缩放
beta = 1/betainv;

end
