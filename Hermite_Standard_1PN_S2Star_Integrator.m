%定义dst函数，便于后续求解微分方程
function dst=rigrd(t,st,GM,c_AU_year)
x=st(1);
y=st(2);
z=st(3);
vx=st(4);
vy=st(5);
vz=st(6);

%基础几何量计算
r=sqrt(x^2+y^2+z^2);
r_min=1e6/1.495978707e11;%阈值—转换单位为AU
eps_rag=1e-6;%正则化小量
r=max(max(r,r_min),eps_rag);%防止出现奇点截断的数值奇异问题
r3=r^3;
r2=r^2;
v2=vx^2+vy^2+vz^2;
vr=(x*vx+y*vy+z*vz)/r;

%牛顿加速度计算
ax=-GM*x/r3;
ay=-GM*y/r3;
az=-GM*z/r3;
c2=c_AU_year^2;
GM_c2r=GM/(c2*r);
GM_c2r2=GM/(c2*r2);

%1PN修正加速度计算
ax_1PN=ax*(1-4*GM_c2r+v2/c2)+4*GM_c2r2*vr*vx;
ay_1PN=ay*(1-4*GM_c2r+v2/c2)+4*GM_c2r2*vr*vy;
az_1PN=az*(1-4*GM_c2r+v2/c2)+4*GM_c2r2*vr*vz;

%微分方程输出
dxdt=vx;
dydt=vy;
dzdt=vz;
dvxdt=ax_1PN;
dvydt=ay_1PN;
dvzdt=az_1PN;
dst=[dxdt;dydt;dzdt;dvxdt;dvydt;dvzdt];
end

%加加速度推导
function [a,j]=calc_a_j(st,GM,c_AU_year)
x=st(1);
y=st(2);
z=st(3);
vx=st(4);
vy=st(5);
vz=st(6);

%调用rigrd函数获取1PN加速度
dst=rigrd(0,st,GM,c_AU_year);
ax_1PN=dst(4);
ay_1PN=dst(5);
az_1PN=dst(6);
a=[ax_1PN;ay_1PN;az_1PN];

%基础几何量计算
r=sqrt(x^2+y^2+z^2);
r_min=1e6/1.495978707e11;%阈值—转换单位为AU
eps_rag=1e-6;%正则化小量
r=max(max(r,r_min),eps_rag);%防止出现奇点截断的数值奇异问题
r2=r^2;  
r3=r^3;  
r5=r^5;
v2=vx^2+vy^2+vz^2;
r_dot_v=x*vx+y*vy+z*vz;
vr=r_dot_v/r;
c2=c_AU_year^2;

%牛顿加加速度计算
adotx_N=-GM*(vx/r3-3*r_dot_v*x/r5);
adoty_N=-GM*(vy/r3-3*r_dot_v*y/r5);
adotz_N=-GM*(vz/r3-3*r_dot_v*z/r5);
adot_N=[adotx_N;adoty_N;adotz_N];

%定义矢量
r_vec=[x;y;z];
v_vec=[vx;vy;vz];

%加加速度求导中间量
a_dot_r=dot(a,r_vec);
a_dot_v=dot(a,v_vec);
dv2_dt=2*a_dot_v;

%牛顿加速度计算
ax_N=-GM*x/r3;
ay_N=-GM*y/r3;
az_N=-GM*z/r3;
a_N=[ax_N;ay_N;az_N];

%1PN加加速度计算
A=1-4*GM/(c2*r)+v2/c2;
dot_A=4*GM*vr/(c2*r2)+dv2_dt/c2;
d_Part1=adot_N*A+a_N*dot_A;

C=4*GM/c2;
temp_factor=((a_dot_r+v2)-3*vr^2)/r3;
term_deriv_vr_r3_vvec=temp_factor*v_vec;
term_deriv_r_dot_v_a=(r_dot_v/r3)*a;
d_Part2=C*(term_deriv_vr_r3_vvec+term_deriv_r_dot_v_a);

%总加加速度计算
j_total=d_Part1+d_Part2;

%各个方向1PN加加速度计算
adot_x=j_total(1);
adot_y=j_total(2);
adot_z=j_total(3);

%输出加速度和加加速度
j=[adot_x;adot_y;adot_z];
end

%标准4阶Hermite积分器
function[t_out,st_out,step_history]=hermite_standard(t_span,st0,GM_sgrA,c_AU_year,opts,e_val)
t_start=t_span(1);
t_end=t_span(2);

%参数检验
if isfield(opts,'MaxStep')
    MaxStep=opts.MaxStep;
else
    MaxStep=1e3; 
end
if isfield(opts,'RelTol')
    RelTol=opts.RelTol;
else
    RelTol=1e-6;
end
if isfield(opts,'AbsTol')
    AbsTol=opts.AbsTol;
else
    AbsTol=1e-9;
end
if isfield(opts, 'InitialStep')
    h=opts.InitialStep;
else
    %自动估算初始步长：基于轨道周期
    r0=norm(st0(1:3));
    v0=norm(st0(4:6));
    a=1/(2/r0-v0^2/GM_sgrA);
    T=2*pi*sqrt(a^3/GM_sgrA);
    h=T/200;
end

%最小步长（避免过小）
MinStep=1e-6;
h=max(h,MinStep);

%安全系数和最大步长变化
safety=0.9;
max_h_change=2.0;
min_h_change=0.5;

%增加循环保护（防止无限迭代）
max_iter=2e7; 
max_step_num=1e7;
iter_count=0;
step_num=0;

%初始化输出
t_out=zeros(1, max_step_num);
st_out=zeros(max_step_num, 6);

%步长历史记录
step_history=zeros(1, max_step_num);
step_idx=2;
t_out(1)=t_start; 
st_out(1,:)=st0';

%状态初始化
t=t_start;
st=st0;
[a0,j0]=calc_a_j(st,GM_sgrA,c_AU_year);

%循环积分
while t<t_end && iter_count<max_iter && step_num<max_step_num
    iter_count=iter_count +1;
    step_num=step_num +1;

    %限制步长不超过剩余时间
    h=min(h,t_end - t);
    h=min(h,MaxStep);

    %提取当前状态
    r_curr=st(1:3);
    v_curr=st(4:6);
    
    %预测下一步状态行,泰勒展开至3阶，近似下一步状态
    r_pred=r_curr+v_curr*h+0.5*a0*h^2+(1/6)*j0*h^3;
    v_pred=v_curr+a0*h+0.5*j0*h^2;
    st_pred=[r_pred;v_pred];
    
    %计算t+h时刻的真实加速度，用于误差估计
    [a_pred, j_pred] = calc_a_j(st_pred,GM_sgrA,c_AU_year);
    
    %4阶Hermite校正步公式
    v_corr=v_curr+(a0+a_pred)*h/2-(j_pred-j0)*h^2/12;
    r_corr=r_curr+(v_curr+v_corr)*h/2+(a0-a_pred)*h^2/12;
    st_corr=[r_corr;v_corr];
       
    %误差估计（基于预测与校正的差异）
    err_pos=norm(r_corr - r_pred) / (AbsTol + RelTol * norm(r_corr));
    err_vel=norm(v_corr - v_pred) / (AbsTol + RelTol * norm(v_corr));
    err=max(err_pos, err_vel);
    
    %步长控制
    if err > 1

        %步长过大，拒绝
        h=h*max(min_h_change,safety*(1/err)^(1/5));
        h=max(h,MinStep);
        continue;
    else
        % 接受步长
        t=t + h;
        st=st_corr;
        [a0,j0]=calc_a_j(st,GM_sgrA,c_AU_year);
        
        %记录输出
        if step_idx<=max_step_num
            t_out(step_idx)=t;
            st_out(step_idx,:)=st';
            step_history(step_idx)=h;
            step_idx=step_idx+1;
        end
        
        % 计算下一步步长
        h_next=h*min(max_h_change,safety*(1/err)^(1/5));
        h=max(MinStep,min(MaxStep,h_next));
    end
end

%截断步长历史到实际长度
t_out=t_out(1:step_idx-1);
st_out=st_out(1:step_idx-1,:);
step_history=step_history(1:step_idx-1);
end

%定义开普勒轨道要素计算函数（为后续近心点进动角做基础）
function[a,e,i,Omega,omega,f]=kepler_elements(st,GM)%a为半长轴,e为偏心率,i为倾角,Omega为升交点赤经,omega为近心点幅角,f为真近点角
x=st(1);
y=st(2);
z=st(3);
vx=st(4);
vy=st(5);
vz=st(6);

%位置和速度矢量与模长
r=[x;y;z];
v=[vx;vy;vz];
r_mag=norm(r);
v_mag=norm(v);
if r_mag<1e-6||v_mag<1e-6
    a=NaN;
    e=NaN;
    i=NaN;
    Omega=NaN;
    omega=NaN;
    f=NaN;
    return;
end

%角动量矢量
L=cross(r,v);
L_mag=norm(L);
if L_mag<1e-6
    a=inf;
    e=1;
    i=0;
    Omega=0;
    omega=0;
    f=0;
    return;
end

%偏心率矢量
e_vec=(cross(v,L)/GM)-(r/r_mag);
e=norm(e_vec);
e=max(e,1e-6);

%半长轴
a=1/(2/r_mag-v_mag^2/GM);

%轨道倾角:h与z轴夹角
i=acos(clamp(L(3)/L_mag,-1+1e-12,1-1e-12));%留余量防acos(±1)数值警告

%升交点赤经:升交点矢量n与偏心率矢量e_vec的夹角
K=[0;0;1]; 
n=cross(K,L);
n_mag=norm(n);
if n_mag<1e-12
   Omega=0;%赤道轨道，升交点赤经为0
else
    Omega=acos(clamp(n(1)/n_mag,-1+1e-12,1-1e-12));%留余量防数值警告
    if n(2)<0
        Omega=2*pi-Omega;
    end
end

%近心点辐角:升交点矢量n与偏心率矢量e_vec的夹角
if n_mag<1e-12
   omega=acos(clamp(e_vec(1)/e,-1+1e-12,1-1e-12));%留余量防数值警告
   if e_vec(2)<0
       omega=2*pi-omega; 
   end
else
    omega=acos(clamp(dot(n,e_vec)/(n_mag*e),-1+1e-12,1-1e-12));%留余量防数值警告
    if e_vec(3)<0
        omega=2*pi-omega; 
    end
end

%圆轨道近心点幅角无意义
if e<1e-12
    omega=0; 
end 

%真近点角:位置矢量r与偏心率矢量e_vec的夹角
f=acos(clamp(dot(e_vec,r)/(e*r_mag),-1+1e-12,1-1e-12));%留余量防数值警告
if dot(r,v)<0
    f=2*pi-f;
end

%圆轨道真近点角无意义
if e<1e-12
    f=0; 
end 
end

%辅助函数:数值截断（防止acos输入超出[-1,1]）
function val=clamp(val,min_val,max_val)
val=max(val, min_val);
val=min(val, max_val);
end

%定义物理常数（国际单位）
AU=1.495978707e11;
G=6.67430e-11;
M_sun=1.9891e30;
year=365.25*24*3600;
c=299792458;

%国际单位换算成天文单位AU-year
c_AU_year=c*year/AU;
M1=4.297e6*M_sun;
M2=14*M_sun;
GM_sgrA=G*(M1+M2)*(year^2)/(AU^3);
GM=G*(M1+M2);
GM_c2_m=GM/c^2;
GM_c2_AU=GM_c2_m/AU;
c2_AU_year=c_AU_year^2;
GM_sgrA_c2_AU_year=GM_sgrA/c2_AU_year;

%三维轨道角度的定义
i=deg2rad(133.818);
Omega=deg2rad(227.85);
omega=deg2rad(66.13);

%旋转矩阵
Rx_i=[1 0 0;0 cos(i) sin(i);0 -sin(i) cos(i)];
Rz_Omega=[cos(Omega) -sin(Omega) 0;sin(Omega) cos(Omega) 0;0 0 1];
Rz_omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
R_total=Rz_omega*Rx_i*Rz_Omega;

%定义轨道常数与函数
a_AU=971;%谐和坐标系(由观测天球坐标系值乘以KL=1.00105 得到)
e_array=[0 0.2 0.4 0.5 0.6 0.7 0.8 0.8842 0.9];
sample_rate=5;
opts=struct();
opts.RelTol=1e-6;
opts.AbsTol=1e-9;
opts.MaxStep=0;

%初始化S2星数据存储变量
s2_hermite_data=struct();
s2_hermite_data.t=[];%时间序列
s2_hermite_data.h=[];%步长序列
s2_hermite_data.r=[];%轨道半径序列
s2_hermite_data.f=[];%真近点角序列

%预分配误差数组
num_e=length(e_array);
wuca_yuan_array=zeros(1,num_e);   
wuca_jing_array=zeros(1,num_e);    
wuca_zhouqi_array=zeros(1,num_e); 
wuca_E_array=zeros(1,num_e);     
wuca_L_array=zeros(1,num_e);   
L_jingdong_lilun_array=zeros(1,num_e); 
wuca_jingdong_array=NaN(1,num_e);
wuca_jingdong_obs_array=zeros(1,num_e);
time_cost_e=zeros(1,num_e);
step_num_e=zeros(1,num_e);
pianyi_E_array=zeros(1,num_e);       
pianyi_L_array=zeros(1,num_e);    
L_jiaodong_deg_array=zeros(1,num_e); 

%初始化S2星收敛数据存储
T_nodes=[1,3,5,7,10,50,100,500,1000];%周期节点
num_nodes=length(T_nodes);

%收敛数据存储（针对e=0.8842）
conv_steps=[];%累积积分步数
conv_prec_err=[];%近心点进动角相对误差
conv_energy_err=[];%机械能相对误差

%临时存储每个周期节点的进动角和机械能相对误差 
node_prec_data=zeros(length(T_nodes),2);%[进动角相对误差，机械能相对误差]
node_step_data=zeros(length(T_nodes),1);%步数

%预分配S2星近心点数据存储变量
s2_peri_data=struct();

%存储1PN史瓦西相对论参考轨道（理论值，含近心点进动）
s2_peri_data.r_peri=0; %近心点距离（AU）
s2_peri_data.x_ref=[]; %1PN相对论参考轨道X坐标
s2_peri_data.y_ref=[]; %1PN相对论参考轨道Y坐标

%存储Hermite计算数据
s2_peri_data.hermite_x=[]; %Hermite轨道X坐标（近心点区域）
s2_peri_data.hermite_y=[]; %Hermite轨道Y坐标（近心点区域）
s2_peri_data.hermite_steps=[]; %Hermite积分步长
s2_peri_data.hermite_sample_idx=[]; %采样点索引（用于绘图标记）

%圆轨道的特殊处理（误差设为0）
wuca_yuan_array(1)=0;
wuca_jing_array(1)=0;
wuca_zhouqi_array(1)=0;
wuca_E_array(1)=0;
wuca_L_array(1)=0;

%三维圆轨道求解绘图
e_val=0;
r0_AU=a_AU;
v0_AU_year=sqrt(GM_sgrA/r0_AU);
r0_yuanchu=[r0_AU;0;0];
v0_yuanchu=[0;v0_AU_year;0];
r0_yuanchu3d=R_total*r0_yuanchu;
v0_yuanchu3d=R_total*v0_yuanchu;
yuanchu=[r0_yuanchu3d;v0_yuanchu3d];

%1PN参考周期
T_Newton=2*pi*sqrt(a_AU^3/GM_sgrA);
beta_1PN=(15*GM_sgrA)./(8*c2_AU_year*a_AU);
gamma_1PN=(3*GM_sgrA)./(c2_AU_year*a_AU.*(1-e_val.^2));
T_1PN=T_Newton.*(1+beta_1PN+gamma_1PN);
T=T_1PN;
opts.MaxStep=T/200;
t_span=[0,1000*T];

%计算并运行积分器
[t0,s0,step_history0]=hermite_standard(t_span,yuanchu,GM_sgrA,c_AU_year,opts,e_val);
x0=s0(1:sample_rate:end,1);
y0=s0(1:sample_rate:end,2);
z0=s0(1:sample_rate:end,3);
figure(1)
hold on; 
grid on
plot3(x0,y0,z0,'r-','LineWidth',1.5,'DisplayName','圆轨道 (e=0)');
plot3(0,0,0,'ro','MarkerSize',8,'MarkerFaceColor','r','DisplayName','Sgr A*黑洞');

%不同偏心率的三维椭圆轨道求解与绘图
color_map=jet(length(e_array)-1);
for idx = 2:length(e_array)
    e_val=e_array(idx);
    r_min=a_AU*(1-e_val);
    r_max=a_AU*(1+e_val);
    v_tuo=sqrt(GM_sgrA*(1+e_val)/r_min);
    r_tuoyuanchu=[r_min;0;0];
    v_tuoyuanchu=[0;v_tuo;0];
    r_tuoyuanchu3d=R_total*r_tuoyuanchu;
    v_tuoyuanchu3d=R_total*v_tuoyuanchu;
    tuoyuanchu=[r_tuoyuanchu3d;v_tuoyuanchu3d];

    %1PN参考周期（含近心点回归）
    T_Newton=2*pi*sqrt(a_AU^3/GM_sgrA);
    beta_1PN=(15*GM_sgrA)./(8*c2_AU_year*a_AU);
    gamma_1PN=(3*GM_sgrA)./(c2_AU_year*a_AU.*(1-e_val.^2));
    T_1PN=T_Newton.*(1+beta_1PN+gamma_1PN);
    T=T_1PN;
    t_span=[0,1000*T];
    opts.MaxStep=T/200;
    if e_val>=0.9
    opts.MaxStep=T/500;
    end

    %S2星实测数据对标
    delta_phi_1000T_OBS=NaN;
    if abs(e_val-0.8842)<1e-4
    delta_phi_1000T_OBS=2.017;     
    end

    %高偏心率减小初始步长,提高积分稳定性
    if e_val>=0.9
        InitialStep=T/5000;
    elseif e_val>=0.88
        InitialStep=T/2000;
    elseif e_val>=0.8
        InitialStep=T/1000;
    else
        InitialStep=T/200;
    end 
    opts.InitialStep=InitialStep;

    %最小步长保护
    MinStep=max(T/5000, 1e-6);
    new_InitialStep=max(opts.InitialStep,MinStep);
    opts.InitialStep=new_InitialStep;

    %计算并运行积分器
    tic;
    [t1,s1,step_history]=hermite_standard(t_span,tuoyuanchu,GM_sgrA,c_AU_year,opts,e_val);
    time_cost_e(idx)=toc;%记录耗时
    step_num_e(idx)=length(t1);%积分步数=时间序列长度

    %采样数据,仅第一次运行绘图
    sample_rate_temp = max(2, min(sample_rate, floor(length(t1)/500)));
    x1=s1(1:sample_rate_temp:end,1);
    y1=s1(1:sample_rate_temp:end,2);
    z1=s1(1:sample_rate_temp:end,3);
    vx=s1(1:sample_rate_temp:end,4);
    vy=s1(1:sample_rate_temp:end,5);
    vz=s1(1:sample_rate_temp:end,6);
    x=s1(:,1);
    y=s1(:,2);
    z=s1(:,3);
    vx1=s1(:,4);
    vy1=s1(:,5);
    vz1=s1(:,6);
    plot3(x1,y1,z1,'LineWidth',1.2,...
        'Color', color_map(idx-1, :), 'DisplayName', ['椭圆轨道 (e=', num2str(e_val), ')']);

        %提取各周期节点的积分步数
        if abs(e_val-0.8842)<1e-4
            for n=1:length(T_nodes)
                T_target=T_nodes(n)*T;
                idx_target=find(t1<=T_target,1,'last');%找到节点对应索引
                if ~isempty(idx_target)
                    node_step_data(n)=idx_target;%存储积分步数
                else
                    node_step_data(n)=NaN;
                end
            end

            %提取S2星（e=0.8842）近心点区域数据 
            r_all=sqrt(x.^2 + y.^2 + z.^2);%计算轨道半径，定位近心点位置 
            [r_peri_val,peri_idx]=min(r_all);%近心点距离+索引 
            s2_peri_data.r_peri=r_peri_val; 
            peri_x=x(peri_idx); 
            peri_y=y(peri_idx); 
            peri_z=z(peri_idx); 

            %定义近心点放大范围（r_peri的1%） 
            delta_r=0.01*r_peri_val;%放大窗口大小 
            peri_win_idx=find((x >= peri_x-delta_r) & (x <= peri_x+delta_r) & ... 
                               (y >= peri_y-delta_r) & (y <= peri_y+delta_r) & ... 
                               (z >= peri_z-delta_r) & (z <= peri_z+delta_r)); 

            % 处理空数据情况 
            if isempty(peri_win_idx) 
                
                %扩大搜索范围 
                delta_r=0.02*r_peri_val; 
                peri_win_idx=find((x >= peri_x-delta_r) & (x <= peri_x+delta_r) & ... 
                                   (y >= peri_y-delta_r) & (y <= peri_y+delta_r) & ... 
                                   (z >= peri_z-delta_r) & (z <= peri_z+delta_r)); 
                
                if isempty(peri_win_idx) 
                    warning('无法找到近心点区域的数据点'); 
                    s2_peri_data.hermite_x=[]; 
                    s2_peri_data.hermite_y=[]; 
                    s2_peri_data.hermite_steps=[]; 
                    s2_peri_data.hermite_sample_idx=[]; 
                else 
                    %提取近心点区域的hermite轨道坐标 
                    s2_peri_data.hermite_x=x(peri_win_idx); 
                    s2_peri_data.hermite_y=y(peri_win_idx); 
                    
                    %提取对应步长（确保长度匹配） 
                    if ~isempty(step_history) && length(step_history) >= peri_win_idx(end)
                         s2_peri_data.hermite_steps = step_history(peri_win_idx);
                    else
                        s2_peri_data.hermite_steps = diff(t1([peri_win_idx, peri_win_idx(end)+1]));
                        warning('步长历史记录长度不足，已通过diff(t1)计算近心点区域步长'); 
                    end 
                    
                    %标记采样点（每隔5步标记一个) 
                    s2_peri_data.hermite_sample_idx=1:5:length(peri_win_idx); 
                end 
            else 
                %提取近心点区域的hermite轨道坐标 
                s2_peri_data.hermite_x=x(peri_win_idx); 
                s2_peri_data.hermite_y=y(peri_win_idx); 
                
                %提取对应步长（确保长度匹配） 
                if ~isempty(step_history) && max(peri_win_idx) <= length(step_history)
                    s2_peri_data.hermite_steps = step_history(peri_win_idx);
                else
                    s2_peri_data.hermite_steps = diff(t1([peri_win_idx, peri_win_idx(end)+1]));
                    warning('步长历史记录长度不足/索引越界，已通过diff(t1)计算近心点区域步长'); 
                end
                
                %标记采样点（每隔5步标记一个) 
                s2_peri_data.hermite_sample_idx=1:5:length(peri_win_idx); 
            end

            %计算1PN史瓦西相对论解析轨道（理论值，含近心点进动） 
            theta=linspace(-pi/4, pi/4, 1000);%近心点附近极角 
            delta_phi_1PN_rad=6*pi*GM_sgrA/(c2_AU_year*a_AU*(1-e_val^2));%每圈进动角（弧度）
            phi_1PN=theta +(delta_phi_1PN_rad/(2*pi))*theta;
            r_1PNr=a_AU*(1-e_val^2)./(1 + e_val*cos(phi_1PN));%1PN相对论轨道极径（含进动）
            x_1PN=r_1PNr .* cos(theta); 
            y_1PN=r_1PNr .* sin(theta); 
            z_1PN=zeros(size(x_1PN));

            %旋转到真实轨道平面 
            kepler_3d=R_total * [x_1PN; y_1PN;z_1PN]; 

            %提取近心点区域的参考轨道 
            ref_win_idx=find((kepler_3d(1,:) >= peri_x - delta_r) & (kepler_3d(1,:) <= peri_x + delta_r) & ... 
                               (kepler_3d(2,:) >= peri_y - delta_r) & (kepler_3d(2,:) <= peri_y + delta_r)); 

            %处理空数据情况 
            if isempty(ref_win_idx) 
                warning('无法找到参考轨道在近心点区域的数据点'); 
                s2_peri_data.x_ref=[]; 
                s2_peri_data.y_ref=[]; 
            else 
                s2_peri_data.x_ref=kepler_3d(1, ref_win_idx); 
                s2_peri_data.y_ref=kepler_3d(2, ref_win_idx); 
            end 
            
            %提取S2星（e=0.8842）的t/h/r/f完整数据
            s2_hermite_data.t=t1; %Hermite积分器输出的完整时间序列

             %步长序列h
             if ~isempty(step_history) && length(step_history) == length(t1)-1
                 s2_hermite_data.h = step_history(:);
                 fprintf('步长生成：t1长度=%d，h长度=%d（使用积分器原生step_history）\n',length(t1),length(s2_hermite_data.h));
             elseif length(t1) >= 2
                 s2_hermite_data.h = diff(t1);
                 fprintf('步长生成：t1长度=%d，h长度=%d（通过diff(t1)计算）\n',length(t1),length(s2_hermite_data.h));
             else
                 s2_hermite_data.h = []';
                 warning('t1长度不足（<2），无法生成步长');
             end

             %轨道半径序列r
             s2_hermite_data.r=sqrt(x.^2 + y.^2 + z.^2);  % 逐点计算轨道半径
            
             %真近点角序列f
             f_series=zeros(length(t1), 1);%预分配真近点角数组
             for k=1:length(t1)
                 st_k = [x(k); y(k); z(k); vx1(k); vy1(k); vz1(k)];
                 [~,~,~,~,~,f_rad]=kepler_elements(st_k, GM_sgrA);
                 f_series(k)=f_rad;
             end
             s2_hermite_data.f=f_series;%真近点角序列

             %保存S2星Hermite数据到mat文件
             save('s2_hermite_stepdata.mat', 's2_hermite_data'); 
        end 

    %重置节点数据
    node_prec_data=zeros(length(T_nodes),2);
    node_step_data=zeros(length(T_nodes),1);

    %计算远心点/近心点位置与相对误差
    r=sqrt(x.^2+y.^2+z.^2);
    if e_val>=0.9
        r_smooth=smoothdata(r,'movmean',3);
    elseif e_val>=0.88
        r_smooth=smoothdata(r,'movmean',5);
    else
        r_smooth=smoothdata(r,'movmean',10);
    end

    r_yuan=r_max;
    r_jing=r_min;
    wuca_yuan=0;
    wuca_jing=0;
    locs_yuan=[];
    bizhi_jing_yuan=0;

    if e_val>=0.88
        min_peak_dist=max(3,round(length(r_smooth)/800));
    else
        min_peak_dist=max(5,round(length(r_smooth)/500));%减小最小峰距,提高检测率
    end
    min_peak_dist=min(min_peak_dist,length(r_smooth)-1);

    if length(r_smooth)>=10
        is_max=islocalmax(r_smooth,'MinSeparation',min_peak_dist);%远心点（半径极大值）用 islocalmax
        locs_yuan=find(is_max);
        is_min=islocalmin(r_smooth,'MinSeparation',min_peak_dist);%近心点（半径极小值）用 islocalmin
        locs_jing=find(is_min);

    %远心点计算
    if ~isempty(locs_yuan)
        r_yuan=r(locs_yuan(1));
        wuca_yuan=abs(r_yuan-r_max)/r_max*100;          
    end

    %近心点计算
    if ~isempty(locs_jing)
        r_jing=r(locs_jing(1));
        wuca_jing=abs(r_jing - r_min)/r_min*100;
    end 

    %误差比值（防除0）
    if wuca_yuan>1e-6 && wuca_jing>1e-6 
            bizhi_jing_yuan=wuca_jing/wuca_yuan;
    elseif wuca_yuan<=1e-6 && wuca_jing>1e-6
            bizhi_jing_yuan=0;
    else
        bizhi_jing_yuan=0;
    end
    end

    %存储远/近心点误差
    wuca_yuan_array(idx)=wuca_yuan;
    wuca_jing_array(idx)=wuca_jing;

    %计算轨道周期与相对误差
    T_tuo=NaN;
    wuca_zhouqi=NaN;
    r_original=sqrt(x.^2+y.^2+z.^2);

    %平滑前判断r_original长度
    r_smooth=smoothdata(r_original, 'movmean', 5);

    %使用 islocalmin 检测局部最小值（近心点）
    is_min=islocalmin(r_smooth, 'MinSeparation', 50);
    locs_jing=find(is_min);

    %如果检测到的近心点过少，放宽间隔
    if length(locs_jing) < 10
        is_min=islocalmin(r_smooth, 'MinSeparation', 20);
        locs_jing=find(is_min);
    end
    
    if length(locs_jing)>=3
        T_tuo_array=diff(t1(locs_jing));

        %去除首尾各一个近心点，避免初始瞬态影响
        if length(T_tuo_array) > 2
           T_mean=median(T_tuo_array);
        else
            T_mean=mean(T_tuo_array);
        end
        T_std=std(T_tuo_array);
        T_tuo_array=T_tuo_array(abs(T_tuo_array-T_mean)<=3*T_std);

        if ~isempty(T_tuo_array)
            T_tuo=mean(T_tuo_array);
            wuca_zhouqi=abs(T_tuo-T)/T*100;
        else
            T_tuo=mean(diff(t1(locs_jing)));
            wuca_zhouqi=abs(T_tuo-T)/T*100;
        end

    elseif length(locs_jing)==2%若近心点数量为2，直接计算单周期
        T_tuo=diff(t1(locs_jing));
        wuca_zhouqi=abs(T_tuo-T)/T*100;
    end
    wuca_zhouqi_array(idx)=wuca_zhouqi;

    %计算机械能相对误差和偏移量（守恒分析）
    if ~isempty(vx1) && ~isempty(vy1) && ~isempty(vz1) && ~isempty(r)
        r=sqrt(x.^2+y.^2+z.^2);
        v2=vx1.^2+vy1.^2+vz1.^2;
        vr=(x.*vx1+y.*vy1+z.*vz1)./r;
        vr2=vr.^2;
        E=0.5*v2-GM_sgrA./r;
        delta_E_1PN=(GM_sgrA_c2_AU_year./r).*((1/2)*GM_sgrA./r+(3/2)*v2-0.5*vr2);
        E_1PN=E+delta_E_1PN;
        E_chu=E_1PN(1);
        E_end=E_1PN(end);
        wuca_E=abs(max(E_1PN)-min(E_1PN))/abs(E_chu)*100;%守恒波动误差
        pianyi_E=abs(E_end-E_chu)/abs(E_chu)*100;%终点偏移误差

        %对S2星提取各周期节点的机械能相对误差
        if abs(e_val-0.8842)<1e-4
            for n=1:length(T_nodes)

                %重新计算当前运行的节点索引
                T_target=T_nodes(n)*T;
                idx_target=find(t1<=T_target,1,'last');

                if ~isempty(idx_target) && idx_target<=length(E_1PN)

                    % 存储当前运行的节点步数
                    node_step_data(n)=idx_target;
                    E_node=E_1PN(1:idx_target);
                    node_prec_data(n,2)=abs(max(E_node)-min(E_node))/abs(E_node(1))*100;
                else
                    node_prec_data(n,2)=NaN;
                    node_step_data(n)=NaN;
                end
            end
        end

    else
        wuca_E=NaN;
        pianyi_E=NaN;
    end
    wuca_E_array(idx)=wuca_E;
    pianyi_E_array(idx)=pianyi_E;

    %计算角动量相对误差和偏移量（守恒分析）
    if ~isempty(x) && ~isempty(y) && ~isempty(z) && ~isempty(vx1) && ~isempty(vy1) && ~isempty(vz1)
        Lx=y.*vz1-z.*vy1;
        Ly=z.*vx1-x.*vz1;
        Lz=x.*vy1-y.*vx1;
        L=sqrt(Lx.^2+Ly.^2+Lz.^2); 
        L_chu=L(1);
        L_end=L(end);
        wuca_L=abs(max(L)-min(L))/abs(L_chu)*100;%角动量大小波动误差
        pianyi_L=abs(L_end-L_chu)/L_chu*100;%角动量大小偏移误差
    else
        wuca_L=NaN;
        pianyi_L=NaN;
    end
    wuca_L_array(idx)=wuca_L;
    pianyi_L_array(idx)=pianyi_L;
    
    %计算角动量方向偏移角
    if ~isempty(Lx) && ~isempty(Ly) && ~isempty(Lz)
        Lx1=Lx(1);
        Ly1=Ly(1);
        Lz1=Lz(1);
        Lx_end=Lx(end);
        Ly_end=Ly(end);
        Lz_end=Lz(end);
        dot_L=Lx1*Lx_end+Ly1*Ly_end+Lz1*Lz_end;
        theta_cos=dot_L/(L_chu*L_end);
        theta_cos=max(min(theta_cos,1-1e-12),-1+1e-12);%数值稳定性修正,留余量防数值警告
        L_jiaodong_rad=acos(theta_cos);
        L_jiaodong_deg=L_jiaodong_rad*180/pi; %转换为角度
    else
        L_jiaodong_deg=NaN;
    end
    L_jiaodong_deg_array(idx)=L_jiaodong_deg;

    %计算近心点进动角
    L_jingdong=NaN;
    wuca_jingdong=NaN;
    L_jingdong_lilun_1000T_deg=NaN;

    %计算轨道半径序列，用于近心点极值点识别
    r_all=sqrt(x.^2+y.^2+z.^2);

    %动态设置近心点识别参数
    if e_val>=0.9
        min_peak_dist_jd=max(2,round(length(r_all)/2000));%极高偏心率:更小峰距
    else
        min_peak_dist_jd=max(3,round(length(r_all)/500));
    end
    min_peak_dist_jd=min(min_peak_dist_jd,length(r_all)-1);
    
    %使用islocalmin检测近心点
    r_smooth=smoothdata(r_all,'movmean',5);
    is_min=islocalmin(r_all,'MinSeparation',min_peak_dist_jd);
    locs_jing=find(is_min);

    %提取近心点时刻和对应的状态量
    if length(locs_jing)>=3

        %插值参数
        N=3;%每侧点数
        t_jing_exact=zeros(length(locs_jing),1);%近心点时刻序列
        omega_jing_exact=zeros(length(locs_jing),1);%近心点幅角序列

        %只计算近心点时刻的开普勒要素
        for k=1:length(locs_jing)
            peakIdx=locs_jing(k);%当前近心点在原序列中的索引

            %确定插值区间
            start_idx=max(1,peakIdx-N);
            end_idx=min(length(t1),peakIdx+N);

            %局部时间、半径序列
            t_local=t1(start_idx:end_idx);
            r_local=r_all(start_idx:end_idx);

            %三次样条插值得到精细时间-半径曲线
            t_fine=linspace(t_local(1),t_local(end),100);
            r_fine=interp1(t_local,r_local,t_fine,'spline');

            %找到半径极小值点（精确近心点时刻）
            [~,min_idx]=min(r_fine);
            t_peri=t_fine(min_idx);
            t_jing_exact(k)=t_peri;

            %位置（用样条插值保持光滑）
            x_local=x(start_idx:end_idx);
            y_local=y(start_idx:end_idx);
            z_local=z(start_idx:end_idx);
            x_peri=interp1(t_local,x_local,t_peri,'spline');
            y_peri=interp1(t_local,y_local,t_peri,'spline');
            z_peri=interp1(t_local,z_local,t_peri,'spline');

            %速度（用样条插值保持光滑）
            vx_local=vx1(start_idx:end_idx);
            vy_local=vy1(start_idx:end_idx);
            vz_local=vz1(start_idx:end_idx);
            vx_peri=interp1(t_local,vx_local,t_peri,'spline');
            vy_peri=interp1(t_local,vy_local,t_peri,'spline');
            vz_peri=interp1(t_local,vz_local,t_peri,'spline');

            %组成状态向量
            st_peri=[x_peri;y_peri;z_peri;vx_peri;vy_peri;vz_peri];

            % 计算近心点幅角
            [~,~,~,~,omega_rad,~]=kepler_elements(st_peri,GM_sgrA);
            if ~isnan(omega_rad)
                omega_jing_exact(k)=rad2deg(omega_rad);
            else
                omega_jing_exact(k)=0;%容错处理
            end
        end

        %将插值结果赋给原变量
        t_jing=t_jing_exact;
        omega_jing=omega_jing_exact;

        %对幅角解缠绕
        omega_jing_unwrap=zeros(size(omega_jing));
        omega_jing_unwrap(1)=omega_jing(1);

        for k=2:length(omega_jing)

            %计算当前点和上一个点的差值
            delta_candidate=omega_jing(k)-omega_jing_unwrap(k-1);
            
            %如果差值的绝对值大于180度，就减去或加上360度
            if delta_candidate>180
                omega_jing_unwrap(k)=omega_jing(k)-360;
            elseif delta_candidate<-180
                omega_jing_unwrap(k)=omega_jing(k)+360;
            else
                omega_jing_unwrap(k)=omega_jing(k);
            end
            
            %确保是连续递增的
            if omega_jing_unwrap(k)<omega_jing_unwrap(k-1)
                omega_jing_unwrap(k)=omega_jing_unwrap(k)+360;
            end
        end

        %有效数据筛选（剔除前1个近心点，避免暂态）
        valid_idx_jd=~isnan(omega_jing_unwrap) & ~isnan(t_jing);
        valid_idx_jd(1)=false;%剔除第一个暂态点
        t_valid=t_jing(valid_idx_jd);
        omega_valid=omega_jing_unwrap(valid_idx_jd);

        %一阶线性拟合
     if length(t_valid)>=2
        p=polyfit(t_valid,omega_valid,1);
        omega_dot=p(1);%进动率

        %计算1000周期总进动角
        t_total_1000T=1000*T;
        L_jingdong=abs(omega_dot*t_total_1000T);%数值进动角（度）

        %计算理论近心点进动角
        L_jingdong_lilun_1PN_rad=6*pi*GM_sgrA/(c2_AU_year*a_AU*(1-e_val^2));
        L_jingdong_lilun_1PN_deg=rad2deg(L_jingdong_lilun_1PN_rad*1000);
        L_jingdong_lilun_1000T_deg=L_jingdong_lilun_1PN_deg;

        %近心点进动角相对误差
        if L_jingdong_lilun_1000T_deg>1e-6
            wuca_jingdong=abs(L_jingdong-L_jingdong_lilun_1000T_deg)/L_jingdong_lilun_1000T_deg*100;
        end
        
        %S2星进动角实测对标相对误差
        if abs(e_val-0.8842)<1e-4 && ~isnan(delta_phi_1000T_OBS)
            wuca_jingdong_obs_array(idx)=abs(L_jingdong-delta_phi_1000T_OBS)/delta_phi_1000T_OBS*100;
        end

        %对S2星提取各周期节点的近心点进动角相对误差
        if abs(e_val-0.8842)<=1e-4
            for n=1:length(T_nodes)
                T_target=T_nodes(n)*T;

                %找到该周期节点前的所有近心点
                idx_jing_target=find(t_jing<=T_target,1,'last');
                
                if ~isempty(idx_jing_target) && idx_jing_target>=2

                    %存储总步数
                    idx_target_global=find(t1<=T_target,1,'last');
       
                    %提取局部数据
                    t_local=t_jing(1:idx_jing_target);
                    omega_local=omega_jing_unwrap(1:idx_jing_target);

                    % 剔除前1个点
                    valid_idx_local=true(size(t_local));
                    if length(t_local)>=3
                        valid_idx_local(1)=false;
                    end
                    t_local_valid=t_local(valid_idx_local);
                    omega_local_valid=omega_local(valid_idx_local);

                    %判断局部有效数据量是否满足拟合要求
                    if length(t_local_valid)>=2
                        %对局部数据单独线性拟合，得到局部进动率
                        p_local=polyfit(t_local_valid,omega_local_valid, 1);
                        omega_dot_local=p_local(1);%该节点内的实时进动率

                        %计算该节点的累积进动角
                        L_jingdong_node=abs(omega_dot_local*T_nodes(n)*T);

                        %计算该节点的理论进动角
                        L_jingdong_lilun_node=rad2deg(6*pi*GM_sgrA/(c2_AU_year*a_AU*(1-e_val^2))*T_nodes(n));

                        %计算该节点的实时误差并存储
                        if L_jingdong_lilun_node>1e-6
                            node_prec_data(n,1)=abs(L_jingdong_node-L_jingdong_lilun_node)/L_jingdong_lilun_node*100;              
                        else
                            node_prec_data(n,1)=NaN;
                        end
                    else
                        node_prec_data(n,1)=NaN;
                    end
                else
                    node_prec_data(n,1)=NaN;
                end
            end
        end
        wuca_jingdong_array(idx)=wuca_jingdong; 
        L_jingdong_lilun_array(idx)=L_jingdong_lilun_1000T_deg;
    end
    end

    %整合收敛数据
    if abs(e_val-0.8842)<1e-4
        conv_steps=node_step_data;%步数
        conv_prec_err=node_prec_data(:,1);%进动角误差
        conv_energy_err=node_prec_data(:,2);%机械能误差
    end

    %输出结果
    fprintf(['偏心率e=%.4f\n' ...
        '远心点=%.6fAU,近心点=%.6fAU,远点相对误差=%.6f%%,近点相对误差=%.6f%%,近与远误差比值=%.6f\n' ...
        '周期=%.6f,周期相对误差=%.6f%%,机械能相对误差=%.6f%%，角动量误差=%.6f%%\n' ...
        '机械能偏移=%.6f%%,角动量偏移=%.6f%%,角动量方向偏移角度=%.6f度\n' ...
        '近心点进动角=%.6f度,理论值=%.6f度，近心点1000周期进动角相对误差=%.6f%%,计算耗时=%.6fs,积分步数=%d\n\n'],...
        e_val,r_yuan,r_jing,wuca_yuan_array(idx),wuca_jing_array(idx),bizhi_jing_yuan,T_tuo,...
        wuca_zhouqi_array(idx),wuca_E_array(idx),wuca_L_array(idx),...
        pianyi_E,pianyi_L,L_jiaodong_deg,...
        L_jingdong,L_jingdong_lilun_1000T_deg,wuca_jingdong_array(idx),time_cost_e(idx),step_num_e(idx));
end

%保存S2星近心点数据到.mat文件（供后续绘图调用）
save('s2_peri_hermite_data.mat','s2_peri_data');

fprintf('\n===== S2星（e=0.8842）收敛图核心数据（Hermite，单次运行） =====\n');
fprintf('周期节点\t累计积分步数\t近心点进动角相对误差(%%)\t机械能相对误差(%%)\n');
for n=1:length(T_nodes)
    fprintf('%dT\t\t%.0f\t\t\t%.6f\t\t\t\t%.6f\n',...
        T_nodes(n),conv_steps(n),conv_prec_err(n),conv_energy_err(n));
end

%图注设置
xlabel('x(AU)','FontSize',12);
ylabel('y(AU)','FontSize',12);
zlabel('z(AU)','FontSize',12);
title('Sgr A*黑洞-类S2恒星系统轨道模拟（1000T）', 'FontSize', 14, 'FontWeight', 'bold');
view(45,30)
rotate3d on
legend('Location','best','FontSize',10);
axis equal;
axis([-2*a_AU,2*a_AU,-2*a_AU,2*a_AU,-2*a_AU,2*a_AU]);

%验证提取的数据维度
fprintf('\n===== 最终验证：Hermite积分器提取的数据维度 =====\n');
fprintf('t维度: %s\n', mat2str(size(s2_hermite_data.t)));
fprintf('h维度: %s\n', mat2str(size(s2_hermite_data.h)));
fprintf('r维度: %s\n', mat2str(size(s2_hermite_data.r)));
fprintf('f维度: %s\n', mat2str(size(s2_hermite_data.f)));
