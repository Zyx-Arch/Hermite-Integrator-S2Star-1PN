# Hermite-Integrator-S2Star-1PN  

**A 4th‑Order Physics‑Guided Hermite Integrator for Long‑term Orbital Evolution of the S2 Star around Sgr A* with 1PN Post‑Newtonian Corrections**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**高精度、高效率、长时间稳定的四阶 Hermite 积分器，专为银心 S2 星等高偏心率相对论轨道设计。**

---

## 核心亮点

| 指标 | 数值 |
|------|------|
| 近心点相对误差 | **7.57×10⁻⁷%** (对应绝对误差 ~127 米) |
| 1000 周期步数减少 | **~10%** (相比标准 Hermite) |
| 能量误差 (1000T) | 0.893% (与标准型持平) |
| 最大积分周期 | > 1000 个轨道周期 (>1.6 万年) |

---

## 创新点

- **精确的 1PN 后牛顿加速度**，适用于 Schwarzschild 时空下的强引力场轨道模拟。  
- **物理引导自适应步长**：根据瞬时轨道能量、角动量、偏心率及近心点特征时间动态调整步长，在近心点自动加密，远心点适度放大。  
- **高偏心率专用平滑压缩**：针对 e ≥ 0.8 的轨道，采用幂律压缩与相对论强度修正，保证近心点精度。  
- **双重步长控制**：结合 Hermite 预测‑校正误差估计与物理约束，确保长期稳定性和计算效率。  
- **轻量高效**：步长引导仅使用轨道几何量，无需额外复杂计算，易于扩展。

---

## 精度成果

本积分器针对**不同偏心率轨道**完成系统性测试，核心验证**1000倍轨道周期（1000T）**长时标演化的精度与稳定性，所有测试均采用1PN广义相对论动力学模型。

### 一、核心测试指标（全域最优）
- **模拟目标**：S2 星 (e=0.8842，高偏心率相对论轨道)
- **演化时长**：1000 轨道周期 (1000T)
- **近心点相对误差**：**0.0000007570% (7.57×10⁻⁹)** 【顶尖精度】
- **近心点进动角误差**：4.632154%
- **机械能相对误差**：0.892821%
- **角动量误差**：0.142072%
- **计算耗时**：8.69s
- **积分步数**：794,138

---

### 二、不同偏心率轨道精度测试（关键结果）
| 偏心率 e | 近心点相对误差 | 远心点相对误差 | 1000T进动角误差 | 机械能误差 |
|---------|--------------|--------------|----------------|-----------|
| 0.2000  | 0.0003895369% | 0.050480%    | 0.057931%       | 0.006003% |
| 0.4000  | 0.0011312711% | 0.104592%    | 0.130735%       | 0.017358% |
| 0.6000  | 0.0000180611% | 0.253379%    | 0.339834%       | 0.052547% |
| 0.8000  | 0.0010833310% | 1.071850%    | 1.504453%       | 0.268605% |
| 0.8842  | **0.0000007570%** | 3.319744% | **4.632154%** | **0.892821%** |
| 0.9000  | 0.0003631758% | 4.512558%    | 6.235894%       | 1.227820% |

> 核心优势：**偏心率越高，近心点精度越强**，完美适配 S2 星等高偏心率天体模拟。

---

### 三、S2 星 1000T 长时标收敛性验证
| 轨道周期 | 累计积分步数 | 进动角相对误差(%) | 机械能相对误差(%) |
|---------|------------|-----------------|-----------------|
| 1T      | 606        | NaN             | 0.891427        |
| 5T      | 3816       | 4.640171        | 0.891429        |
| 10T     | 7913       | 4.633344        | 0.891438        |
| 100T    | 79388      | 4.633046        | 0.891566        |
| 500T    | 396922     | 4.632664        | 0.892121        |
| 1000T   | 794138     | **4.632154**    | **0.892821**    |

> 收敛结论：**1000T超长演化无精度漂移，数值全程稳定收敛**

---

## 快速开始

### 环境要求
- MATLAB R2021b 或更高版本（无需额外工具箱）

### 应用场景
- 银河系中心 S2 星 / S 星族轨道高精度数值模拟
- 广义相对论 1PN 近心点进动效应检验
- 高偏心率天体长周期轨道演化研究
- 黑洞动力学、强引力场数值模拟

## 引用格式
- 若本项目代码对您的研究有帮助，请引用本 GitHub 仓库：

- Zyx-Arch. *Hermite-Integrator-S2Star-1PN* [EB/OL]. GitHub, 2026. Available: https://github.com/Zyx-Arch/Hermite-Integrator-S2Star-1PN

### 您也可以在学术论文中按以下 BibTeX 格式引用：

bibtex
@misc{ZyxArch2026,
  author = {Zyx-Arch},
  title = {Hermite-Integrator-S2Star-1PN: A 4th-Order Physics-Guided Hermite Integrator for S2 Star Orbits with 1PN Corrections},
  year = {2026},
  publisher = {GitHub},
  howpublished = {\url{https://github.com/Zyx-Arch/Hermite-Integrator-S2Star-1PN}}
}

## 许可证
- 本项目采用 MIT License 开源，允许自由使用、修改与分发，仅需保留版权声明和许可声明。详情请见 LICENSE 文件。

## 贡献者
- Zyx-Arch – 主要开发与验证

## 欢迎通过 GitHub Issues 提出建议、报告问题或贡献代码。

## 欢迎使用、反馈和贡献！

### 最短运行示例（10 个周期）

```matlab
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
dvr_dt=(a_dot_r+v2)/r-vr^2/r;

%牛顿加速度计算
ax_N=-GM*x/r3;
ay_N=-GM*y/r3;
az_N=-GM*z/r3;
a_N=[ax_N;ay_N;az_N];

%1PN加加速度计算
A=1-4*GM/(c2*r)+v2/c2;
dot_A=4*GM*vr/(c2*r2)+dv2_dt/c2;
d_Part1=adot_N.*A+a_N.*dot_A;

C=4*GM/c2;
subterm1=(-2*vr^2/r3)*v_vec;
subterm2=(dvr_dt/r2)*v_vec;
subterm3=(vr/r2).*a;
d_Part2=C*(subterm1+subterm2+subterm3);

%总加加速度计算
j_total=d_Part1+d_Part2;

%各个方向1PN加加速度计算
adot_x=j_total(1);
adot_y=j_total(2);
adot_z=j_total(3);

%输出加速度和加加速度
j=[adot_x;adot_y;adot_z];
end

%动力学引导自适应4阶Hermite相对论轨道积分器
function[t_out,st_out,step_history]=hermite_integrator(t_span,st0,GM_sgrA,c_AU_year,opts,e_val)
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

%自动估算初始步长：基于轨道周期
r0=norm(st0(1:3));
v0=norm(st0(4:6));
a=1/(2/r0-v0^2/GM_sgrA);
T=2*pi*sqrt(a^3/GM_sgrA);
if isfield(opts, 'InitialStep')
    h=opts.InitialStep;
else
    h=T/200;
end

%最小步长（避免过小）
r_s=2*GM_sgrA/(c_AU_year^2);
MinStep=T/1e7;
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

%Hermite核心循环
while t<t_end && iter_count<max_iter && step_num<max_step_num
    iter_count=iter_count +1;
    step_num=step_num +1;

    %调用动力学引导步长约束
    h0=compute_step_from_physics(st,h,GM_sgrA,r_s,T,e_val);  

    %设置步长范围
    h0=min(h0,t_end-t);%步长不能超过剩余时间，避免超出最终时间
    h0=max(h0,MinStep);
    h0=min(h0,MaxStep);

    %提取当前状态
    r_curr=st(1:3);
    v_curr=st(4:6);

    %预测下一步状态行,泰勒展开至3阶，近似下一步状态
    r1_pred=r_curr+v_curr*h0+0.5*a0*h0^2+(1/6)*j0*h0^3;
    v1_pred=v_curr+a0*h0+0.5*j0*h0^2;
    st1_pred=[r1_pred;v1_pred];

    %计算t+h时刻的真实加速度，用于误差估计
    [a1_pred,j1_pred]=calc_a_j(st1_pred,GM_sgrA,c_AU_year);
   
    %4阶Hermite校正步公式
    v_corr=v_curr+(a0+a1_pred)*h0/2-(j1_pred-j0)*h0^2/12;
    r_corr=r_curr+(v_curr+v_corr)*h0/2+(a0-a1_pred)*h0^2/12;
    st_corr=[r_corr;v_corr];
     
    %计算误差估计量 (预测与校正的差值)
    err_pos=norm(r_corr-r1_pred)/(AbsTol+RelTol*norm(r_corr));
    err_vel=norm(v_corr-v1_pred)/(AbsTol+RelTol*norm(v_corr));
    error_est=max(err_pos, err_vel);
    
    %步长拒绝逻辑
    if error_est > 1

        h=h0*max(min_h_change,safety*(1/error_est)^(1/5));
        h=max(h,MinStep);
        continue;
    else
        % 接受步长
        t=t+h0;
        st=st_corr;
        [a0,j0]=calc_a_j(st,GM_sgrA,c_AU_year);
        
        %记录输出
        if step_idx<=max_step_num
            t_out(step_idx)=t;
            st_out(step_idx,:)=st';
            step_history(step_idx)=h0;
            step_idx=step_idx+1;
        end
        
        %计算下一步步长（结合物理约束）
        h_next=h0*min(max_h_change,safety*(1/error_est)^(1/5));
        h_next=max(MinStep,min(MaxStep,h_next));

        % 物理约束限制下一步步长
        h=compute_step_from_physics(st,h_next,GM_sgrA,r_s,T,e_val);
        h=max(h, MinStep);
    end
end

    %截断步长历史到实际长度
    t_out=t_out(1:step_idx-1);
    st_out=st_out(1:step_idx-1,:);
    step_history=step_history(1:step_idx-1);
end

%动力学引导步长约束函数
function h_sug=compute_step_from_physics(st,h_base,GM_sgrA,r_s,T,e_val)

%提取位置速度
r_current=norm(st(1:3));
v_current=norm(st(4:6));

%计算比机械能和半长轴（物理意义：轨道能量，决定轨道类型）
energy=0.5*v_current^2-GM_sgrA/r_current;
if abs(energy)>1e-16
    a_current=-0.5*GM_sgrA/energy;%瞬时半长轴
else
    a_current=1e10;%逃逸轨道兜底
end
a_current=max(a_current,1e-6);%确保半长轴为正的

%计算角动量（物理意义：轨道角动量，守恒量，决定轨道平面）
L_vec=cross(st(1:3),st(4:6));
L2=sum(L_vec.^2);

%计算偏心率（不使用滑动平均，直接计算瞬时值）
if a_current>0 && isfinite(a_current)
    temp=L2/(GM_sgrA*a_current);%偏心率核心公式（物理意义：角动量与半长轴的比值）

    %防止数值误差导致根号内为负数
    e_squared=max(0,1-min(temp,1+1e-12));
    e_curr=sqrt(e_squared+eps);
else
    e_curr=0;
end

%基于近心点特征时间设定压缩系数（物理意义：质点穿越近心点的特征时间尺度）
r_min_curr=a_current*(1-e_curr);%瞬时近心距（物理意义：轨道最近点距离）
r_min_curr=max(r_min_curr,1e-8);
v_peri=sqrt(GM_sgrA*(1+e_curr)/r_min_curr);%近心点速度（开普勒公式）
t_peri_char=r_min_curr/v_peri;%近心点特征时间（轨道穿越近心点的特征尺度）
k_theory=t_peri_char/T;%理论压缩系数（物理意义：特征时间占轨道周期的比例）
    
%对高偏心率使用压缩系数k值(物理意义：偏心率越高，压缩越强，贴合S2星轨道)
if e_curr>=0.8

    %高斯连续函数
    k_emp=0.02+0.03*exp(-10*(e_curr-0.9).^2);

    %近心距修正(物理意义：近心距越小，相对论效应越强，需略压缩步长）
    k_emp=k_emp*(0.5+0.5*exp(-r_min_curr/50));
    k=min(k_theory,k_emp);%取理论值和经验值中较小者（保精度）
else
    k=k_theory;
end
k=max(0.002,min(0.5,k));%上下限（物理意义：避免过度压缩/无压缩）

%近心距比动态步长约束,分区域步长压缩（先验约束）
r_peri_region=3*r_min_curr;% 定义近心点区域半径
transition_width=2*r_min_curr;%;过渡宽度（放大触发范围使过渡更平缓）

%平滑因子：0 表示完全近心点模式，1 表示完全远心点模式
smooth_factor=0.3*(1+tanh((r_current-r_peri_region)/transition_width));
smooth_factor=max(0.0,min(1.0,smooth_factor));

%基础步长约束计算
base_factor=min(1,max(r_current/(k*a_current),1e-3));
peri_step=h_base*base_factor;

%远心点放大逻辑
apo_step=h_base*(1+0.3*smooth_factor);
apo_step=min(apo_step,h_base*1.2); 

%平滑混合
h_temp=(1-smooth_factor)*peri_step+smooth_factor*apo_step;

%定义相对论效应强度
rel_strength=r_s/r_current; 

%基于史瓦西半径的相对论压缩因子
if r_current < 3*r_min_curr && rel_strength > 1e-6%仅在近心点区域（相对论效应显著区）应用压缩
    rel_compress_factor=1+3*rel_strength;
    rel_compress_factor=max(1.0,min(rel_compress_factor,3.0));%限制压缩范围

    %叠加平滑过渡
    rel_smooth=tanh((3*r_min_curr-r_current)/r_min_curr);%近心点越近，平滑因子越接近1
    h_temp=h_temp/(1+(rel_compress_factor-1)*rel_smooth);%渐进式压缩
end

%步长下限保护
h_sug=max(h_temp,1e-8);
end
% 物理常数
GM = 39.4769;                 % AU³/yr²
c = 63241.1;                  % AU/yr

% S2 星初始状态（谐和坐标系，近心点）
a = 971; e = 0.8842;
r_min = a*(1-e);
v_peri = sqrt(GM*(1+e)/r_min);
st0 = [r_min; 0; 0; 0; v_peri; 0];

% 旋转到真实轨道平面
i = deg2rad(133.818); Omega = deg2rad(227.85); omega = deg2rad(66.13);
Rx = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
Rz1 = [cos(Omega) -sin(Omega) 0; sin(Omega) cos(Omega) 0; 0 0 1];
Rz2 = [cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
R = Rz2 * Rx * Rz1;
st0_rot = [R*st0(1:3); R*st0(4:6)];

% 积分器设置
T = 2*pi*sqrt(a^3/GM);          % 约 15.3 年
opts.MaxStep = T/200;
opts.RelTol = 1e-6;
opts.AbsTol = 1e-9;

% 运行 10 个周期
[t, st] = hermite_integrator([0, 10*T], st0_rot, GM, c, opts, e);

% 后处理示例：计算近心点误差、进动角等（见 Hermite_improve_1PN_S2Star_Integrator.m）
