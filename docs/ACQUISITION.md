# 仪器配置与历史脚本

## 安装与配置

先按README建立Python 3.13虚拟环境，再安装 `requirements-acquisition-lock.txt`。PyVISA是Python接口，USB设备还需要适配的厂商VISA runtime和操作系统USB驱动。归档环境使用系统默认VISA；纯Python后端没有在本项目实测。

将 `config.example.json` 复制为 `config.local.json`。填写本机PyVISA识别到的示波器与发生器资源字符串；配置文件被Git忽略。可在本机交互环境中执行以下命令列出资源，输出可能包含序列号，应留在本地：

```python
import pyvisa
rm = pyvisa.ResourceManager()
print(rm.list_resources())
rm.close()
```

配置项：

| 键 | 配置方式 |
|---|---|
| `scope_address` | 填入本机SDS824X HD的完整VISA资源字符串；仓库模板为空 |
| `generator_address` | 填入本机SDG1062X Plus的完整VISA资源字符串；仓库模板为空 |
| `workspace` | 必须是 `local_runs/` 下的子目录；建议每组新实验选择一个新目录 |
| `visa_library` | 空字符串表示系统默认；需要时在本地指定已安装的VISA backend |

脚本没有自动选择“第一台USB设备”，需要两台仪器明确配置。检查配置：

```powershell
.\.venv\Scripts\python.exe acquisition/run.py --check
```

该命令只检查配置格式，不打开VISA，也不能证明设备已连接。执行 `test_instruments.py` 才会实际访问仪器。

## 连接

发生器使用CH2；示波器CH1测量电路实际输入，CH2测量标称5.1 Ω外接电阻的电压。输入均设置1 MΩ、DC、1×，发生器按Hi-Z负载设置。不要把发生器设定幅度当成CH1的Vin。采集流程只适用于所记录的电路和仪器协议；运行前检查实际接线与脚本参数。

示例：

```powershell
.\.venv\Scripts\python.exe acquisition/run.py test_instruments.py
.\.venv\Scripts\python.exe acquisition/run.py characterize_rlc.py
```

启动器将工作目录切换到私有 `local_runs/` 子目录。历史代码中的相对输出也写入该处，公开 `data/` 不会作为采集结果输出。部分步骤需要已有校准或共振文件，启动器仅复制明确的先决参考输入，不自动填入可续采的最终结果表。

## 脚本用途

| 脚本 | 用途与注意事项 |
|---|---|
| `test_instruments.py`、`inspect_instrument_state.py` | 查看身份与状态；身份输出可能包含序列号 |
| `read_scope.py`、`probe_capture.py`、`test_fft.py` | 波形传输/FFT诊断；`test_*` 是硬件诊断脚本，不是自动化单元测试 |
| `characterize_rlc.py`、`refine_resonance.py` | 共振扫描、精细扫描和复数RLC拟合；refine在私有参考副本上继续处理 |
| `characterize_drive_amplitude.py` | 19个幅度的升/降扫描，有续采行为；新实验应使用新的workspace |
| `summarize_drive_amplitude.py` | 早期幅度数据重新校正与绘图；在私有副本上运行，仍需本地配置上下文；不连接硬件 |
| `calibrate_scope_vertical_ranges.py`、`check_scope_scale_dependence.py`、`pilot_low_drive_scale.py` | 量程校准与检查，会设置发生器与示波器 |
| `check_amplitude_linearity.py`、`compare_drive_levels.py` | 幅度/频率对照诊断 |
| `diagnose_am_1khz.py`、`diagnose_am_1khz_500khz.py`、`diagnose_am_1khz_low_drive.py` | 1 kHz AM对照诊断 |
| `sweep_low_drive_snr.py` | 低驱动15频点、3次捕获的SNR筛选；带续采逻辑 |
| `verify_low_drive_repeatability.py`、`verify_low_drive_repeatability_bab.py` | 不同量程/顺序的重复性诊断 |
| `collect_formal_am_dataset.py` | 15个频率 × 5次正式sweep，已完成数据会拒绝覆盖 |
| `collect_formal_am_full_rerun.py` | 31个频率 × 5次正式sweep；以早期汇总作为历史重复性阈值参考 |
| `collect_drive_dependence_support.py` | 0.2–5 Vpp，各3次新测量 |
| `collect_drive_dependence_extreme.py` | **历史恢复专用**：从4 Vpp重复2继续已存在的中断实验，不是完整的新采集入口 |

极限恢复脚本没有配套公开的私有部分运行工作目录，因此默认会拒绝从空目录启动。已公开的21条完成结果与14张图可直接使用；离线脚本不需要执行恢复程序。

## 行为与验证范围

归档脚本保留原测量流程、硬编码的实验幅度/频率和异常恢复行为，仅将设备地址、输出根目录改为本地配置，并增加通过启动器执行的检查。

部分旧诊断脚本在结束时恢复到 **5 Vpp、约493.9 kHz、AM OFF，输出可能仍开启**。`characterize_rlc.py` 等脚本也可能保留最后输出状态。不要假设所有脚本结束都会关闭输出。12 Vpp只用于当时的极限诊断，恢复脚本中保留的检查不能替代核对实际电路。

本次验证包括Python语法、配置拒绝/接受、私有目录隔离、离线统计及发布内容扫描；未连接仪器重新采集，因此不把配置检查或脚本语法通过宣称为硬件联调通过。
