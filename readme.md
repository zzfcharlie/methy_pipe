##  项目结构

```bash
.
├── config.R                 # 主配置文件，定义参数和路径
├── main.R                  # 主分析脚本，运行所有分析流程
├── Dockerfile              
├── Methy_pipeline/         # 包含所有分析函数和模板文件
├── data/
│   ├── 208933290067_R01C01_Grn.idat # 示例数据
│   ├── 208933290067_R01C01_Red.idat # 示例数据
│   ├── sample_mapping.xlsx # 样本编号与芯片文件映射关系
│   └── sample_info.csv     # 样本采集时间与信息
├── output/                 # 分析结果输出路径
```

##  使用方式



```bash
docker pull zzfcharlie/methy_pipeline:latest


docker run --rm \
  -v $(pwd):/project \
  -w /project \
  zzfcharlie/methy_pipeline:latest \
  --config=config.R \
  --idat-dir=./data \
  --sample-mapping-path=./data/sample_mapping_test.xlsx \
  --sample-info-path=./data/sample_info_test.csv \
  --output-dir=./output
```

### 参数说明

- `--config`: 配置文件路径(config.R)
- `--idat-dir`: 原始idat数据目录
- `--sample-mapping-path`: 样本文件对应关系（sample_mapping_test.xlsx）
- `--sample-info-path`: 样本信息表(sample_info.csv)
- `--output-dir`: 分析结果的输出目录