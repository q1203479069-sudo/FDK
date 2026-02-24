# FDK MATLAB 重建实现

本仓库提供一个基于 MATLAB 的 FDK（Feldkamp-Davis-Kress）锥束 CT 重建实现。

## 文件位置

- **主函数**：`fdk_reconstruction.m`
  - 入口函数：`volume = fdk_reconstruction(proj, geo)`
- **说明文档**：`README.md`（当前文件）

## 快速使用

```matlab
% proj: [nu, nv, nBeta]
% geo:  几何参数结构体
volume = fdk_reconstruction(proj, geo);
```

## geo 参数要求

必需字段：

- `DSO`, `DSD`
- `du`, `dv`
- `nu`, `nv`
- `nx`, `ny`, `nz`
- `dx`, `dy`, `dz`
- `beta`（弧度制，长度等于投影数）

可选字段：

- `filterType`：`'ram-lak'`（默认）、`'shepp-logan'`、`'cosine'`

## 说明

实现流程包含：

1. 距离加权
2. 沿探测器 `u` 方向滤波
3. 锥束反投影与角度积分
