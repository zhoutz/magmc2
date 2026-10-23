# StepperMagnus：偏振步进与 dense output

`src/magnus.hpp` 求解 `S' = omega(l) × S`，其中 `S = (Q,U,V)`，不演化强度
`I`。`omega(l)` 返回 `double3`，单位为路径长度的倒数，分量与 Stokes 使用同一
坐标系。正的 `omega.z` 将正 Q 向正 U 旋转。

积分器保存回调的常量引用，回调对象必须比积分器存活更久；不能传临时 lambda。
无需预先计算 `dS/dl`。偏振度由旋转保持，不强制将输入归一化；部分偏振和零向量
均可使用。磁场对齐坐标变换、公共相位、绝热模式切换和冻结判据由调用方负责。

## 接口与生命周期

```cpp
auto omega = [&](double l) -> double3 {
  auto geometry = geodesic.interp(l);
  return compute_omega(geometry); // 上层提供物理系数
};
StepperMagnus polarization(omega, atol, rtol);
polarization.init(l_old, h_initial, l_new, stokes);
```

- `init(x_init, h_init, x_end, y_init)` 初始化当前测地线区间并清除旧缓存；只支持
  正向积分，允许空区间。`h_init` 必须为有限正数。
- `do_step(max_step)` 接受一个偏振步并准备缓存，返回 `true`。默认不额外限制步长，
  但始终受 `x_end` 限制。当前区间全部提交后返回 `false`。
- 接受后，`[x_old,x_new]` 为当前偏振步，`h_old=x_new-x_old`，`y_old/y_new`
  为两端状态；`h_new` 是下一步的建议步长，可能超过下个测地线区间长度。
- `dense_out(l)` 只在当前接受步内有效，返回任意位置的 Stokes；端点直接返回缓存状态。
  查询不会再次计算测地线、磁场或 `omega`。
- `update_old()` 提交已接受步，令 `x_old=x_new`、`y_old=y_new`、`h_old=h_new`，
  并使 dense output 失效。必须提交后才能开始下一步。公开状态供读取，不应在步进中手动改写。

`atol` 和 `rtol` 必须有限、非负且至少一个大于零。接受条件是

```text
norm(y_err) <= atol + rtol * max(norm(y_old), norm(y_new))
y_err = (S_two_half_steps - S_one_whole_step) / 15
```

接受解为两半步结果，不作 Richardson 外推。容差约束端点的局部误差估计，**不是**
全局误差上界，也没有独立的步内误差控制。光滑生成元下，端点全局四阶、步内局部
误差为五阶；不光滑系数或任意大的变化相位不保证该渐近表现。

## 测地线 → 偏振 → 光深

下面是接入 `evolve()` 的流程示意；`compute_omega`、`integrate_opacity`、
`locate_scattering` 及结果处理函数由上层提供，不属于积分器。一次散射前的剩余
光深 `tau_remaining` 在外层初始化一次，跨测地线步和偏振步累计扣减。

```cpp
for (;;) {
  bool absorbed = geodesic.step();
  polarization.init(geodesic.l_old(), h_pol, geodesic.l_new(), stokes);

  while (polarization.do_step()) {
    double a = polarization.x_old;
    double b = polarization.x_new;
    auto opacity = [&](double l) {
      auto geometry = geodesic.interp(l);
      auto s = polarization.dense_out(l);
      return compute_opacity(geometry, s);
    };
    // integrate_opacity 在积分区间内采样 opacity。
    // 建议在 std::midpoint(a,b) 分段；物理共振边界也由上层分段。
    double dtau = integrate_opacity(opacity, a, b);
    if (dtau >= tau_remaining) {
      // 求解 integral(a,l,opacity) == tau_remaining。
      double l = locate_scattering(opacity, a, b, tau_remaining);
      auto geometry = geodesic.interp(l);
      auto incident_stokes = polarization.dense_out(l);
      return handle_scattering(l, geometry, incident_stokes);
    }
    tau_remaining -= dtau;
    polarization.update_old();
  }

  stokes = polarization.y_old;
  h_pol = polarization.h_new;
  // geodesic.step() 可能已截断到吸收点；先检查沿途散射，再处理终点吸收。
  if (absorbed)
    return handle_absorption();
  // 上层的逃逸终点也应在此处理。
  geodesic.stepper.update_old();
}
```

整个内层循环及散射定位期间，必须保留当前测地线的 dense cache；不能提前调用
`geodesic.stepper.update_old()` 或下一次 `geodesic.step()`。散射发生时，不应先提交
包含散射点的整个偏振步，而应立即读取该点状态并交给上层，随后重新初始化新轨道。

没有接受步、提交后或重新初始化后查询 dense output 会抛 `std::logic_error`；
越界或 NaN 查询抛 `std::out_of_range`。非法初始化参数或容差抛
`std::invalid_argument`；非有限系数、数值溢出或无法表示两个半步时抛
`std::runtime_error`。空区间的 `do_step()` 返回 `false`，不调用生成元。

## 可运行验证

`tests/magnus.cpp` 包含一个完整的三层嵌套例子：使用实际 `Geodesic` 的平直空间
径向轨道，取 `omega=(0,0,1)`、`S(0)=(1,0,0)` 和示例光深率 `2+Q/2`，其解析
累计光深为 `2*l+sin(l)/2`。测试跨多个测地线步定位 `l=0.37` 的散射事件。
这里的光深率用于接口验证，不代表项目的物理散射模型。

在仓库根目录运行：

```sh
mkdir -p build
g++-16 -std=c++23 -O2 -Wall -Wextra -pedantic -I src tests/magnus.cpp -o build/test_magnus
build/test_magnus
```

测试还覆盖常生成元的大相位解析解、旋转轴变化的解析解与收敛阶、拒步和容差收紧、
偏振度守恒、两半步连续性、无回调查询、跨区间状态传递及异常处理。
