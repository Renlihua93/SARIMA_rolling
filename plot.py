import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 基础工具函数
# ============================================================

def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def safe_to_numeric(df: pd.DataFrame, cols):
    """
    把指定列转成数值。
    遇到 nan、inf、非法字符时转成 NaN。
    """
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
            df[col] = df[col].replace([np.inf, -np.inf], np.nan)
    return df


def set_chinese_font():
    """
    尝试设置中文字体。
    如果系统没有这些字体，也不影响程序运行。
    """
    plt.rcParams["font.sans-serif"] = [
        "Microsoft YaHei",
        "SimHei",
        "Arial Unicode MS",
        "DejaVu Sans"
    ]
    plt.rcParams["axes.unicode_minus"] = False


def get_segment_ylim(values, margin_ratio=0.15):
    """
    根据当前窗口数据动态设置 y 轴范围。
    避免鉴相值整体很大但变化很小时图像看不清。
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if len(arr) == 0:
        return None

    low = np.nanpercentile(arr, 1)
    high = np.nanpercentile(arr, 99)

    if not np.isfinite(low) or not np.isfinite(high):
        return None

    if abs(high - low) < 1e-12:
        center = (high + low) / 2
        return center - 1, center + 1

    margin = (high - low) * margin_ratio
    return low - margin, high + margin


# ============================================================
# 读取数据
# ============================================================

def load_data(mean_file: Path, summary_file: Path, detail_file: Path, points_per_day: int):
    """
    读取 C 程序输出的三个文件。
    """

    if not mean_file.exists():
        raise FileNotFoundError(f"找不到小时平均数据文件: {mean_file}")

    if not summary_file.exists():
        raise FileNotFoundError(f"找不到汇总文件: {summary_file}")

    if not detail_file.exists():
        raise FileNotFoundError(f"找不到详细结果文件: {detail_file}")

    # 小时平均后的原始曲线
    mean_data = pd.read_csv(
        mean_file,
        header=None,
        names=["actual_ns"]
    )

    mean_data["actual_ns"] = pd.to_numeric(mean_data["actual_ns"], errors="coerce")
    mean_data["actual_ns"] = mean_data["actual_ns"].replace([np.inf, -np.inf], np.nan)
    mean_data["day"] = np.arange(len(mean_data)) / points_per_day

    # 每个窗口汇总
    summary = pd.read_csv(summary_file)

    numeric_summary_cols = [
        "case_id",
        "start_day",
        "train_start_day",
        "train_end_day",
        "forecast_days",
        "n_test",
        "linear_mae",
        "linear_rmse",
        "linear_max_abs",
        "linear_exceed_100",
        "sarima_mae",
        "sarima_rmse",
        "sarima_max_abs",
        "sarima_exceed_100",
        "linear_sarima_mae",
        "linear_sarima_rmse",
        "linear_sarima_max_abs",
        "linear_sarima_exceed_100",
    ]

    summary = safe_to_numeric(summary, numeric_summary_cols)

    # 每个预测点详细结果
    detail = pd.read_csv(detail_file)

    numeric_detail_cols = [
        "case_id",
        "start_day",
        "forecast_point",
        "forecast_day",
        "actual",
        "linear_pred",
        "linear_error",
        "sarima_pred",
        "sarima_error",
        "linear_sarima_pred",
        "linear_sarima_error",
    ]

    detail = safe_to_numeric(detail, numeric_detail_cols)

    return mean_data, summary, detail


# ============================================================
# 绘制单个用例：原始曲线 + 预测曲线
# ============================================================

def plot_case_forecast(
    case_id,
    case_summary,
    case_detail,
    mean_data,
    outdir,
    points_per_day=24,
    show_train=True,
    show_full_context=False
):
    """
    绘制某个滚动窗口的:
        原始曲线
        Linear预测
        SARIMA预测
        Linear_SARIMA预测
    """

    train_start_day = float(case_summary["train_start_day"])
    train_end_day = float(case_summary["train_end_day"])

    if len(case_detail) == 0:
        print(f"[跳过] case {case_id}: detail为空")
        return

    forecast_start_day = train_end_day
    forecast_end_day = float(case_detail["forecast_day"].max())

    # 只画当前窗口附近的数据
    if show_full_context:
        raw_seg = mean_data.copy()
    else:
        raw_seg = mean_data[
            (mean_data["day"] >= train_start_day) &
            (mean_data["day"] <= forecast_end_day)
        ].copy()

    train_seg = raw_seg[
        (raw_seg["day"] >= train_start_day) &
        (raw_seg["day"] < train_end_day)
    ]

    test_seg = raw_seg[
        (raw_seg["day"] >= train_end_day) &
        (raw_seg["day"] <= forecast_end_day)
    ]

    fig, ax = plt.subplots(figsize=(14, 6))

    # 训练区间原始曲线
    if show_train and len(train_seg) > 0:
        ax.plot(
            train_seg["day"],
            train_seg["actual_ns"],
            linewidth=1.4,
            label="训练区间原始曲线"
        )

    # 测试区间真实曲线
    if len(test_seg) > 0:
        ax.plot(
            test_seg["day"],
            test_seg["actual_ns"],
            linewidth=1.8,
            label="预测区间真实曲线"
        )

    # 如果 detail 中 actual 更精确，也叠加一下预测区间真实值
    ax.plot(
        case_detail["forecast_day"],
        case_detail["actual"],
        linewidth=1.2,
        linestyle="--",
        label="预测区间真实值(detail)"
    )

    # 三种算法预测曲线
    ax.plot(
        case_detail["forecast_day"],
        case_detail["linear_pred"],
        linewidth=1.3,
        label="Linear预测"
    )

    ax.plot(
        case_detail["forecast_day"],
        case_detail["sarima_pred"],
        linewidth=1.3,
        label="SARIMA预测"
    )

    ax.plot(
        case_detail["forecast_day"],
        case_detail["linear_sarima_pred"],
        linewidth=1.3,
        label="Linear_SARIMA预测"
    )

    # 标记训练结束位置
    ax.axvline(
        train_end_day,
        linestyle=":",
        linewidth=1.5,
        label="训练结束"
    )

    # 训练区间阴影
    ax.axvspan(
        train_start_day,
        train_end_day,
        alpha=0.08,
        label="训练窗口"
    )

    # y轴动态范围
    ylim_values = []

    for col in [
        "actual_ns",
    ]:
        if col in raw_seg.columns:
            ylim_values.extend(raw_seg[col].dropna().values.tolist())

    for col in [
        "actual",
        "linear_pred",
        "sarima_pred",
        "linear_sarima_pred",
    ]:
        if col in case_detail.columns:
            ylim_values.extend(case_detail[col].dropna().values.tolist())

    ylim = get_segment_ylim(ylim_values)
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_title(
        f"Case {case_id:03d}  原始曲线与预测曲线对比\n"
        f"训练: {train_start_day:.2f}d - {train_end_day:.2f}d, "
        f"预测到: {forecast_end_day:.2f}d"
    )
    ax.set_xlabel("时间 / 天")
    ax.set_ylabel("相对鉴相值 / ns")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    fig.tight_layout()

    outpath = outdir / f"case_{case_id:03d}_start_{train_start_day:.2f}d_forecast.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


# ============================================================
# 绘制单个用例：误差曲线
# ============================================================

def plot_case_error(
    case_id,
    case_summary,
    case_detail,
    outdir
):
    """
    绘制某个滚动窗口的误差曲线:
        error = actual - prediction
    """

    train_start_day = float(case_summary["train_start_day"])
    train_end_day = float(case_summary["train_end_day"])

    if len(case_detail) == 0:
        print(f"[跳过] case {case_id}: detail为空")
        return

    forecast_end_day = float(case_detail["forecast_day"].max())

    fig, ax = plt.subplots(figsize=(14, 5))

    ax.plot(
        case_detail["forecast_day"],
        case_detail["linear_error"],
        linewidth=1.3,
        label="Linear误差"
    )

    ax.plot(
        case_detail["forecast_day"],
        case_detail["sarima_error"],
        linewidth=1.3,
        label="SARIMA误差"
    )

    ax.plot(
        case_detail["forecast_day"],
        case_detail["linear_sarima_error"],
        linewidth=1.3,
        label="Linear_SARIMA误差"
    )

    # ±100 ns 阈值线
    ax.axhline(100, linestyle="--", linewidth=1.2, label="+100 ns")
    ax.axhline(-100, linestyle="--", linewidth=1.2, label="-100 ns")
    ax.axhline(0, linestyle=":", linewidth=1.0)

    # y轴动态范围
    ylim_values = []
    for col in [
        "linear_error",
        "sarima_error",
        "linear_sarima_error",
    ]:
        if col in case_detail.columns:
            ylim_values.extend(case_detail[col].dropna().values.tolist())

    ylim = get_segment_ylim(ylim_values)
    if ylim is not None:
        low, high = ylim
        low = min(low, -120)
        high = max(high, 120)
        ax.set_ylim(low, high)

    ax.set_title(
        f"Case {case_id:03d}  预测误差曲线\n"
        f"训练: {train_start_day:.2f}d - {train_end_day:.2f}d, "
        f"预测到: {forecast_end_day:.2f}d"
    )
    ax.set_xlabel("时间 / 天")
    ax.set_ylabel("误差 / ns")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    fig.tight_layout()

    outpath = outdir / f"case_{case_id:03d}_start_{train_start_day:.2f}d_error.png"
    fig.savefig(outpath, dpi=160)
    plt.close(fig)


# ============================================================
# 绘制总览图：每个用例的最大误差
# ============================================================

def plot_summary_max_abs(summary, outdir):
    """
    绘制所有用例的最大绝对误差对比。
    """

    fig, ax = plt.subplots(figsize=(14, 5))

    ax.plot(
        summary["start_day"],
        summary["linear_max_abs"],
        marker="o",
        linewidth=1.3,
        label="Linear最大误差"
    )

    ax.plot(
        summary["start_day"],
        summary["sarima_max_abs"],
        marker="o",
        linewidth=1.3,
        label="SARIMA最大误差"
    )

    ax.plot(
        summary["start_day"],
        summary["linear_sarima_max_abs"],
        marker="o",
        linewidth=1.3,
        label="Linear_SARIMA最大误差"
    )

    ax.axhline(100, linestyle="--", linewidth=1.2, label="100 ns阈值")

    ax.set_title("所有滚动窗口最大绝对误差对比")
    ax.set_xlabel("训练起点 / 天")
    ax.set_ylabel("最大绝对误差 / ns")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    fig.tight_layout()

    outpath = outdir / "summary_max_abs_error_by_case.png"
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def plot_summary_mae(summary, outdir):
    """
    绘制所有用例的 MAE 对比。
    """

    fig, ax = plt.subplots(figsize=(14, 5))

    ax.plot(
        summary["start_day"],
        summary["linear_mae"],
        marker="o",
        linewidth=1.3,
        label="Linear MAE"
    )

    ax.plot(
        summary["start_day"],
        summary["sarima_mae"],
        marker="o",
        linewidth=1.3,
        label="SARIMA MAE"
    )

    ax.plot(
        summary["start_day"],
        summary["linear_sarima_mae"],
        marker="o",
        linewidth=1.3,
        label="Linear_SARIMA MAE"
    )

    ax.set_title("所有滚动窗口 MAE 对比")
    ax.set_xlabel("训练起点 / 天")
    ax.set_ylabel("MAE / ns")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")

    fig.tight_layout()

    outpath = outdir / "summary_mae_by_case.png"
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


def plot_average_max_abs_bar(summary, outdir):
    """
    绘制所有用例最大误差平均值柱状图。
    """

    models = [
        "Linear",
        "SARIMA",
        "Linear_SARIMA",
    ]

    avg_values = [
        summary["linear_max_abs"].dropna().mean(),
        summary["sarima_max_abs"].dropna().mean(),
        summary["linear_sarima_max_abs"].dropna().mean(),
    ]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.bar(models, avg_values)

    for i, v in enumerate(avg_values):
        if np.isfinite(v):
            ax.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=10)

    ax.set_title("所有用例最大绝对误差的平均值")
    ax.set_xlabel("模型")
    ax.set_ylabel("平均最大绝对误差 / ns")
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()

    outpath = outdir / "summary_avg_of_max_abs_error.png"
    fig.savefig(outpath, dpi=180)
    plt.close(fig)


# ============================================================
# 主函数
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="绘制滚动窗口预测结果图：原始曲线 vs 各算法预测曲线，以及误差曲线。"
    )

    parser.add_argument(
        "--mean",
        type=str,
        default="data/IPCLK_mean_data.txt",
        help="小时平均后的原始曲线文件，默认 data/IPCLK_mean_data.txt"
    )

    parser.add_argument(
        "--summary",
        type=str,
        default="data/rolling_summary.csv",
        help="滚动窗口汇总文件，默认 data/rolling_summary.csv"
    )

    parser.add_argument(
        "--detail",
        type=str,
        default="data/rolling_detail.csv",
        help="滚动窗口详细预测文件，默认 data/rolling_detail.csv"
    )

    parser.add_argument(
        "--outdir",
        type=str,
        default="data/rolling_plots",
        help="图片输出目录，默认 data/rolling_plots"
    )

    parser.add_argument(
        "--points-per-day",
        type=int,
        default=24,
        help="每天采样点数。小时级数据为24，默认24。"
    )

    parser.add_argument(
        "--case-start",
        type=int,
        default=None,
        help="只绘制 case_id >= case_start 的用例。默认全部。"
    )

    parser.add_argument(
        "--case-end",
        type=int,
        default=None,
        help="只绘制 case_id <= case_end 的用例。默认全部。"
    )

    parser.add_argument(
        "--full-context",
        action="store_true",
        help="如果设置，则每张图都显示完整原始曲线；默认只显示当前窗口附近曲线。"
    )

    args = parser.parse_args()

    set_chinese_font()

    mean_file = Path(args.mean)
    summary_file = Path(args.summary)
    detail_file = Path(args.detail)
    outdir = Path(args.outdir)

    forecast_dir = outdir / "forecast_compare"
    error_dir = outdir / "error_curve"
    summary_dir = outdir / "summary"

    ensure_dir(forecast_dir)
    ensure_dir(error_dir)
    ensure_dir(summary_dir)

    mean_data, summary, detail = load_data(
        mean_file=mean_file,
        summary_file=summary_file,
        detail_file=detail_file,
        points_per_day=args.points_per_day
    )

    # 去掉 case_id 为空的行
    summary = summary.dropna(subset=["case_id"]).copy()
    detail = detail.dropna(subset=["case_id"]).copy()

    summary["case_id"] = summary["case_id"].astype(int)
    detail["case_id"] = detail["case_id"].astype(int)

    # 筛选用例范围
    if args.case_start is not None:
        summary = summary[summary["case_id"] >= args.case_start]

    if args.case_end is not None:
        summary = summary[summary["case_id"] <= args.case_end]

    case_ids = summary["case_id"].drop_duplicates().sort_values().tolist()

    print(f"读取 mean_data 点数: {len(mean_data)}")
    print(f"读取 summary 用例数: {len(summary)}")
    print(f"读取 detail 行数: {len(detail)}")
    print(f"准备绘制用例数: {len(case_ids)}")

    # 按 case_id 分组，避免每次循环重复筛选全表
    detail_grouped = dict(tuple(detail.groupby("case_id")))

    for idx, case_id in enumerate(case_ids, start=1):
        case_summary_df = summary[summary["case_id"] == case_id]

        if len(case_summary_df) == 0:
            continue

        case_summary = case_summary_df.iloc[0]
        case_detail = detail_grouped.get(case_id, pd.DataFrame()).copy()

        if len(case_detail) == 0:
            print(f"[跳过] case {case_id}: 没有 detail 数据")
            continue

        case_detail = case_detail.sort_values("forecast_day").copy()

        plot_case_forecast(
            case_id=case_id,
            case_summary=case_summary,
            case_detail=case_detail,
            mean_data=mean_data,
            outdir=forecast_dir,
            points_per_day=args.points_per_day,
            show_train=True,
            show_full_context=args.full_context
        )

        plot_case_error(
            case_id=case_id,
            case_summary=case_summary,
            case_detail=case_detail,
            outdir=error_dir
        )

        print(f"[{idx}/{len(case_ids)}] 已绘制 case {case_id}")

    # 总览图
    plot_summary_max_abs(summary, summary_dir)
    plot_summary_mae(summary, summary_dir)
    plot_average_max_abs_bar(summary, summary_dir)

    print("\n绘图完成。")
    print(f"预测对比图目录: {forecast_dir}")
    print(f"误差曲线图目录: {error_dir}")
    print(f"汇总图目录: {summary_dir}")


if __name__ == "__main__":
    main()