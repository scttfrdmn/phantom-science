#!/usr/bin/env python3
"""
Performance regression detection for Phantom Science Enterprise
Compares current performance metrics against established baselines
"""

import json
import sys
import argparse
from typing import Dict, Any, List, Tuple
from dataclasses import dataclass
from enum import Enum

class RegressionSeverity(Enum):
    NONE = "none"
    MINOR = "minor"
    MAJOR = "major" 
    CRITICAL = "critical"

@dataclass
class PerformanceMetric:
    name: str
    current_value: float
    baseline_value: float
    unit: str
    higher_is_better: bool = True
    
    @property
    def percentage_change(self) -> float:
        if self.baseline_value == 0:
            return 0.0
        return ((self.current_value - self.baseline_value) / self.baseline_value) * 100
    
    @property
    def regression_severity(self) -> RegressionSeverity:
        change = self.percentage_change
        
        # Invert logic if lower values are better
        if not self.higher_is_better:
            change = -change
            
        if change >= -5:  # Within 5% tolerance
            return RegressionSeverity.NONE
        elif change >= -15:  # 5-15% regression
            return RegressionSeverity.MINOR
        elif change >= -30:  # 15-30% regression
            return RegressionSeverity.MAJOR
        else:  # >30% regression
            return RegressionSeverity.CRITICAL

def load_performance_data(file_path: str) -> Dict[str, Any]:
    """Load performance data from JSON file"""
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Warning: Performance file {file_path} not found")
        return {}
    except json.JSONDecodeError as e:
        print(f"Error: Invalid JSON in {file_path}: {e}")
        sys.exit(1)

def extract_metrics(data: Dict[str, Any]) -> List[PerformanceMetric]:
    """Extract standardized metrics from performance data"""
    metrics = []
    
    # Standard benchmark metrics
    benchmark_mapping = {
        "small": ("Small Dataset Throughput", "MB/s", True),
        "medium": ("Medium Dataset Throughput", "MB/s", True), 
        "large": ("Large Dataset Throughput", "MB/s", True),
        "xlarge": ("XLarge Dataset Throughput", "MB/s", True),
        "gpu-stress": ("GPU Stress Test Throughput", "MB/s", True)
    }
    
    for test_name, test_data in data.items():
        if test_name in benchmark_mapping and isinstance(test_data, dict):
            metric_name, unit, higher_better = benchmark_mapping[test_name]
            
            if "throughput" in test_data:
                metrics.append(PerformanceMetric(
                    name=metric_name,
                    current_value=float(test_data["throughput"]),
                    baseline_value=0.0,  # Will be set when comparing
                    unit=unit,
                    higher_is_better=higher_better
                ))
                
            if "duration" in test_data:
                metrics.append(PerformanceMetric(
                    name=f"{metric_name} Duration",
                    current_value=float(test_data["duration"]),
                    baseline_value=0.0,
                    unit="seconds",
                    higher_is_better=False  # Lower duration is better
                ))
    
    return metrics

def compare_performance(baseline_file: str, current_file: str, threshold: float) -> Tuple[List[PerformanceMetric], bool]:
    """Compare current performance against baseline"""
    baseline_data = load_performance_data(baseline_file)
    current_data = load_performance_data(current_file)
    
    current_metrics = extract_metrics(current_data)
    baseline_metrics = extract_metrics(baseline_data)
    
    # Create baseline lookup
    baseline_lookup = {metric.name: metric.current_value for metric in baseline_metrics}
    
    # Update current metrics with baseline values
    compared_metrics = []
    for metric in current_metrics:
        if metric.name in baseline_lookup:
            metric.baseline_value = baseline_lookup[metric.name]
            compared_metrics.append(metric)
    
    # Check if any critical regressions exist
    has_regression = any(
        metric.regression_severity in [RegressionSeverity.MAJOR, RegressionSeverity.CRITICAL]
        for metric in compared_metrics
    )
    
    return compared_metrics, has_regression

def generate_report(metrics: List[PerformanceMetric], output_file: str = None):
    """Generate performance comparison report"""
    report_lines = [
        "# Performance Regression Report",
        "",
        f"**Analysis Date**: {__import__('datetime').datetime.now().isoformat()}",
        "",
        "## Summary",
        ""
    ]
    
    # Count regressions by severity
    severity_counts = {severity: 0 for severity in RegressionSeverity}
    for metric in metrics:
        severity_counts[metric.regression_severity] += 1
    
    report_lines.extend([
        f"- **Total Metrics**: {len(metrics)}",
        f"- **No Regression**: {severity_counts[RegressionSeverity.NONE]}",
        f"- **Minor Regressions**: {severity_counts[RegressionSeverity.MINOR]}",
        f"- **Major Regressions**: {severity_counts[RegressionSeverity.MAJOR]}",
        f"- **Critical Regressions**: {severity_counts[RegressionSeverity.CRITICAL]}",
        "",
        "## Detailed Results",
        "",
        "| Metric | Current | Baseline | Change | Status |",
        "|--------|---------|----------|---------|--------|"
    ])
    
    for metric in sorted(metrics, key=lambda m: m.regression_severity.value, reverse=True):
        status_emoji = {
            RegressionSeverity.NONE: "✅",
            RegressionSeverity.MINOR: "⚠️",
            RegressionSeverity.MAJOR: "❌", 
            RegressionSeverity.CRITICAL: "🚨"
        }
        
        change_str = f"{metric.percentage_change:+.1f}%"
        status = f"{status_emoji[metric.regression_severity]} {metric.regression_severity.value.title()}"
        
        report_lines.append(
            f"| {metric.name} | {metric.current_value:.2f} {metric.unit} | "
            f"{metric.baseline_value:.2f} {metric.unit} | {change_str} | {status} |"
        )
    
    report_lines.extend([
        "",
        "## Performance Insights",
        ""
    ])
    
    # Add performance insights
    if severity_counts[RegressionSeverity.CRITICAL] > 0:
        report_lines.append("⚠️ **CRITICAL**: Major performance degradation detected. Immediate investigation required.")
    elif severity_counts[RegressionSeverity.MAJOR] > 0:
        report_lines.append("⚠️ **WARNING**: Significant performance regression detected. Review required.")
    elif severity_counts[RegressionSeverity.MINOR] > 0:
        report_lines.append("ℹ️ **INFO**: Minor performance variations detected. Monitoring recommended.")
    else:
        report_lines.append("✅ **SUCCESS**: No performance regressions detected. All metrics within acceptable ranges.")
    
    report_content = "\n".join(report_lines)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_content)
        print(f"Performance report saved to: {output_file}")
    else:
        print(report_content)

def main():
    parser = argparse.ArgumentParser(description="Compare performance metrics against baseline")
    parser.add_argument("--baseline", required=True, help="Baseline performance JSON file")
    parser.add_argument("--current", required=True, help="Current performance JSON file") 
    parser.add_argument("--threshold", type=float, default=5.0, help="Regression threshold percentage (default: 5.0)")
    parser.add_argument("--output", help="Output report file (optional)")
    parser.add_argument("--fail-on-regression", action="store_true", help="Exit with error code if regression detected")
    
    args = parser.parse_args()
    
    # Compare performance
    metrics, has_regression = compare_performance(args.baseline, args.current, args.threshold)
    
    if not metrics:
        print("Warning: No comparable metrics found between baseline and current data")
        return
    
    # Generate report
    generate_report(metrics, args.output)
    
    # Exit with appropriate code
    if args.fail_on_regression and has_regression:
        print("\n❌ Performance regression detected!")
        sys.exit(1)
    else:
        print("\n✅ Performance analysis complete")
        sys.exit(0)

if __name__ == "__main__":
    main()