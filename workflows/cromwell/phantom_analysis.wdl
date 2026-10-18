version 1.0

# Phantom Science FASTQ report workflow
#
# A WDL demonstration that reads FASTQ files written by `phantom-science
# generate` and writes text, HTML and JSON reports: line-count and N-content
# checks, quality-character counts, and file sizes and read counts.
#
# It does not detect or verify fingerprints (use `phantom-science verify`),
# and it does not measure performance: the "estimated processing time" assumes
# a fixed 50 MB/s. customer_id, analysis_type, min_quality_score,
# fingerprint_strategy, enable_steganography and batch_size are only echoed
# into the reports; memory_limit_gb and cpu_cores set the runtime of the
# PerformanceAnalysis task.
#
# Checked with `womtool validate`; see README.md for what has been run.

workflow PhantomAnalysis {
    input {
        # Input files and parameters
        Array[File] fastq_files
        String output_dir = "./cromwell_output"
        String customer_id = "phantom-science"
        String analysis_type = "comprehensive"
        
        # Validation parameters
        Int min_quality_score = 20
        Int max_n_percentage = 5
        
        # Fingerprinting parameters
        String fingerprint_strategy = "balanced"
        Boolean enable_steganography = true
        
        # Performance parameters
        Int batch_size = 1000
        Float memory_limit_gb = 8.0
        Int cpu_cores = 4
    }
    
    # Step 1: Validate FASTQ files
    call BatchValidation {
        input:
            fastq_files = fastq_files,
            min_quality_score = min_quality_score,
            max_n_percentage = max_n_percentage,
            output_dir = output_dir
    }
    
    # Step 2: Steganographic fingerprint analysis
    call FingerprintAnalysis {
        input:
            fastq_files = fastq_files,
            validation_report = BatchValidation.validation_report,
            strategy = fingerprint_strategy,
            enable_steganography = enable_steganography,
            customer_id = customer_id,
            output_dir = output_dir
    }
    
    # Step 3: Performance analysis
    call PerformanceAnalysis {
        input:
            fastq_files = fastq_files,
            fingerprint_report = FingerprintAnalysis.fingerprint_report,
            batch_size = batch_size,
            memory_limit_gb = memory_limit_gb,
            cpu_cores = cpu_cores,
            output_dir = output_dir
    }
    
    # Step 4: Generate comprehensive report
    call GenerateReport {
        input:
            validation_report = BatchValidation.validation_report,
            fingerprint_report = FingerprintAnalysis.fingerprint_report,
            performance_report = PerformanceAnalysis.performance_report,
            customer_id = customer_id,
            analysis_type = analysis_type,
            output_dir = output_dir
    }
    
    output {
        File validation_report = BatchValidation.validation_report
        File fingerprint_report = FingerprintAnalysis.fingerprint_report
        File performance_report = PerformanceAnalysis.performance_report
        File comprehensive_report = GenerateReport.html_report
        File workflow_summary = GenerateReport.summary_json
    }
}

# Task: Validate FASTQ files for quality and format compliance
task BatchValidation {
    input {
        Array[File] fastq_files
        Int min_quality_score
        Int max_n_percentage
        String output_dir
    }
    
    command <<<
        set -eo pipefail
        
        # Create output directory
        mkdir -p ~{output_dir}/validation
        
        # Initialize validation report
        echo "# FASTQ Validation Report" > ~{output_dir}/validation/validation_report.txt
        echo "Generated: $(date)" >> ~{output_dir}/validation/validation_report.txt
        echo "" >> ~{output_dir}/validation/validation_report.txt
        
        # Validate each FASTQ file
        for fastq in ~{sep=' ' fastq_files}; do
            echo "Validating: $fastq" >> ~{output_dir}/validation/validation_report.txt
            
            # Check file format
            if [[ "$fastq" == *.gz ]]; then
                cat_cmd="zcat"
            else
                cat_cmd="cat"
            fi
            
            # Basic FASTQ validation
            lines=$($cat_cmd "$fastq" | wc -l)
            if (( lines % 4 != 0 )); then
                echo "  ERROR: Invalid FASTQ format - line count not divisible by 4" >> ~{output_dir}/validation/validation_report.txt
                continue
            fi
            
            # Quality scores are not checked (min_quality_score is not used)
            echo "  Quality scores: not checked" >> ~{output_dir}/validation/validation_report.txt
            
            # N percentage check
            sequences=$($cat_cmd "$fastq" | awk 'NR%4==2 && NR<=4000')
            n_count=$(echo "$sequences" | tr -cd 'N' | wc -c)
            total_bases=$(echo "$sequences" | tr -d '\n' | wc -c)
            
            if (( total_bases > 0 )); then
                n_percentage=$(( (n_count * 100) / total_bases ))
                if (( n_percentage <= ~{max_n_percentage} )); then
                    echo "  N percentage: ${n_percentage}% (PASS)" >> ~{output_dir}/validation/validation_report.txt
                else
                    echo "  N percentage: ${n_percentage}% (FAIL)" >> ~{output_dir}/validation/validation_report.txt
                fi
            fi
            
            echo "  Status: checked (line count; N content of the first 1000 reads)" >> ~{output_dir}/validation/validation_report.txt
            echo "" >> ~{output_dir}/validation/validation_report.txt
        done
        
        echo "Validation completed: $(date)" >> ~{output_dir}/validation/validation_report.txt
    >>>
    
    output {
        File validation_report = "${output_dir}/validation/validation_report.txt"
    }
    
    runtime {
        docker: "ubuntu:22.04"
        memory: "2 GB"
        cpu: 2
        disks: "local-disk 20 HDD"
        continueOnReturnCode: [0, 141]
    }
}

# Task: Count quality characters. This task does not read fingerprints; it
# reports the quality-string statistics and points to `phantom-science verify`.
task FingerprintAnalysis {
    input {
        Array[File] fastq_files
        File validation_report
        String strategy
        Boolean enable_steganography
        String customer_id
        String output_dir
    }
    
    command <<<
        set -eo pipefail
        
        # Create output directory
        mkdir -p ~{output_dir}/fingerprinting
        
        # Initialize fingerprint report
        echo "# Quality Character Report" > ~{output_dir}/fingerprinting/fingerprint_report.txt
        echo "Generated: $(date)" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        echo "Strategy: ~{strategy}" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        echo "Customer ID: ~{customer_id}" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        echo "Steganography Enabled: ~{enable_steganography}" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        echo "" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        
        # Count quality characters in each FASTQ file
        for fastq in ~{sep=' ' fastq_files}; do
            echo "Counting quality characters in: $fastq" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
            
            # Check file format
            if [[ "$fastq" == *.gz ]]; then
                cat_cmd="zcat"
            else
                cat_cmd="cat"
            fi
            
            # Quality lines of the first 1000 reads
            quality_lines=$($cat_cmd "$fastq" | awk 'NR%4==0 && NR<=4000')
            
            # Calculate quality score statistics
            total_qualities=$(echo "$quality_lines" | tr -d '\n' | wc -c)
            unique_qualities=$(echo "$quality_lines" | tr -d '\n' | fold -w1 | sort | uniq | wc -l)
            
            echo "  Total quality characters: $total_qualities" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
            echo "  Unique quality characters: $unique_qualities" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
            
            # Fingerprints are encrypted and can only be extracted with the
            # stegano library and the key; this workflow does neither.
            echo "  Fingerprint: not checked (run phantom-science verify --key-file KEY FILE)" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
            
            echo "  Status: counted" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
            echo "" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
        done
        
        echo "Quality character report completed: $(date)" >> ~{output_dir}/fingerprinting/fingerprint_report.txt
    >>>
    
    output {
        File fingerprint_report = "${output_dir}/fingerprinting/fingerprint_report.txt"
    }
    
    runtime {
        docker: "ubuntu:22.04"
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 HDD"
        continueOnReturnCode: [0, 141]
    }
}

# Task: Report file sizes and read counts. The processing time and throughput
# are computed from an assumed 50 MB/s, not measured.
task PerformanceAnalysis {
    input {
        Array[File] fastq_files
        File fingerprint_report
        Int batch_size
        Float memory_limit_gb
        Int cpu_cores
        String output_dir
    }
    
    command <<<
        set -eo pipefail
        
        # Create output directory
        mkdir -p ~{output_dir}/performance
        
        # Initialize performance report
        echo "# Performance Analysis Report" > ~{output_dir}/performance/performance_report.txt
        echo "Generated: $(date)" >> ~{output_dir}/performance/performance_report.txt
        echo "Batch Size: ~{batch_size}" >> ~{output_dir}/performance/performance_report.txt
        echo "Memory Limit: ~{memory_limit_gb} GB" >> ~{output_dir}/performance/performance_report.txt
        echo "CPU Cores: ~{cpu_cores}" >> ~{output_dir}/performance/performance_report.txt
        echo "" >> ~{output_dir}/performance/performance_report.txt
        
        # Analyze performance for each FASTQ file
        total_files=0
        total_size=0
        total_reads=0
        
        for fastq in ~{sep=' ' fastq_files}; do
            echo "Performance analysis for: $fastq" >> ~{output_dir}/performance/performance_report.txt
            
            # File size analysis
            file_size=$(stat -c%s "$fastq" 2>/dev/null || stat -f%z "$fastq")
            file_size_mb=$((file_size / 1024 / 1024))
            
            # Read count analysis
            if [[ "$fastq" == *.gz ]]; then
                cat_cmd="zcat"
            else
                cat_cmd="cat"
            fi
            
            read_count=$($cat_cmd "$fastq" | wc -l)
            read_count=$((read_count / 4))  # 4 lines per read
            
            # Processing time simulation
            processing_time=$(awk "BEGIN {printf \"%.2f\", $file_size_mb / 50}")  # ~50 MB/s processing
            
            echo "  File size: ${file_size_mb} MB" >> ~{output_dir}/performance/performance_report.txt
            echo "  Read count: $read_count" >> ~{output_dir}/performance/performance_report.txt
            echo "  Estimated processing time at an assumed 50 MB/s: ${processing_time} seconds" >> ~{output_dir}/performance/performance_report.txt
            
            # Throughput calculation
            if (( $(awk "BEGIN {print ($processing_time > 0)}") )); then
                throughput=$(awk "BEGIN {printf \"%.2f\", $read_count / $processing_time}")
                echo "  Throughput at the assumed rate: ${throughput} reads/second" >> ~{output_dir}/performance/performance_report.txt
            fi
            
            echo "  Status: ANALYZED" >> ~{output_dir}/performance/performance_report.txt
            echo "" >> ~{output_dir}/performance/performance_report.txt
            
            total_files=$((total_files + 1))
            total_size=$((total_size + file_size_mb))
            total_reads=$((total_reads + read_count))
        done
        
        # Summary statistics
        echo "## Performance Summary" >> ~{output_dir}/performance/performance_report.txt
        echo "Total files processed: $total_files" >> ~{output_dir}/performance/performance_report.txt
        echo "Total data size: ${total_size} MB" >> ~{output_dir}/performance/performance_report.txt
        echo "Total reads: $total_reads" >> ~{output_dir}/performance/performance_report.txt
        
        if (( total_files > 0 )); then
            avg_file_size=$((total_size / total_files))
            avg_reads_per_file=$((total_reads / total_files))
            echo "Average file size: ${avg_file_size} MB" >> ~{output_dir}/performance/performance_report.txt
            echo "Average reads per file: ${avg_reads_per_file}" >> ~{output_dir}/performance/performance_report.txt
        fi
        
        echo "Performance analysis completed: $(date)" >> ~{output_dir}/performance/performance_report.txt
    >>>
    
    output {
        File performance_report = "${output_dir}/performance/performance_report.txt"
    }
    
    runtime {
        docker: "ubuntu:22.04"
        memory: "${memory_limit_gb} GB"
        cpu: cpu_cores
        disks: "local-disk 20 HDD"
        continueOnReturnCode: [0, 141]
    }
}

# Task: Combine the reports into HTML and JSON
task GenerateReport {
    input {
        File validation_report
        File fingerprint_report
        File performance_report
        String customer_id
        String analysis_type
        String output_dir
    }
    
    command <<<
        set -eo pipefail
        
        # Create output directory
        mkdir -p ~{output_dir}/reports
        
        # Generate HTML report
        cat > ~{output_dir}/reports/phantom_analysis_report.html << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Phantom Science Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background-color: #f0f0f0; padding: 20px; border-radius: 5px; }
        .section { margin: 20px 0; padding: 15px; border: 1px solid #ddd; border-radius: 5px; }
        .success { color: green; }
        .error { color: red; }
        .info { color: blue; }
        pre { background-color: #f8f8f8; padding: 10px; border-radius: 3px; overflow-x: auto; }
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 Phantom Science Analysis Report</h1>
        <p><strong>Customer ID:</strong> ~{customer_id}</p>
        <p><strong>Analysis Type:</strong> ~{analysis_type}</p>
        <p><strong>Generated:</strong> $(date)</p>
    </div>
    
    <div class="section">
        <h2>📊 Validation Results</h2>
        <pre>
EOF
        
        # Include validation report content
        cat ~{validation_report} >> ~{output_dir}/reports/phantom_analysis_report.html
        
        cat >> ~{output_dir}/reports/phantom_analysis_report.html << 'EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>Quality characters</h2>
        <pre>
EOF
        
        # Include fingerprint report content
        cat ~{fingerprint_report} >> ~{output_dir}/reports/phantom_analysis_report.html
        
        cat >> ~{output_dir}/reports/phantom_analysis_report.html << 'EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>⚡ Performance Analysis</h2>
        <pre>
EOF
        
        # Include performance report content
        cat ~{performance_report} >> ~{output_dir}/reports/phantom_analysis_report.html
        
        cat >> ~{output_dir}/reports/phantom_analysis_report.html << 'EOF'
        </pre>
    </div>
    
    <div class="section">
        <h2>📋 Summary</h2>
        <p>This report includes:</p>
        <ul>
            <li>FASTQ line-count and N-content checks</li>
            <li>Quality-character counts (fingerprints are not checked; use phantom-science verify)</li>
            <li>File sizes, read counts and an estimate at an assumed 50 MB/s</li>
        </ul>
        <p><em>Generated by Phantom Science Cromwell Workflow</em></p>
    </div>
</body>
</html>
EOF
        
        # Generate JSON summary
        cat > ~{output_dir}/reports/workflow_summary.json << EOF
{
    "workflow_name": "PhantomAnalysis",
    "customer_id": "~{customer_id}",
    "analysis_type": "~{analysis_type}",
    "generated_timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "reports": {
        "validation_report": "~{validation_report}",
        "fingerprint_report": "~{fingerprint_report}",
        "performance_report": "~{performance_report}",
        "html_report": "~{output_dir}/reports/phantom_analysis_report.html"
    },
    "status": "completed"
}
EOF
        
        echo "Report generation completed: $(date)"
    >>>
    
    output {
        File html_report = "${output_dir}/reports/phantom_analysis_report.html"
        File summary_json = "${output_dir}/reports/workflow_summary.json"
    }
    
    runtime {
        docker: "ubuntu:22.04"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 20 HDD"
        continueOnReturnCode: [0, 141]
    }
}