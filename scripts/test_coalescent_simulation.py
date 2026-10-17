#!/usr/bin/env python3
"""
Coalescent Simulation Testing and Memory Profiling
Phase 1, Week 1 - Testing memory requirements with available hardware
"""

import msprime
import tskit
import numpy as np
import psutil
import time
import sys
from typing import Dict, List, Any

def get_memory_usage():
    """Get current memory usage in MB"""
    process = psutil.Process()
    return process.memory_info().rss / 1024 / 1024

def simulate_population_demographic_model(population: str, sample_size: int = 100) -> Dict[str, Any]:
    """
    Simulate demographic model for a specific population
    Based on established demographic parameters from literature
    """
    print(f"🧬 Simulating {population} demographic model...")
    
    start_memory = get_memory_usage()
    start_time = time.time()
    
    # Demographic parameters based on population genetics literature
    demographic_models = {
        "YRI": {
            "effective_population": 12300,
            "growth_rate": 0.0166,  # African expansion
            "mutation_rate": 1.25e-8,
            "recombination_rate": 1.13e-8,
            "description": "Yoruba - African expansion ~70kya"
        },
        "CEU": {
            "effective_population": 11418,
            "growth_rate": 0.0038,  # European expansion
            "mutation_rate": 1.25e-8,
            "recombination_rate": 1.13e-8,
            "description": "CEPH European - Out of Africa ~40kya"
        },
        "CHB": {
            "effective_population": 9025,
            "growth_rate": 0.0048,  # East Asian expansion
            "mutation_rate": 1.25e-8,
            "recombination_rate": 1.13e-8,
            "description": "Han Chinese - East Asian expansion ~50kya"
        }
    }
    
    if population not in demographic_models:
        population = "YRI"  # Default to YRI
    
    params = demographic_models[population]
    
    # Create demographic events (simplified model)
    demography = msprime.Demography()
    demography.add_population(name=population, initial_size=params["effective_population"])
    
    # Add population growth
    demography.add_population_parameters_change(
        time=0, initial_size=params["effective_population"], growth_rate=params["growth_rate"], population=population
    )
    
    # Ancient bottleneck (simplified)
    demography.add_population_parameters_change(
        time=2000, initial_size=1000, growth_rate=0, population=population  # ~70kya bottleneck
    )
    
    try:
        # Simulate chromosome 22 (smallest autosome for testing)
        sequence_length = 50_000_000  # 50Mb chromosome 22
        
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=params["recombination_rate"],
            random_seed=42
        )
        
        # Add mutations
        ts = msprime.sim_mutations(
            ts, 
            rate=params["mutation_rate"],
            random_seed=42
        )
        
        end_time = time.time()
        end_memory = get_memory_usage()
        
        # Extract basic statistics
        num_variants = ts.num_mutations
        num_trees = ts.num_trees
        num_samples = ts.num_samples
        
        result = {
            "population": population,
            "sample_size": sample_size,
            "sequence_length": sequence_length,
            "num_variants": num_variants,
            "num_trees": num_trees,
            "num_samples": num_samples,
            "simulation_time": end_time - start_time,
            "memory_used_mb": end_memory - start_memory,
            "peak_memory_mb": end_memory,
            "demographic_params": params,
            "tree_sequence": ts  # For further analysis
        }
        
        print(f"  ✅ {population}: {num_variants:,} variants, {num_trees:,} trees")
        print(f"     Memory: {result['memory_used_mb']:.1f}MB, Time: {result['simulation_time']:.1f}s")
        
        return result
        
    except Exception as e:
        print(f"  ❌ Error simulating {population}: {e}")
        return {"error": str(e), "population": population}

def test_memory_scaling():
    """Test memory requirements for different sample sizes"""
    print("\n📊 Testing Memory Scaling...")
    
    sample_sizes = [50, 100, 200, 500, 1000]
    results = []
    
    for sample_size in sample_sizes:
        print(f"\n🔬 Testing sample size: {sample_size}")
        result = simulate_population_demographic_model("YRI", sample_size)
        
        if "error" not in result:
            results.append({
                "sample_size": sample_size,
                "memory_mb": result["memory_used_mb"],
                "time_seconds": result["simulation_time"],
                "variants": result["num_variants"]
            })
            
            # Memory check
            if result["peak_memory_mb"] > 32000:  # 32GB threshold
                print(f"  ⚠️  High memory usage: {result['peak_memory_mb']:.1f}MB")
            
            if result["peak_memory_mb"] > 64000:  # 64GB threshold
                print(f"  🚨 Exceeds 64GB limit!")
                break
    
    return results

def analyze_variant_patterns(tree_sequence: tskit.TreeSequence) -> Dict[str, Any]:
    """Analyze variant patterns for authenticity validation"""
    print("🔍 Analyzing variant patterns...")
    
    # Calculate allele frequencies
    variants = []
    for variant in tree_sequence.variants():
        allele_count = np.sum(variant.genotypes)
        total_alleles = len(variant.genotypes)
        frequency = allele_count / total_alleles
        variants.append(frequency)
    
    variants = np.array(variants)
    
    # Site frequency spectrum
    freq_bins = np.linspace(0, 1, 21)  # 20 bins
    sfs, _ = np.histogram(variants, bins=freq_bins)
    
    # Basic statistics
    stats = {
        "total_variants": len(variants),
        "mean_frequency": np.mean(variants),
        "rare_variants": np.sum(variants < 0.05),  # MAF < 5%
        "common_variants": np.sum(variants > 0.05),
        "site_frequency_spectrum": sfs.tolist(),
        "diversity_pi": np.mean(variants * (1 - variants)) * 2  # Simplified pi
    }
    
    print(f"  📊 {stats['total_variants']:,} variants analyzed")
    print(f"     Rare variants (MAF<5%): {stats['rare_variants']:,} ({stats['rare_variants']/stats['total_variants']*100:.1f}%)")
    print(f"     Diversity (π): {stats['diversity_pi']:.6f}")
    
    return stats

def main():
    """Main testing function"""
    print("🧬 COALESCENT SIMULATION TESTING")
    print("=================================")
    print(f"💻 Available Memory: {psutil.virtual_memory().total / 1024**3:.1f} GB")
    print(f"🧠 Current Memory Usage: {get_memory_usage():.1f} MB")
    
    # Test single population simulation
    print("\n📍 Testing Single Population Simulation...")
    result_yri = simulate_population_demographic_model("YRI", 100)
    
    if "error" not in result_yri:
        # Analyze variant patterns
        variant_stats = analyze_variant_patterns(result_yri["tree_sequence"])
        
        # Test different populations
        print("\n🌍 Testing Multiple Populations...")
        populations = ["YRI", "CEU", "CHB"]
        population_results = []
        
        for pop in populations:
            result = simulate_population_demographic_model(pop, 100)
            if "error" not in result:
                population_results.append(result)
    
    # Memory scaling test
    scaling_results = test_memory_scaling()
    
    # Summary
    print("\n📋 TESTING SUMMARY")
    print("==================")
    print(f"✅ msprime version: {msprime.__version__}")
    print(f"✅ tskit version: {tskit.__version__}")
    
    if scaling_results:
        max_samples = max(r["sample_size"] for r in scaling_results)
        max_memory = max(r["memory_mb"] for r in scaling_results)
        print(f"✅ Maximum sample size tested: {max_samples}")
        print(f"✅ Peak memory usage: {max_memory:.1f} MB")
        
        # Feasibility assessment
        if max_memory < 32000:
            print("🎉 FEASIBLE: Well within 64GB memory limit")
        elif max_memory < 64000:
            print("⚠️  CAUTION: Approaching 64GB limit")
        else:
            print("🚨 PROBLEM: Exceeds 64GB memory limit")
    
    print("\n🚀 Ready for Phase 1 implementation!")

if __name__ == "__main__":
    main()