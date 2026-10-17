#!/usr/bin/env python3
"""
Haplotype Simulation CLI
Command-line interface for haplotype simulation that can be called from Go
"""

import argparse
import json
import sys
import numpy as np
from haplotypes import HaplotypeSimulator, Population

def serialize_haplotype_result(result):
    """Convert HaplotypeSimulationResult to JSON-serializable format"""
    
    # Convert numpy arrays to lists for JSON serialization
    haplotype_matrix = None
    if hasattr(result, 'haplotypes') and result.haplotypes is not None:
        haplotype_matrix = result.haplotypes.tolist()
    
    return {
        'population': result.population.value,
        'chromosome': result.chromosome,
        'num_variants': len(result.variants),
        'num_trees': result.tree_sequence.num_trees,
        'num_samples': result.tree_sequence.num_samples,
        'sequence_length': result.tree_sequence.sequence_length,
        'simulation_time': result.simulation_time,
        'memory_used_mb': result.memory_used_mb,
        'variants': result.variants,
        'recombination_map': {
            'chromosome': result.recombination_map.chromosome,
            'positions': result.recombination_map.positions,
            'rates': result.recombination_map.rates,
            'map_length': result.recombination_map.map_length,
            'physical_length': result.recombination_map.physical_length
        },
        'ld_blocks': [
            {
                'start_position': block.start_position,
                'end_position': block.end_position,
                'length': block.length,
                'variant_count': block.variant_count,
                'r_squared_mean': block.r_squared_mean,
                'population': block.population.value
            }
            for block in result.ld_blocks
        ],
        'statistics': result.statistics,
        'haplotype_matrix': haplotype_matrix
    }

def main():
    parser = argparse.ArgumentParser(description='Run haplotype simulation with LD structure')
    parser.add_argument('--population', required=True, choices=[p.value for p in Population],
                        help='Population to simulate')
    parser.add_argument('--chromosome', default='chr22',
                        help='Chromosome to simulate (default: chr22)')
    parser.add_argument('--sample-size', type=int, default=50,
                        help='Number of diploid individuals to sample')
    parser.add_argument('--sequence-length', type=int, default=20_000_000,
                        help='Length of sequence to simulate in base pairs')
    parser.add_argument('--random-seed', type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--output-format', choices=['json', 'summary'], default='summary',
                        help='Output format')
    parser.add_argument('--include-haplotype-matrix', action='store_true',
                        help='Include full haplotype matrix in JSON output (warning: can be large)')
    
    args = parser.parse_args()
    
    try:
        # Create simulator and run haplotype simulation
        simulator = HaplotypeSimulator()
        population = Population(args.population)
        
        result = simulator.simulate_haplotypes(
            population=population,
            chromosome=args.chromosome,
            sample_size=args.sample_size,
            sequence_length=args.sequence_length,
            random_seed=args.random_seed
        )
        
        if args.output_format == 'json':
            # Prepare JSON output
            output_data = serialize_haplotype_result(result)
            
            # Remove haplotype matrix if not requested (saves memory and bandwidth)
            if not args.include_haplotype_matrix:
                output_data['haplotype_matrix'] = None
            
            print(json.dumps(output_data, indent=2))
            
        else:
            # Summary output
            print(f"🧬 {result.population.value} Haplotype Simulation Results")
            print("=" * 60)
            print(f"Population: {result.population.value}")
            print(f"Chromosome: {result.chromosome}")
            print()
            print(f"Simulation Parameters:")
            print(f"  Sample Size: {result.tree_sequence.num_samples} chromosomes ({args.sample_size} individuals)")
            print(f"  Sequence Length: {int(result.tree_sequence.sequence_length):,} bp")
            print()
            print(f"Results:")
            print(f"  Variants Generated: {len(result.variants):,}")
            print(f"  Trees in Genealogy: {result.tree_sequence.num_trees:,}")
            print(f"  Haplotype Blocks Identified: {len(result.ld_blocks):,}")
            print(f"  Simulation Time: {result.simulation_time:.3f} seconds")
            print(f"  Memory Used: {result.memory_used_mb:.1f} MB")
            print()
            print(f"Recombination Map:")
            print(f"  Genetic Map Length: {result.recombination_map.map_length:.1f} cM")
            print(f"  Effective Recombination Rate: {result.statistics['effective_recomb_rate']:.2f} cM/Mb")
            print(f"  Recombination Intervals: {len(result.recombination_map.positions):,}")
            print()
            print(f"Linkage Disequilibrium Statistics:")
            print(f"  Mean r²: {result.statistics['mean_r_squared']:.6f}")
            print(f"  Median r²: {result.statistics['median_r_squared']:.6f}")
            print(f"  High LD Fraction (r²>0.8): {result.statistics['high_ld_fraction']:.3f}")
            print(f"  LD Pairs Analyzed: {result.statistics['ld_pairs_analyzed']:,}")
            print()
            print(f"Population Genetics Statistics:")
            print(f"  Diversity (π): {result.statistics['diversity_pi']:.6f}")
            print(f"  Segregating Sites: {result.statistics['segregating_sites']:,}")
            print(f"  Haplotype Diversity: {result.statistics['haplotype_diversity']:.6f}")
            print()
            print(f"Haplotype Block Summary:")
            if result.ld_blocks:
                total_block_length = sum(block.length for block in result.ld_blocks)
                mean_block_size = total_block_length / len(result.ld_blocks)
                mean_block_r_squared = sum(block.r_squared_mean for block in result.ld_blocks) / len(result.ld_blocks)
                
                print(f"  Total Blocks: {len(result.ld_blocks):,}")
                print(f"  Mean Block Size: {mean_block_size:,.0f} bp")
                print(f"  Mean Block r²: {mean_block_r_squared:.3f}")
                print(f"  Genome Coverage: {total_block_length/result.tree_sequence.sequence_length*100:.1f}%")
            else:
                print(f"  No haplotype blocks identified")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()