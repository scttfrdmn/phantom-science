#!/usr/bin/env python3
"""
Demographic Simulation CLI
Command-line interface for demographic models that can be called from Go
"""

import argparse
import json
import sys
from demographics import PopulationSimulator, Population

def main():
    parser = argparse.ArgumentParser(description='Run demographic population simulation')
    parser.add_argument('--population', required=True, choices=[p.value for p in Population],
                        help='Population to simulate')
    parser.add_argument('--sample-size', type=int, default=100,
                        help='Number of diploid individuals to sample')
    parser.add_argument('--sequence-length', type=int, default=50_000_000,
                        help='Length of sequence to simulate in base pairs')
    parser.add_argument('--random-seed', type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--output-format', choices=['json', 'summary'], default='summary',
                        help='Output format')
    parser.add_argument('--include-variants', action='store_true',
                        help='Include detailed variant data in JSON output')
    
    args = parser.parse_args()
    
    try:
        # Create simulator and run simulation
        simulator = PopulationSimulator()
        population = Population(args.population)
        
        result = simulator.simulate_population(
            population=population,
            sample_size=args.sample_size,
            sequence_length=args.sequence_length,
            random_seed=args.random_seed
        )
        
        # Calculate statistics
        stats = simulator.calculate_population_statistics(result)
        
        if args.output_format == 'json':
            # Prepare JSON output
            output_data = {
                'population': result.population.value,
                'num_variants': result.num_variants,
                'num_trees': result.num_trees,
                'num_samples': result.num_samples,
                'sequence_length': result.sequence_length,
                'simulation_time': result.simulation_time,
                'memory_used_mb': result.memory_used_mb,
                'demographic_params': {
                    'population': result.demographic_params.population.value,
                    'effective_population': result.demographic_params.effective_population,
                    'growth_rate': result.demographic_params.growth_rate,
                    'mutation_rate': result.demographic_params.mutation_rate,
                    'recombination_rate': result.demographic_params.recombination_rate,
                    'bottleneck_time': result.demographic_params.bottleneck_time,
                    'bottleneck_size': result.demographic_params.bottleneck_size,
                    'founding_time_years': result.demographic_params.founding_time_years,
                    'founding_time_generations': result.demographic_params.founding_time_generations,
                    'super_population': result.demographic_params.super_population,
                    'region': result.demographic_params.region,
                    'description': result.demographic_params.description
                },
                'statistics': stats
            }
            
            # Include variants if requested (can be large)
            if args.include_variants:
                output_data['variants'] = simulator.extract_variants_for_go_integration(result)
            
            print(json.dumps(output_data, indent=2))
            
        else:
            # Summary output
            print(f"🧬 {result.population.value} Demographic Simulation Results")
            print("=" * 50)
            print(f"Population: {result.demographic_params.description}")
            print(f"Super Population: {result.demographic_params.super_population}")
            print(f"Region: {result.demographic_params.region}")
            print()
            print(f"Simulation Parameters:")
            print(f"  Sample Size: {result.num_samples} individuals")
            print(f"  Sequence Length: {result.sequence_length:,} bp")
            print(f"  Effective Population: {result.demographic_params.effective_population:,}")
            print(f"  Growth Rate: {result.demographic_params.growth_rate:.4f}")
            print()
            print(f"Results:")
            print(f"  Variants Generated: {result.num_variants:,}")
            print(f"  Trees in Genealogy: {result.num_trees:,}")
            print(f"  Simulation Time: {result.simulation_time:.3f} seconds")
            print(f"  Memory Used: {result.memory_used_mb:.1f} MB")
            print()
            print(f"Population Genetics Statistics:")
            print(f"  Diversity (π): {stats['diversity_pi']:.6f}")
            print(f"  Watterson's θ: {stats['watterson_theta']:.6f}")
            print(f"  Tajima's D: {stats['tajima_d']:.3f}")
            print(f"  Segregating Sites: {stats['segregating_sites']:,}")
            print(f"  Rare Variants (MAF<5%): {stats['rare_variants_count']:,} ({stats['rare_variants_count']/stats['num_variants']*100:.1f}%)")
            print(f"  Mean Allele Frequency: {stats['mean_allele_frequency']:.6f}")
            
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()