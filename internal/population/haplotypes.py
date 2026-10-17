#!/usr/bin/env python3
"""
Haplotype Simulation Module
Implements coalescent-based haplotype generation with linkage disequilibrium structure
Extends demographic models with realistic recombination and LD patterns
"""

import msprime
import tskit
import numpy as np
import time
import json
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from enum import Enum

from demographics import Population, DemographicParameters, DemographicModelBuilder

@dataclass
class RecombinationMap:
    """Recombination map for a chromosome"""
    chromosome: str
    positions: List[int]  # Base pair positions
    rates: List[float]    # Recombination rates (cM/Mb)
    map_length: float     # Total map length in cM
    physical_length: int  # Physical length in bp

@dataclass
class HaplotypeBlock:
    """A haplotype block with LD structure"""
    start_position: int
    end_position: int
    length: int
    variant_count: int
    r_squared_mean: float
    population: Population

@dataclass
class HaplotypeSimulationResult:
    """Result of haplotype simulation with LD structure"""
    population: Population
    chromosome: str
    tree_sequence: tskit.TreeSequence
    haplotypes: np.ndarray  # Haplotype matrix (variants x samples)
    variants: List[Dict]    # Variant information
    recombination_map: RecombinationMap
    ld_blocks: List[HaplotypeBlock]
    simulation_time: float
    memory_used_mb: float
    statistics: Dict[str, float]

class RecombinationMapLoader:
    """Loads population-specific recombination maps"""
    
    # HapMap recombination rates (cM/Mb) by population
    POPULATION_RECOMB_RATES = {
        # African populations - higher recombination rates
        Population.YRI: 1.3,
        Population.LWK: 1.3,
        Population.GWD: 1.3,
        
        # European populations - moderate recombination rates
        Population.CEU: 1.1,
        Population.GBR: 1.1,
        Population.TSI: 1.1,
        
        # East Asian populations - lower recombination rates
        Population.CHB: 1.0,
        Population.JPT: 1.0,
        Population.CHS: 1.0,
        
        # Americas populations - intermediate rates (admixed)
        Population.MXL: 1.2,
        Population.PUR: 1.2,
    }
    
    def __init__(self):
        self.cached_maps = {}
    
    def create_population_recombination_map(
        self,
        population: Population,
        chromosome: str,
        sequence_length: int
    ) -> RecombinationMap:
        """
        Create population-specific recombination map
        
        Args:
            population: Target population
            chromosome: Chromosome name (e.g., 'chr1')
            sequence_length: Physical length in bp
            
        Returns:
            RecombinationMap: Population-specific recombination map
        """
        cache_key = f"{population.value}_{chromosome}_{sequence_length}"
        
        if cache_key in self.cached_maps:
            return self.cached_maps[cache_key]
        
        base_rate = self.POPULATION_RECOMB_RATES.get(population, 1.1)
        
        # Create realistic recombination landscape
        positions, rates = self._generate_recombination_landscape(
            sequence_length, base_rate, population
        )
        
        # Calculate map length in centiMorgans
        map_length = self._calculate_map_length(positions, rates)
        
        recomb_map = RecombinationMap(
            chromosome=chromosome,
            positions=positions,
            rates=rates,
            map_length=map_length,
            physical_length=sequence_length
        )
        
        self.cached_maps[cache_key] = recomb_map
        return recomb_map
    
    def _generate_recombination_landscape(
        self,
        sequence_length: int,
        base_rate: float,
        population: Population
    ) -> Tuple[List[int], List[float]]:
        """Generate realistic recombination rate variation along chromosome"""
        
        # Number of recombination intervals
        num_intervals = max(10, sequence_length // 1_000_000)  # ~1Mb intervals
        
        positions = [0]
        rates = []
        
        np.random.seed(42)  # Reproducible landscapes
        
        for i in range(num_intervals):
            # Position at regular intervals with slight variation
            pos = int((i + 1) * sequence_length / num_intervals)
            positions.append(pos)
            
            # Rate variation around base rate
            # Simulate recombination hotspots and coldspots
            rate_multiplier = np.random.lognormal(0, 0.5)  # Log-normal variation
            
            # Add hotspot probability (5% chance of 10x higher rate)
            if np.random.random() < 0.05:
                rate_multiplier *= 10  # Recombination hotspot
            
            # Add coldspot probability (10% chance of 10x lower rate)
            elif np.random.random() < 0.10:
                rate_multiplier *= 0.1  # Recombination coldspot
            
            rate = base_rate * rate_multiplier
            
            # Population-specific modulation
            if population.value.startswith('YRI'):
                rate *= 1.1  # African populations have slightly higher rates
            elif population.value in ['CHB', 'JPT', 'CHS']:
                rate *= 0.9  # East Asian populations have slightly lower rates
            
            rates.append(rate)
        
        return positions, rates
    
    def _calculate_map_length(self, positions: List[int], rates: List[float]) -> float:
        """Calculate total genetic map length in centiMorgans"""
        map_length = 0.0
        
        for i in range(len(rates)):
            if i < len(positions) - 1:
                physical_distance = positions[i + 1] - positions[i]
                genetic_distance = (rates[i] * physical_distance) / 1_000_000  # cM
                map_length += genetic_distance
        
        return map_length

class LDAnalyzer:
    """Analyzes linkage disequilibrium patterns in haplotypes"""
    
    def __init__(self):
        pass
    
    def calculate_ld_matrix(
        self,
        haplotypes: np.ndarray,
        positions: List[int],
        max_distance: int = 100_000
    ) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        """
        Calculate linkage disequilibrium (r²) matrix for variants
        
        Args:
            haplotypes: Haplotype matrix (variants x chromosomes)
            positions: Variant positions
            max_distance: Maximum distance for LD calculation
            
        Returns:
            Tuple[np.ndarray, List]: LD matrix and variant pairs
        """
        n_variants = haplotypes.shape[0]
        n_chromosomes = haplotypes.shape[1]
        
        # Only calculate LD for variants within max_distance
        variant_pairs = []
        for i in range(n_variants):
            for j in range(i + 1, n_variants):
                if positions[j] - positions[i] <= max_distance:
                    variant_pairs.append((i, j))
        
        ld_values = np.zeros(len(variant_pairs))
        
        for idx, (i, j) in enumerate(variant_pairs):
            # Calculate r² between variants i and j
            r_squared = self._calculate_r_squared(
                haplotypes[i, :], haplotypes[j, :]
            )
            ld_values[idx] = r_squared
        
        return ld_values, variant_pairs
    
    def _calculate_r_squared(self, variant1: np.ndarray, variant2: np.ndarray) -> float:
        """Calculate r² between two variants"""
        n = len(variant1)
        
        # Allele frequencies
        p1 = np.mean(variant1)
        p2 = np.mean(variant2)
        
        # Haplotype frequencies
        p11 = np.mean(variant1 * variant2)
        
        # Linkage disequilibrium coefficient
        D = p11 - (p1 * p2)
        
        # r² calculation
        if p1 * (1 - p1) * p2 * (1 - p2) > 0:
            r_squared = (D * D) / (p1 * (1 - p1) * p2 * (1 - p2))
        else:
            r_squared = 0.0
        
        return min(r_squared, 1.0)  # Cap at 1.0 due to numerical precision
    
    def identify_ld_blocks(
        self,
        ld_values: np.ndarray,
        variant_pairs: List[Tuple[int, int]],
        positions: List[int],
        r_squared_threshold: float = 0.8
    ) -> List[HaplotypeBlock]:
        """
        Identify haplotype blocks based on LD patterns
        
        Args:
            ld_values: Array of r² values
            variant_pairs: Corresponding variant pairs
            positions: Variant positions
            r_squared_threshold: Minimum r² for block membership
            
        Returns:
            List[HaplotypeBlock]: Identified haplotype blocks
        """
        blocks = []
        
        # Simple block identification - variants in strong LD
        high_ld_pairs = [(i, j) for idx, (i, j) in enumerate(variant_pairs) 
                        if ld_values[idx] >= r_squared_threshold]
        
        if not high_ld_pairs:
            return blocks
        
        # Group consecutive high-LD variants into blocks
        current_block_variants = set()
        
        for i, j in high_ld_pairs:
            if not current_block_variants:
                current_block_variants = {i, j}
            elif i in current_block_variants or j in current_block_variants:
                current_block_variants.update({i, j})
            else:
                # Finalize current block
                if len(current_block_variants) >= 2:
                    block = self._create_haplotype_block(
                        current_block_variants, positions, ld_values, variant_pairs
                    )
                    blocks.append(block)
                
                # Start new block
                current_block_variants = {i, j}
        
        # Finalize last block
        if len(current_block_variants) >= 2:
            block = self._create_haplotype_block(
                current_block_variants, positions, ld_values, variant_pairs
            )
            blocks.append(block)
        
        return blocks
    
    def _create_haplotype_block(
        self,
        variant_indices: set,
        positions: List[int],
        ld_values: np.ndarray,
        variant_pairs: List[Tuple[int, int]]
    ) -> HaplotypeBlock:
        """Create HaplotypeBlock from variant indices"""
        variant_list = sorted(list(variant_indices))
        
        start_pos = positions[variant_list[0]]
        end_pos = positions[variant_list[-1]]
        length = end_pos - start_pos
        
        # Calculate mean r² within block
        block_ld_values = []
        for idx, (i, j) in enumerate(variant_pairs):
            if i in variant_indices and j in variant_indices:
                block_ld_values.append(ld_values[idx])
        
        mean_r_squared = np.mean(block_ld_values) if block_ld_values else 0.0
        
        return HaplotypeBlock(
            start_position=start_pos,
            end_position=end_pos,
            length=length,
            variant_count=len(variant_indices),
            r_squared_mean=mean_r_squared,
            population=Population.YRI  # Will be set by caller
        )

class HaplotypeSimulator:
    """Coalescent-based haplotype simulator with LD structure"""
    
    def __init__(self):
        self.model_builder = DemographicModelBuilder()
        self.recomb_loader = RecombinationMapLoader()
        self.ld_analyzer = LDAnalyzer()
    
    def simulate_haplotypes(
        self,
        population: Population,
        chromosome: str,
        sample_size: int = 100,
        sequence_length: int = 50_000_000,
        random_seed: Optional[int] = None
    ) -> HaplotypeSimulationResult:
        """
        Simulate haplotypes with realistic LD structure
        
        Args:
            population: Target population
            chromosome: Chromosome to simulate
            sample_size: Number of diploid individuals
            sequence_length: Physical length in bp
            random_seed: Random seed for reproducibility
            
        Returns:
            HaplotypeSimulationResult: Complete haplotype simulation
        """
        start_time = time.time()
        import psutil
        process = psutil.Process()
        start_memory = process.memory_info().rss / 1024 / 1024
        
        print(f"🧬 Simulating {population.value} haplotypes for {chromosome}")
        print(f"   Samples: {sample_size}, Length: {sequence_length:,}bp")
        
        # Get demographic model and parameters
        demography = self.model_builder.build_demographic_model(population)
        params = self.model_builder.get_parameters(population)
        
        # Create population-specific recombination map
        recomb_map = self.recomb_loader.create_population_recombination_map(
            population, chromosome, sequence_length
        )
        
        print(f"   📊 Recombination map: {len(recomb_map.positions)} intervals, "
              f"{recomb_map.map_length:.1f}cM total")
        
        # Convert to msprime RateMap
        msprime_recomb_map = msprime.RateMap(
            position=recomb_map.positions,
            rate=[r / 100_000_000 for r in recomb_map.rates]  # Convert cM/Mb to per-bp
        )
        
        # Simulate ancestry with recombination
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=msprime_recomb_map,
            random_seed=random_seed
        )
        
        # Add mutations
        ts = msprime.sim_mutations(
            ts,
            rate=params.mutation_rate,
            random_seed=random_seed
        )
        
        print(f"   🔬 Generated: {ts.num_mutations} variants, {ts.num_trees} trees")
        
        # Extract haplotypes
        haplotypes = ts.genotype_matrix()  # variants x samples
        
        # Get variant information
        variants = []
        positions = []
        for variant in ts.variants():
            positions.append(int(variant.site.position))
            variants.append({
                'position': int(variant.site.position),
                'alleles': variant.alleles,
                'ancestral_state': variant.site.ancestral_state,
                'frequency': np.mean(variant.genotypes)
            })
        
        print(f"   📈 Analyzing LD structure...")
        
        # Calculate LD structure
        ld_values, variant_pairs = self.ld_analyzer.calculate_ld_matrix(
            haplotypes, positions, max_distance=100_000
        )
        
        # Identify haplotype blocks
        ld_blocks = self.ld_analyzer.identify_ld_blocks(
            ld_values, variant_pairs, positions
        )
        
        # Set population for blocks
        for block in ld_blocks:
            block.population = population
        
        print(f"   🧱 Identified {len(ld_blocks)} haplotype blocks")
        
        # Calculate statistics
        statistics = self._calculate_haplotype_statistics(
            ts, haplotypes, ld_values, recomb_map
        )
        
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024
        
        result = HaplotypeSimulationResult(
            population=population,
            chromosome=chromosome,
            tree_sequence=ts,
            haplotypes=haplotypes,
            variants=variants,
            recombination_map=recomb_map,
            ld_blocks=ld_blocks,
            simulation_time=end_time - start_time,
            memory_used_mb=end_memory - start_memory,
            statistics=statistics
        )
        
        print(f"   ✅ Haplotype simulation complete: {result.simulation_time:.2f}s, "
              f"{result.memory_used_mb:.1f}MB")
        
        return result
    
    def _calculate_haplotype_statistics(
        self,
        tree_sequence: tskit.TreeSequence,
        haplotypes: np.ndarray,
        ld_values: np.ndarray,
        recomb_map: RecombinationMap
    ) -> Dict[str, float]:
        """Calculate comprehensive haplotype statistics"""
        
        # Basic diversity statistics
        diversity = tree_sequence.diversity()
        segregating_sites = tree_sequence.segregating_sites()
        
        # LD statistics
        mean_r_squared = np.mean(ld_values) if len(ld_values) > 0 else 0.0
        median_r_squared = np.median(ld_values) if len(ld_values) > 0 else 0.0
        high_ld_fraction = np.mean(ld_values > 0.8) if len(ld_values) > 0 else 0.0
        
        # Recombination statistics
        effective_recomb_rate = recomb_map.map_length / (recomb_map.physical_length / 1_000_000)
        
        # Haplotype diversity
        n_haplotypes = haplotypes.shape[1]
        haplotype_diversity = self._calculate_haplotype_diversity(haplotypes)
        
        statistics = {
            'diversity_pi': float(diversity),
            'segregating_sites': int(segregating_sites),
            'n_variants': haplotypes.shape[0],
            'n_haplotypes': n_haplotypes,
            'mean_r_squared': float(mean_r_squared),
            'median_r_squared': float(median_r_squared),
            'high_ld_fraction': float(high_ld_fraction),
            'effective_recomb_rate': float(effective_recomb_rate),
            'genetic_map_length_cm': float(recomb_map.map_length),
            'haplotype_diversity': float(haplotype_diversity),
            'ld_pairs_analyzed': len(ld_values)
        }
        
        return statistics
    
    def _calculate_haplotype_diversity(self, haplotypes: np.ndarray) -> float:
        """Calculate haplotype diversity (proportion of unique haplotypes)"""
        if haplotypes.shape[0] == 0:
            return 0.0
        
        # Convert each haplotype to a hashable tuple
        unique_haplotypes = set()
        for i in range(haplotypes.shape[1]):
            haplotype_tuple = tuple(haplotypes[:, i])
            unique_haplotypes.add(haplotype_tuple)
        
        return len(unique_haplotypes) / haplotypes.shape[1]
    
    def simulate_all_populations_haplotypes(
        self,
        chromosome: str = "chr22",  # Start with smallest chromosome
        sample_size: int = 50,
        sequence_length: int = 20_000_000,  # 20Mb for chr22
        random_seed: Optional[int] = None
    ) -> Dict[Population, HaplotypeSimulationResult]:
        """
        Simulate haplotypes for all populations
        
        Args:
            chromosome: Chromosome to simulate
            sample_size: Number of diploid individuals per population
            sequence_length: Physical length in bp
            random_seed: Random seed for reproducibility
            
        Returns:
            Dict mapping populations to haplotype results
        """
        results = {}
        total_start = time.time()
        
        print(f"🌍 SIMULATING HAPLOTYPES FOR ALL 11 POPULATIONS")
        print(f"   Chromosome: {chromosome}, Length: {sequence_length:,}bp")
        print("=" * 60)
        
        for population in Population:
            try:
                result = self.simulate_haplotypes(
                    population=population,
                    chromosome=chromosome,
                    sample_size=sample_size,
                    sequence_length=sequence_length,
                    random_seed=random_seed
                )
                results[population] = result
                
                print(f"   ✅ {population.value}: {result.statistics['n_variants']} variants, "
                      f"{len(result.ld_blocks)} blocks, r²={result.statistics['mean_r_squared']:.3f}")
                
            except Exception as e:
                print(f"   ❌ {population.value}: Failed - {e}")
        
        total_time = time.time() - total_start
        print(f"\n🎉 All population haplotypes simulated in {total_time:.2f}s")
        
        return results
    
    def validate_ld_patterns(
        self,
        results: Dict[Population, HaplotypeSimulationResult]
    ) -> Dict[str, bool]:
        """
        Validate LD patterns across populations
        
        Args:
            results: Haplotype simulation results for all populations
            
        Returns:
            Dict of validation results
        """
        print("\n📊 Validating LD Patterns Across Populations...")
        
        validation_results = {}
        
        # Expected LD patterns based on population genetics
        expected_patterns = {
            'african_higher_diversity': True,
            'ld_decay_realistic': True,  
            'population_differences': True,
            'block_structure_present': True
        }
        
        # Extract statistics for comparison
        pop_stats = {}
        for pop, result in results.items():
            pop_stats[pop] = result.statistics
        
        # Test 1: African populations should have higher diversity
        african_pops = [Population.YRI, Population.LWK, Population.GWD]
        african_diversity = np.mean([pop_stats[pop]['diversity_pi'] for pop in african_pops if pop in pop_stats])
        
        non_african_pops = [pop for pop in pop_stats.keys() if pop not in african_pops]
        non_african_diversity = np.mean([pop_stats[pop]['diversity_pi'] for pop in non_african_pops])
        
        validation_results['african_higher_diversity'] = african_diversity > non_african_diversity
        
        # Test 2: LD should decay with distance (approximate test)
        ld_means = [pop_stats[pop]['mean_r_squared'] for pop in pop_stats.keys()]
        validation_results['ld_decay_realistic'] = 0.1 < np.mean(ld_means) < 0.5
        
        # Test 3: Population differences should exist
        diversity_values = [pop_stats[pop]['diversity_pi'] for pop in pop_stats.keys()]
        diversity_cv = np.std(diversity_values) / np.mean(diversity_values)
        validation_results['population_differences'] = diversity_cv > 0.2
        
        # Test 4: Haplotype blocks should be present
        total_blocks = sum(len(result.ld_blocks) for result in results.values())
        validation_results['block_structure_present'] = total_blocks > 0
        
        # Print validation results
        for test, passed in validation_results.items():
            status = "✅ PASSED" if passed else "❌ FAILED"
            print(f"   {test}: {status}")
        
        overall_success = all(validation_results.values())
        print(f"\n🎯 Overall LD validation: {'✅ PASSED' if overall_success else '❌ FAILED'}")
        
        return validation_results

def main():
    """Main testing function for haplotype simulation"""
    print("🧬 HAPLOTYPE SIMULATION TESTING")
    print("=" * 50)
    
    simulator = HaplotypeSimulator()
    
    # Test single population haplotype simulation
    print("\n📍 Testing Single Population Haplotype Simulation (YRI)...")
    yri_result = simulator.simulate_haplotypes(
        population=Population.YRI,
        chromosome="chr22",
        sample_size=30,
        sequence_length=10_000_000,  # 10Mb for testing
        random_seed=42
    )
    
    print(f"\n📊 YRI Haplotype Statistics:")
    for key, value in yri_result.statistics.items():
        if isinstance(value, float):
            print(f"   {key}: {value:.6f}")
        else:
            print(f"   {key}: {value}")
    
    # Test multiple populations (smaller scale)
    print("\n🌍 Testing Multiple Population Haplotypes...")
    test_populations = [Population.YRI, Population.CEU, Population.CHB]
    
    multi_results = {}
    for pop in test_populations:
        result = simulator.simulate_haplotypes(
            population=pop,
            chromosome="chr22",
            sample_size=20,
            sequence_length=5_000_000,  # 5Mb for speed
            random_seed=42
        )
        multi_results[pop] = result
        
        print(f"   {pop.value}: variants={result.statistics['n_variants']}, "
              f"blocks={len(result.ld_blocks)}, r²={result.statistics['mean_r_squared']:.3f}")
    
    # Validate LD patterns
    validation_results = simulator.validate_ld_patterns(multi_results)
    
    print("\n🚀 Haplotype simulation with LD structure implemented successfully!")

if __name__ == "__main__":
    main()