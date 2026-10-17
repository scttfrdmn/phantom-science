#!/usr/bin/env python3
"""
Population Demographics Module
Implements 11-population demographic models using msprime coalescent simulation
Integrates with existing Go-based authenticity engine
"""

import msprime
import tskit
import numpy as np
import time
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from enum import Enum

class Population(Enum):
    """Enumeration of supported populations"""
    YRI = "YRI"  # Yoruba - West Africa
    LWK = "LWK"  # Luhya - East Africa  
    GWD = "GWD"  # Gambian - West Africa
    CEU = "CEU"  # CEPH European - Northern Europe
    GBR = "GBR"  # British - Northern Europe
    TSI = "TSI"  # Tuscan - Southern Europe
    CHB = "CHB"  # Han Chinese - East Asia
    JPT = "JPT"  # Japanese - East Asia
    CHS = "CHS"  # Southern Han Chinese - East Asia
    MXL = "MXL"  # Mexican - Americas
    PUR = "PUR"  # Puerto Rican - Americas

@dataclass
class DemographicParameters:
    """Demographic parameters for a population"""
    population: Population
    effective_population: int
    growth_rate: float
    mutation_rate: float
    recombination_rate: float
    bottleneck_time: int  # generations ago
    bottleneck_size: int
    founding_time_years: int
    founding_time_generations: int
    super_population: str
    region: str
    description: str

class DemographicModelBuilder:
    """Builder for msprime demographic models"""
    
    # Standard rates from population genetics literature
    MUTATION_RATE = 1.25e-8  # per base per generation
    RECOMBINATION_RATE = 1.13e-8  # per base per generation
    GENERATION_TIME = 25  # years per generation
    
    # Population parameters from demographic inference studies
    PARAMETERS = {
        Population.YRI: DemographicParameters(
            population=Population.YRI,
            effective_population=12300,
            growth_rate=0.0166,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=2000,  # ~70kya bottleneck
            bottleneck_size=1000,
            founding_time_years=70000,
            founding_time_generations=70000 // GENERATION_TIME,
            super_population="African",
            region="West Africa",
            description="Yoruba in Ibadan, Nigeria"
        ),
        Population.LWK: DemographicParameters(
            population=Population.LWK,
            effective_population=10000,
            growth_rate=0.0150,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=2000,
            bottleneck_size=1000,
            founding_time_years=70000,
            founding_time_generations=70000 // GENERATION_TIME,
            super_population="African",
            region="East Africa",
            description="Luhya in Webuye, Kenya"
        ),
        Population.GWD: DemographicParameters(
            population=Population.GWD,
            effective_population=11000,
            growth_rate=0.0160,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=2000,
            bottleneck_size=1000,
            founding_time_years=70000,
            founding_time_generations=70000 // GENERATION_TIME,
            super_population="African",
            region="West Africa",
            description="Gambian in Western Divisions"
        ),
        Population.CEU: DemographicParameters(
            population=Population.CEU,
            effective_population=11418,
            growth_rate=0.0038,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1600,  # ~56kya bottleneck
            bottleneck_size=2000,  # Less severe European bottleneck
            founding_time_years=40000,
            founding_time_generations=40000 // GENERATION_TIME,
            super_population="European",
            region="Northern Europe",
            description="Utah residents with Northern European ancestry"
        ),
        Population.GBR: DemographicParameters(
            population=Population.GBR,
            effective_population=11000,
            growth_rate=0.0040,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1600,
            bottleneck_size=2000,
            founding_time_years=40000,
            founding_time_generations=40000 // GENERATION_TIME,
            super_population="European",
            region="Northern Europe",
            description="British in England and Scotland"
        ),
        Population.TSI: DemographicParameters(
            population=Population.TSI,
            effective_population=10500,
            growth_rate=0.0035,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1600,
            bottleneck_size=2000,
            founding_time_years=40000,
            founding_time_generations=40000 // GENERATION_TIME,
            super_population="European",
            region="Southern Europe",
            description="Toscani in Italia"
        ),
        Population.CHB: DemographicParameters(
            population=Population.CHB,
            effective_population=9025,
            growth_rate=0.0048,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1800,  # ~63kya bottleneck 
            bottleneck_size=1500,  # Severe East Asian bottleneck
            founding_time_years=50000,
            founding_time_generations=50000 // GENERATION_TIME,
            super_population="East Asian",
            region="East Asia",
            description="Han Chinese in Beijing, China"
        ),
        Population.JPT: DemographicParameters(
            population=Population.JPT,
            effective_population=9000,
            growth_rate=0.0050,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1800,
            bottleneck_size=1500,
            founding_time_years=50000,
            founding_time_generations=50000 // GENERATION_TIME,
            super_population="East Asian",
            region="East Asia",
            description="Japanese in Tokyo, Japan"
        ),
        Population.CHS: DemographicParameters(
            population=Population.CHS,
            effective_population=8500,
            growth_rate=0.0052,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=1800,
            bottleneck_size=1500,
            founding_time_years=50000,
            founding_time_generations=50000 // GENERATION_TIME,
            super_population="East Asian",
            region="East Asia",
            description="Southern Han Chinese"
        ),
        Population.MXL: DemographicParameters(
            population=Population.MXL,
            effective_population=7000,
            growth_rate=0.0080,  # Rapid expansion
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=500,  # ~17.5kya bottleneck
            bottleneck_size=800,  # Severe Americas bottleneck
            founding_time_years=15000,
            founding_time_generations=15000 // GENERATION_TIME,
            super_population="Americas",
            region="North America",
            description="Mexican Ancestry from Los Angeles USA"
        ),
        Population.PUR: DemographicParameters(
            population=Population.PUR,
            effective_population=6500,
            growth_rate=0.0090,
            mutation_rate=MUTATION_RATE,
            recombination_rate=RECOMBINATION_RATE,
            bottleneck_time=500,
            bottleneck_size=800,
            founding_time_years=15000,
            founding_time_generations=15000 // GENERATION_TIME,
            super_population="Americas",
            region="Caribbean",
            description="Puerto Ricans from Puerto Rico"
        )
    }
    
    def build_demographic_model(self, population: Population) -> msprime.Demography:
        """
        Build demographic model for specified population
        
        Args:
            population: Target population
            
        Returns:
            msprime.Demography: Configured demographic model
        """
        params = self.PARAMETERS[population]
        
        # Create demography
        demography = msprime.Demography()
        
        # Add population with current parameters
        demography.add_population(
            name=population.value,
            initial_size=params.effective_population
        )
        
        # Add exponential growth from present
        demography.add_population_parameters_change(
            time=0,
            initial_size=params.effective_population,
            growth_rate=params.growth_rate,
            population=population.value
        )
        
        # Add historical bottleneck
        demography.add_population_parameters_change(
            time=params.bottleneck_time,
            initial_size=params.bottleneck_size,
            growth_rate=0.0,  # No growth during bottleneck
            population=population.value
        )
        
        # Add ancestral population size before bottleneck
        demography.add_population_parameters_change(
            time=params.bottleneck_time + 100,  # Ancestral period
            initial_size=params.bottleneck_size * 2,  # Larger ancestral population
            growth_rate=0.0,
            population=population.value
        )
        
        return demography
    
    def get_parameters(self, population: Population) -> DemographicParameters:
        """Get demographic parameters for population"""
        return self.PARAMETERS[population]

@dataclass 
class SimulationResult:
    """Result of demographic simulation"""
    population: Population
    tree_sequence: tskit.TreeSequence
    num_variants: int
    num_trees: int
    num_samples: int
    sequence_length: int
    simulation_time: float
    memory_used_mb: float
    demographic_params: DemographicParameters

class PopulationSimulator:
    """Population genetics simulator using coalescent theory"""
    
    def __init__(self):
        self.model_builder = DemographicModelBuilder()
    
    def simulate_population(
        self,
        population: Population,
        sample_size: int = 100,
        sequence_length: int = 50_000_000,  # 50Mb chromosome
        random_seed: Optional[int] = None
    ) -> SimulationResult:
        """
        Simulate population using coalescent model
        
        Args:
            population: Target population
            sample_size: Number of diploid individuals to sample
            sequence_length: Length of sequence to simulate (bp)
            random_seed: Random seed for reproducibility
            
        Returns:
            SimulationResult: Simulation results with tree sequence
        """
        start_time = time.time()
        import psutil
        process = psutil.Process()
        start_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Get demographic model and parameters
        demography = self.model_builder.build_demographic_model(population)
        params = self.model_builder.get_parameters(population)
        
        print(f"🧬 Simulating {population.value} ({params.description})")
        print(f"   Samples: {sample_size}, Sequence: {sequence_length:,}bp")
        
        # Simulate ancestry
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=params.recombination_rate,
            random_seed=random_seed
        )
        
        # Add mutations
        ts = msprime.sim_mutations(
            ts,
            rate=params.mutation_rate,
            random_seed=random_seed
        )
        
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        result = SimulationResult(
            population=population,
            tree_sequence=ts,
            num_variants=ts.num_mutations,
            num_trees=ts.num_trees,
            num_samples=ts.num_samples,
            sequence_length=sequence_length,
            simulation_time=end_time - start_time,
            memory_used_mb=end_memory - start_memory,
            demographic_params=params
        )
        
        print(f"   ✅ Generated: {result.num_variants:,} variants, {result.num_trees:,} trees")
        print(f"   📊 Time: {result.simulation_time:.2f}s, Memory: {result.memory_used_mb:.1f}MB")
        
        return result
    
    def simulate_all_populations(
        self,
        sample_size: int = 100,
        sequence_length: int = 50_000_000,
        random_seed: Optional[int] = None
    ) -> Dict[Population, SimulationResult]:
        """
        Simulate all 11 populations
        
        Args:
            sample_size: Number of diploid individuals per population
            sequence_length: Length of sequence to simulate (bp)
            random_seed: Random seed for reproducibility
            
        Returns:
            Dict[Population, SimulationResult]: Results for all populations
        """
        results = {}
        total_start = time.time()
        
        print("🌍 SIMULATING ALL 11 POPULATIONS")
        print("=" * 50)
        
        for population in Population:
            result = self.simulate_population(
                population=population,
                sample_size=sample_size,
                sequence_length=sequence_length,
                random_seed=random_seed
            )
            results[population] = result
            print()
        
        total_time = time.time() - total_start
        print(f"🎉 All populations simulated in {total_time:.2f}s")
        
        return results
    
    def extract_variants_for_go_integration(
        self,
        result: SimulationResult
    ) -> List[Dict]:
        """
        Extract variant data for Go authenticity engine integration
        
        Args:
            result: Simulation result
            
        Returns:
            List[Dict]: Variant data for Go integration
        """
        variants = []
        
        for variant in result.tree_sequence.variants():
            # Calculate allele frequency
            allele_count = np.sum(variant.genotypes)
            total_alleles = len(variant.genotypes)
            frequency = allele_count / total_alleles
            
            variant_data = {
                'position': int(variant.site.position),
                'alleles': variant.alleles,
                'genotypes': variant.genotypes.tolist(),
                'frequency': frequency,
                'ancestral_state': variant.site.ancestral_state,
                'num_alleles': len(variant.alleles)
            }
            variants.append(variant_data)
        
        return variants
    
    def calculate_population_statistics(
        self,
        result: SimulationResult
    ) -> Dict:
        """
        Calculate population genetics statistics
        
        Args:
            result: Simulation result
            
        Returns:
            Dict: Population statistics
        """
        ts = result.tree_sequence
        
        # Basic diversity statistics
        diversity = ts.diversity()
        segregating_sites = ts.segregating_sites()
        
        # Allele frequency spectrum
        frequencies = []
        for variant in ts.variants():
            freq = np.sum(variant.genotypes) / len(variant.genotypes)
            frequencies.append(freq)
        
        frequencies = np.array(frequencies)
        
        # Site frequency spectrum bins
        sfs_bins = np.linspace(0, 1, 21)
        sfs, _ = np.histogram(frequencies, bins=sfs_bins)
        
        # Tajima's D (simplified calculation)
        n = ts.num_samples
        a1 = sum(1/i for i in range(1, n))
        a2 = sum(1/(i*i) for i in range(1, n))
        b1 = (n + 1) / (3 * (n - 1))
        b2 = 2 * (n*n + n + 3) / (9 * n * (n - 1))
        c1 = b1 - 1/a1
        c2 = b2 - (n + 2)/(a1 * n) + a2/(a1*a1)
        e1 = c1 / a1
        e2 = c2 / (a1*a1 + a2)
        
        theta_w = segregating_sites / a1
        theta_pi = diversity
        
        if e2 > 0:
            tajima_d = (theta_pi - theta_w) / np.sqrt(e1 * segregating_sites + e2 * segregating_sites * (segregating_sites - 1))
        else:
            tajima_d = 0.0
        
        stats = {
            'population': result.population.value,
            'num_variants': result.num_variants,
            'num_samples': result.num_samples,
            'diversity_pi': float(diversity),
            'watterson_theta': float(theta_w),
            'tajima_d': float(tajima_d),
            'segregating_sites': int(segregating_sites),
            'rare_variants_count': int(np.sum(frequencies < 0.05)),
            'common_variants_count': int(np.sum(frequencies >= 0.05)),
            'site_frequency_spectrum': sfs.tolist(),
            'mean_allele_frequency': float(np.mean(frequencies)),
            'demographic_params': {
                'effective_population': result.demographic_params.effective_population,
                'growth_rate': result.demographic_params.growth_rate,
                'super_population': result.demographic_params.super_population,
                'region': result.demographic_params.region
            }
        }
        
        return stats

def main():
    """Main testing function"""
    print("🧬 POPULATION DEMOGRAPHICS TESTING")
    print("=" * 50)
    
    simulator = PopulationSimulator()
    
    # Test single population
    print("\n📍 Testing YRI Population Model...")
    yri_result = simulator.simulate_population(
        population=Population.YRI,
        sample_size=100,
        sequence_length=10_000_000,  # 10Mb for testing
        random_seed=42
    )
    
    # Calculate statistics
    yri_stats = simulator.calculate_population_statistics(yri_result)
    print(f"\n📊 YRI Population Statistics:")
    print(f"   Diversity (π): {yri_stats['diversity_pi']:.6f}")
    print(f"   Watterson's θ: {yri_stats['watterson_theta']:.6f}")
    print(f"   Tajima's D: {yri_stats['tajima_d']:.3f}")
    print(f"   Rare variants: {yri_stats['rare_variants_count']:,}")
    
    # Test multiple populations (smaller scale)
    print("\n🌍 Testing Multiple Population Models...")
    test_populations = [Population.YRI, Population.CEU, Population.CHB]
    
    for pop in test_populations:
        result = simulator.simulate_population(
            population=pop,
            sample_size=50,
            sequence_length=5_000_000,  # 5Mb for speed
            random_seed=42
        )
        stats = simulator.calculate_population_statistics(result)
        print(f"   {pop.value}: π={stats['diversity_pi']:.6f}, D={stats['tajima_d']:.3f}")
    
    print("\n🚀 Demographic models implemented successfully!")

if __name__ == "__main__":
    main()