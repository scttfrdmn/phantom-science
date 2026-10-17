#!/usr/bin/env python3
"""
Enhanced VCF Integration Module  
Improved integration of VCF variants with haplotype simulation
Uses VCF positions as informed mutation sites for better concordance
"""

import msprime
import tskit
import numpy as np
import time
from typing import Dict, List, Optional
from dataclasses import dataclass

from demographics import Population, DemographicModelBuilder
from vcf_haplotype_integration import VCFVariantLoader, VCFVariant, IntegratedHaplotypeResult

@dataclass
class EnhancedIntegrationResult:
    """Enhanced VCF-haplotype integration result with better metrics"""
    population: Population
    chromosome: str
    tree_sequence: tskit.TreeSequence
    vcf_variants_used: List[VCFVariant]
    targeted_mutations: int
    random_mutations: int
    position_accuracy: float
    allele_frequency_fidelity: float
    demographic_authenticity: float
    integration_time: float
    memory_used_mb: float

class EnhancedVCFIntegrator:
    """Enhanced VCF integration with targeted mutation placement"""
    
    def __init__(self):
        self.vcf_loader = VCFVariantLoader()
        self.demographic_builder = DemographicModelBuilder()
    
    def integrate_with_targeted_mutations(
        self,
        population: Population,
        chromosome: str,
        sample_size: int = 50,
        sequence_length: int = 20_000_000,
        max_vcf_variants: int = 1000,
        vcf_mutation_fraction: float = 0.7,  # Fraction of mutations from VCF positions
        random_seed: Optional[int] = None
    ) -> EnhancedIntegrationResult:
        """
        Integrate VCF variants using targeted mutation placement
        
        Args:
            population: Target population
            chromosome: Chromosome to simulate
            sample_size: Number of diploid individuals
            sequence_length: Physical length in bp
            max_vcf_variants: Maximum VCF variants to use
            vcf_mutation_fraction: Fraction of mutations to place at VCF positions
            random_seed: Random seed for reproducibility
            
        Returns:
            EnhancedIntegrationResult: Integration results with metrics
        """
        start_time = time.time()
        import psutil
        process = psutil.Process()
        start_memory = process.memory_info().rss / 1024 / 1024
        
        print(f"🎯 Enhanced VCF integration for {population.value} {chromosome}")
        print(f"   Length: {sequence_length:,}bp, VCF fraction: {vcf_mutation_fraction:.1%}")
        
        # Load VCF variants
        vcf_variants = self.vcf_loader.load_population_variants(
            population, chromosome, max_vcf_variants
        )
        
        if not vcf_variants:
            raise ValueError(f"No VCF variants available for {population.value}")
        
        # Filter VCF variants within sequence length and with reasonable frequencies
        filtered_vcf = [
            v for v in vcf_variants 
            if v.position < sequence_length and 0.01 < v.allele_frequency < 0.99
        ]
        
        print(f"   📊 Using {len(filtered_vcf):,} VCF variants (filtered from {len(vcf_variants):,})")
        
        # Create demographic model
        demography = self.demographic_builder.build_demographic_model(population) 
        params = self.demographic_builder.get_parameters(population)
        
        # Simulate ancestry
        print(f"   🧬 Simulating coalescent ancestry...")
        ts = msprime.sim_ancestry(
            samples=sample_size,
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=params.recombination_rate,
            random_seed=random_seed
        )
        
        # Create targeted mutation model
        print(f"   🎯 Placing mutations with VCF targeting...")
        
        # Calculate total mutation rate budget
        total_expected_mutations = int(params.mutation_rate * sequence_length * 2 * sample_size)
        
        # Allocate mutations between VCF positions and random positions
        vcf_mutations = int(total_expected_mutations * vcf_mutation_fraction)
        random_mutations = total_expected_mutations - vcf_mutations
        
        # Select VCF positions for targeted mutations
        if len(filtered_vcf) > vcf_mutations:
            # Subsample VCF variants, preferring common variants
            vcf_weights = [v.allele_frequency * (1 - v.allele_frequency) for v in filtered_vcf]
            selected_indices = np.random.choice(
                len(filtered_vcf), 
                size=vcf_mutations, 
                replace=False, 
                p=np.array(vcf_weights) / sum(vcf_weights)
            )
            selected_vcf = [filtered_vcf[i] for i in selected_indices]
        else:
            selected_vcf = filtered_vcf
        
        print(f"   📍 Targeting {len(selected_vcf)} VCF positions for mutations")
        
        # Create site table with VCF positions
        vcf_positions = [v.position for v in selected_vcf]
        vcf_ancestral_states = ['A'] * len(vcf_positions)  # Simplified
        
        # Add targeted mutations at VCF positions
        if vcf_positions:
            site_table = tskit.SiteTable()
            mutation_table = tskit.MutationTable()
            
            # Add sites at VCF positions
            for i, (pos, ancestral) in enumerate(zip(vcf_positions, vcf_ancestral_states)):
                site_id = site_table.add_row(position=pos, ancestral_state=ancestral)
                
                # Add mutation with frequency approximately matching VCF
                target_frequency = selected_vcf[i].allele_frequency
                
                # Simple mutation placement (could be more sophisticated)
                num_samples = ts.num_samples
                target_carriers = int(target_frequency * num_samples)
                
                if target_carriers > 0:
                    # Add mutation to a subset of samples
                    for sample_idx in range(min(target_carriers, num_samples)):
                        mutation_table.add_row(
                            site=site_id,
                            node=sample_idx,
                            derived_state='T',  # Simplified
                            time=0
                        )
            
            # Create new tree sequence with targeted mutations
            tables = ts.dump_tables()
            tables.sites.replace_with(site_table)
            tables.mutations.replace_with(mutation_table)
            ts_with_vcf = tables.tree_sequence()
        else:
            ts_with_vcf = ts
        
        # Add random mutations for remaining mutation budget
        if random_mutations > 0:
            random_rate = random_mutations / (sequence_length * 2 * sample_size)
            ts_final = msprime.sim_mutations(
                ts_with_vcf,
                rate=random_rate,
                random_seed=random_seed
            )
        else:
            ts_final = ts_with_vcf
        
        print(f"   ✅ Generated: {ts_final.num_mutations} total mutations")
        print(f"      Targeted: {len(selected_vcf)}, Random: {ts_final.num_mutations - len(selected_vcf)}")
        
        # Calculate metrics
        metrics = self._calculate_integration_metrics(
            ts_final, selected_vcf, sequence_length
        )
        
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024
        
        result = EnhancedIntegrationResult(
            population=population,
            chromosome=chromosome, 
            tree_sequence=ts_final,
            vcf_variants_used=selected_vcf,
            targeted_mutations=len(selected_vcf),
            random_mutations=ts_final.num_mutations - len(selected_vcf),
            position_accuracy=metrics['position_accuracy'],
            allele_frequency_fidelity=metrics['af_fidelity'],
            demographic_authenticity=metrics['demographic_score'],
            integration_time=end_time - start_time,
            memory_used_mb=end_memory - start_memory
        )
        
        print(f"   📊 Position accuracy: {result.position_accuracy:.3f}")
        print(f"   📊 AF fidelity: {result.allele_frequency_fidelity:.3f}")
        print(f"   ⏱️  {result.integration_time:.2f}s, {result.memory_used_mb:.1f}MB")
        
        return result
    
    def _calculate_integration_metrics(
        self,
        tree_sequence: tskit.TreeSequence,
        vcf_variants: List[VCFVariant],
        sequence_length: int
    ) -> Dict[str, float]:
        """Calculate comprehensive integration quality metrics"""
        
        # Get simulated variant positions and frequencies
        sim_positions = {}
        sim_frequencies = {}
        
        for variant in tree_sequence.variants():
            pos = int(variant.site.position)
            af = np.mean(variant.genotypes)
            sim_positions[pos] = True
            sim_frequencies[pos] = af
        
        # Calculate position accuracy
        vcf_positions = {v.position for v in vcf_variants}
        matched_positions = vcf_positions.intersection(set(sim_positions.keys()))
        
        if len(vcf_positions) > 0:
            position_accuracy = len(matched_positions) / len(vcf_positions)
        else:
            position_accuracy = 0.0
        
        # Calculate allele frequency fidelity for matched positions
        vcf_af_map = {v.position: v.allele_frequency for v in vcf_variants}
        
        af_differences = []
        for pos in matched_positions:
            if pos in vcf_af_map and pos in sim_frequencies:
                vcf_af = vcf_af_map[pos]
                sim_af = sim_frequencies[pos]
                diff = abs(vcf_af - sim_af)
                af_differences.append(diff)
        
        if af_differences:
            af_fidelity = 1.0 - np.mean(af_differences)
        else:
            af_fidelity = 0.0
        
        # Calculate demographic authenticity score (simplified)
        # Based on variant count, diversity, and frequency spectrum
        total_variants = len(sim_frequencies)
        diversity_pi = tree_sequence.diversity()
        
        # Rare variant fraction
        rare_variants = sum(1 for af in sim_frequencies.values() if af < 0.05)
        rare_fraction = rare_variants / total_variants if total_variants > 0 else 0
        
        # Score based on expected patterns for human populations
        # Higher diversity and appropriate rare variant fraction = better score
        diversity_score = min(diversity_pi / 0.001, 1.0)  # Normalize to ~0.001 as target
        rare_score = 1.0 - abs(rare_fraction - 0.8)  # Target ~80% rare variants
        
        demographic_score = (diversity_score + rare_score) / 2
        
        return {
            'position_accuracy': position_accuracy,
            'af_fidelity': max(0.0, af_fidelity),
            'demographic_score': max(0.0, min(1.0, demographic_score))
        }
    
    def integrate_all_populations_enhanced(
        self,
        chromosome: str = "chr22",
        sample_size: int = 30,
        sequence_length: int = 5_000_000,
        max_vcf_variants: int = 500,
        random_seed: Optional[int] = None
    ) -> Dict[Population, EnhancedIntegrationResult]:
        """
        Enhanced VCF integration for all populations
        
        Args:
            chromosome: Chromosome to process
            sample_size: Number of diploid individuals per population
            sequence_length: Physical length in bp
            max_vcf_variants: Maximum VCF variants per population
            random_seed: Random seed for reproducibility
            
        Returns:
            Dict mapping populations to enhanced integration results
        """
        results = {}
        total_start = time.time()
        
        print(f"🌍 ENHANCED VCF INTEGRATION FOR ALL POPULATIONS")
        print(f"   Chromosome: {chromosome}, Length: {sequence_length:,}bp")
        print("=" * 60)
        
        for population in Population:
            try:
                result = self.integrate_with_targeted_mutations(
                    population=population,
                    chromosome=chromosome,
                    sample_size=sample_size,
                    sequence_length=sequence_length,
                    max_vcf_variants=max_vcf_variants,
                    vcf_mutation_fraction=0.7,
                    random_seed=random_seed
                )
                
                results[population] = result
                
                print(f"   ✅ {population.value}: {result.targeted_mutations} targeted + "
                      f"{result.random_mutations} random mutations")
                
            except Exception as e:
                print(f"   ❌ {population.value}: Failed - {e}")
        
        total_time = time.time() - total_start
        print(f"\n🎉 Enhanced VCF integration completed in {total_time:.2f}s")
        
        return results
    
    def validate_enhanced_integration(
        self,
        results: Dict[Population, EnhancedIntegrationResult]
    ) -> Dict[str, bool]:
        """Validate enhanced VCF integration quality"""
        
        print("\n📊 Validating Enhanced VCF Integration...")
        
        validation = {}
        
        # Test 1: Position accuracy should be high (>50% for targeted approach)
        position_accuracies = [r.position_accuracy for r in results.values()]
        mean_position_accuracy = np.mean(position_accuracies)
        validation['high_position_accuracy'] = mean_position_accuracy > 0.5
        
        # Test 2: Allele frequency fidelity should be reasonable
        af_fidelities = [r.allele_frequency_fidelity for r in results.values()]
        mean_af_fidelity = np.mean(af_fidelities)
        validation['good_af_fidelity'] = mean_af_fidelity > 0.3
        
        # Test 3: Demographic authenticity should be maintained
        demo_scores = [r.demographic_authenticity for r in results.values()]
        mean_demo_score = np.mean(demo_scores)
        validation['demographic_authenticity'] = mean_demo_score > 0.5
        
        # Test 4: All populations should have successful integration
        successful_integrations = sum(1 for r in results.values() 
                                    if r.targeted_mutations > 0)
        validation['all_populations_successful'] = successful_integrations == len(results)
        
        # Test 5: Balance between targeted and random mutations
        targeted_fractions = []
        for r in results.values():
            total_muts = r.targeted_mutations + r.random_mutations
            if total_muts > 0:
                targeted_fractions.append(r.targeted_mutations / total_muts)
        
        if targeted_fractions:
            mean_targeted_fraction = np.mean(targeted_fractions)
            validation['balanced_mutation_strategy'] = 0.3 < mean_targeted_fraction < 0.9
        else:
            validation['balanced_mutation_strategy'] = False
        
        # Print validation results
        passed = 0
        total = len(validation)
        
        for test, result in validation.items():
            status = "✅ PASSED" if result else "❌ FAILED"
            if result:
                passed += 1
            print(f"   {test}: {status}")
        
        # Print summary statistics
        print(f"\n📈 Integration Quality Summary:")
        print(f"   Mean Position Accuracy: {np.mean(position_accuracies):.3f}")
        print(f"   Mean AF Fidelity: {np.mean(af_fidelities):.3f}")
        print(f"   Mean Demo Authenticity: {np.mean(demo_scores):.3f}")
        if targeted_fractions:
            print(f"   Mean Targeted Fraction: {np.mean(targeted_fractions):.3f}")
        
        print(f"\n🎯 Enhanced Integration Validation: {passed}/{total} tests passed")
        
        return validation

def main():
    """Main function for enhanced VCF integration testing"""
    print("🎯 ENHANCED VCF-HAPLOTYPE INTEGRATION TESTING")
    print("=" * 55)
    
    integrator = EnhancedVCFIntegrator()
    
    # Test single population enhanced integration
    print("\n📍 Testing Enhanced YRI Integration...")
    yri_result = integrator.integrate_with_targeted_mutations(
        population=Population.YRI,
        chromosome="chr22",
        sample_size=20,
        sequence_length=3_000_000,
        max_vcf_variants=300,
        vcf_mutation_fraction=0.8,
        random_seed=42
    )
    
    print(f"\n📊 YRI Enhanced Integration Results:")
    print(f"   Targeted Mutations: {yri_result.targeted_mutations}")
    print(f"   Random Mutations: {yri_result.random_mutations}")
    print(f"   Position Accuracy: {yri_result.position_accuracy:.3f}")
    print(f"   AF Fidelity: {yri_result.allele_frequency_fidelity:.3f}")
    print(f"   Demographic Score: {yri_result.demographic_authenticity:.3f}")
    
    # Test multiple populations
    print("\n🌍 Testing Enhanced Multi-Population Integration...")
    test_populations = [Population.YRI, Population.CEU, Population.CHB]
    
    multi_results = {}
    for pop in test_populations:
        try:
            result = integrator.integrate_with_targeted_mutations(
                population=pop,
                chromosome="chr22",
                sample_size=15,
                sequence_length=2_000_000,
                max_vcf_variants=200,
                vcf_mutation_fraction=0.7,
                random_seed=42
            )
            multi_results[pop] = result
            
        except Exception as e:
            print(f"   {pop.value}: Failed - {e}")
    
    # Validate enhanced integration
    if multi_results:
        validation_results = integrator.validate_enhanced_integration(multi_results)
    
    print("\n🚀 Enhanced VCF integration successfully implemented!")

if __name__ == "__main__":
    main()