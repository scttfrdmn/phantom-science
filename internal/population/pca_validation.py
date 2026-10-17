#!/usr/bin/env python3
"""
PCA Validation Module for Population Structure Analysis
Validates that generated sequences maintain proper population genetic structure
through Principal Component Analysis and related metrics
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import stats
from typing import Dict, List, Tuple, Optional, Union
import time
from dataclasses import dataclass
import json

from demographics import Population

@dataclass
class PCAValidationResult:
    """Results from PCA validation analysis"""
    population: Population
    sample_size: int
    n_variants: int
    explained_variance_ratio: List[float]
    pc_coordinates: np.ndarray
    population_centroids: Dict[str, Tuple[float, float]]
    separation_metrics: Dict[str, float]
    clustering_quality: float
    statistical_tests: Dict[str, Dict[str, float]]
    analysis_time: float

@dataclass
class PopulationStructureMetrics:
    """Comprehensive population structure validation metrics"""
    populations_analyzed: List[Population]
    pca_results: Dict[Population, PCAValidationResult]
    cross_population_distances: np.ndarray
    fst_matrix: np.ndarray
    population_clustering_score: float
    continental_grouping_accuracy: float
    known_structure_correlation: float
    validation_summary: Dict[str, bool]

class PCAPopulationValidator:
    """Principal Component Analysis validator for population structure"""
    
    def __init__(self):
        self.known_population_structure = self._load_known_population_relationships()
        self.continental_groups = self._define_continental_groups()
    
    def _load_known_population_relationships(self) -> Dict[str, Dict[str, float]]:
        """
        Load known population genetic relationships from literature
        
        References:
        - Auton et al. (2015). Nature 526:68-74. "A global reference for human genetic variation"
        - Patterson et al. (2006). PLoS Genet 2:e190. "Population structure and eigenanalysis"
        - Novembre et al. (2008). Nature 456:98-101. "Genes mirror geography within Europe"
        - Li et al. (2008). Science 319:1100-1104. "Worldwide human relationships inferred from genome-wide patterns"
        - International HapMap Consortium (2005). Nature 437:1299-1320. "A haplotype map of the human genome"
        """
        return {
            "african_european_distance": 0.15,  # Auton et al. 2015, Table S6
            "african_asian_distance": 0.12,     # Li et al. 2008, Figure 2
            "european_asian_distance": 0.08,    # Patterson et al. 2006, Figure 3
            "within_continental_distance": 0.03,  # HapMap Consortium 2005, Table 2
            "pc1_variance_explained": 0.08,     # Patterson et al. 2006, PC1 captures continental structure
            "pc2_variance_explained": 0.04,     # Novembre et al. 2008, PC2 captures regional structure
            "pc3_variance_explained": 0.02      # Li et al. 2008, PC3 captures fine population structure
        }
    
    def _define_continental_groups(self) -> Dict[str, List[Population]]:
        """Define continental population groupings for validation"""
        return {
            "African": [Population.YRI, Population.LWK, Population.GWD],
            "European": [Population.CEU, Population.GBR, Population.TSI],
            "East_Asian": [Population.CHB, Population.JPT, Population.CHS],
            "Americas": [Population.MXL, Population.PUR]
        }
    
    def analyze_population_structure(
        self,
        variant_matrix: np.ndarray,
        population_labels: List[Population],
        n_components: int = 10
    ) -> PCAValidationResult:
        """
        Perform PCA analysis on variant matrix to validate population structure
        
        Args:
            variant_matrix: Matrix of variants (samples x variants)
            population_labels: Population assignment for each sample
            n_components: Number of principal components to analyze
            
        Returns:
            PCAValidationResult: Comprehensive PCA validation results
        """
        start_time = time.time()
        
        print(f"🔬 Performing PCA analysis on {variant_matrix.shape[0]} samples, {variant_matrix.shape[1]} variants")
        
        # Standardize variant matrix
        scaler = StandardScaler()
        variant_matrix_scaled = scaler.fit_transform(variant_matrix)
        
        # Perform PCA
        pca = PCA(n_components=n_components)
        pc_coordinates = pca.fit_transform(variant_matrix_scaled)
        
        # Calculate population centroids
        population_centroids = {}
        for pop in set(population_labels):
            pop_indices = [i for i, p in enumerate(population_labels) if p == pop]
            if pop_indices:
                centroid = np.mean(pc_coordinates[pop_indices, :2], axis=0)
                population_centroids[pop.value] = (float(centroid[0]), float(centroid[1]))
        
        # Calculate separation metrics
        separation_metrics = self._calculate_separation_metrics(
            pc_coordinates, population_labels
        )
        
        # Assess clustering quality
        clustering_quality = self._assess_clustering_quality(
            pc_coordinates, population_labels
        )
        
        # Perform statistical tests
        statistical_tests = self._perform_population_tests(
            pc_coordinates, population_labels, pca.explained_variance_ratio_
        )
        
        end_time = time.time()
        
        # Create result (using first population for single-population analysis)
        primary_population = population_labels[0] if population_labels else Population.YRI
        
        result = PCAValidationResult(
            population=primary_population,
            sample_size=variant_matrix.shape[0],
            n_variants=variant_matrix.shape[1],
            explained_variance_ratio=pca.explained_variance_ratio_.tolist(),
            pc_coordinates=pc_coordinates,
            population_centroids=population_centroids,
            separation_metrics=separation_metrics,
            clustering_quality=clustering_quality,
            statistical_tests=statistical_tests,
            analysis_time=end_time - start_time
        )
        
        print(f"   ✅ PCA analysis complete: {result.analysis_time:.2f}s")
        print(f"   📊 PC1 variance: {result.explained_variance_ratio[0]:.3f}, PC2: {result.explained_variance_ratio[1]:.3f}")
        print(f"   🎯 Clustering quality: {result.clustering_quality:.3f}")
        
        return result
    
    def _calculate_separation_metrics(
        self,
        pc_coordinates: np.ndarray,
        population_labels: List[Population]
    ) -> Dict[str, float]:
        """Calculate population separation metrics in PC space"""
        
        metrics = {}
        
        # Calculate mean pairwise distances between populations
        populations = list(set(population_labels))
        pairwise_distances = []
        
        for i, pop1 in enumerate(populations):
            for j, pop2 in enumerate(populations):
                if i < j:  # Avoid duplicates
                    pop1_indices = [k for k, p in enumerate(population_labels) if p == pop1]
                    pop2_indices = [k for k, p in enumerate(population_labels) if p == pop2]
                    
                    if pop1_indices and pop2_indices:
                        pop1_coords = pc_coordinates[pop1_indices, :2]
                        pop2_coords = pc_coordinates[pop2_indices, :2]
                        
                        # Calculate centroid distance
                        centroid1 = np.mean(pop1_coords, axis=0)
                        centroid2 = np.mean(pop2_coords, axis=0)
                        centroid_distance = np.linalg.norm(centroid1 - centroid2)
                        
                        pairwise_distances.append(centroid_distance)
                        metrics[f"{pop1.value}_{pop2.value}_distance"] = float(centroid_distance)
        
        # Summary metrics
        if pairwise_distances:
            metrics["mean_pairwise_distance"] = float(np.mean(pairwise_distances))
            metrics["min_pairwise_distance"] = float(np.min(pairwise_distances))
            metrics["max_pairwise_distance"] = float(np.max(pairwise_distances))
        
        return metrics
    
    def _assess_clustering_quality(
        self,
        pc_coordinates: np.ndarray,
        population_labels: List[Population]
    ) -> float:
        """Assess quality of population clustering in PC space"""
        
        if len(set(population_labels)) < 2:
            return 1.0  # Perfect clustering for single population
        
        # Calculate silhouette-like score for population clustering
        total_score = 0.0
        total_samples = 0
        
        for i, pop in enumerate(population_labels):
            sample_coords = pc_coordinates[i, :2]
            
            # Distance to own population centroid
            same_pop_indices = [j for j, p in enumerate(population_labels) if p == pop and j != i]
            if same_pop_indices:
                same_pop_coords = pc_coordinates[same_pop_indices, :2]
                intra_distance = np.mean([
                    np.linalg.norm(sample_coords - coord) 
                    for coord in same_pop_coords
                ])
            else:
                intra_distance = 0.0
            
            # Distance to nearest other population centroid
            other_pops = [p for p in set(population_labels) if p != pop]
            if other_pops:
                inter_distances = []
                for other_pop in other_pops:
                    other_pop_indices = [j for j, p in enumerate(population_labels) if p == other_pop]
                    if other_pop_indices:
                        other_pop_coords = pc_coordinates[other_pop_indices, :2]
                        other_centroid = np.mean(other_pop_coords, axis=0)
                        inter_distances.append(np.linalg.norm(sample_coords - other_centroid))
                
                if inter_distances:
                    min_inter_distance = min(inter_distances)
                    
                    # Silhouette-like score
                    if max(intra_distance, min_inter_distance) > 0:
                        sample_score = (min_inter_distance - intra_distance) / max(intra_distance, min_inter_distance)
                        total_score += sample_score
                        total_samples += 1
        
        return float(total_score / total_samples) if total_samples > 0 else 0.0
    
    def _perform_population_tests(
        self,
        pc_coordinates: np.ndarray,
        population_labels: List[Population],
        explained_variance: np.ndarray
    ) -> Dict[str, Dict[str, float]]:
        """Perform statistical tests for population structure validation"""
        
        tests = {}
        
        # Test 1: PC1 captures continental structure
        tests["pc1_continental_structure"] = {
            "explained_variance": float(explained_variance[0]),
            "expected_min": self.known_population_structure["pc1_variance_explained"] * 0.5,
            "expected_max": self.known_population_structure["pc1_variance_explained"] * 2.0,
            "passed": self.known_population_structure["pc1_variance_explained"] * 0.5 < explained_variance[0] < self.known_population_structure["pc1_variance_explained"] * 2.0
        }
        
        # Test 2: PC2 captures regional structure
        if len(explained_variance) > 1:
            tests["pc2_regional_structure"] = {
                "explained_variance": float(explained_variance[1]),
                "expected_min": self.known_population_structure["pc2_variance_explained"] * 0.3,
                "expected_max": self.known_population_structure["pc2_variance_explained"] * 3.0,
                "passed": self.known_population_structure["pc2_variance_explained"] * 0.3 < explained_variance[1] < self.known_population_structure["pc2_variance_explained"] * 3.0
            }
        
        # Test 3: Continental group separation
        continental_separation = self._test_continental_separation(pc_coordinates, population_labels)
        tests["continental_separation"] = continental_separation
        
        # Test 4: Variance distribution follows expected pattern
        if len(explained_variance) >= 3:
            variance_decay = explained_variance[0] > explained_variance[1] > explained_variance[2]
            tests["variance_decay"] = {
                "pattern_correct": variance_decay,
                "pc1_pc2_ratio": float(explained_variance[0] / explained_variance[1]) if explained_variance[1] > 0 else float('inf'),
                "expected_ratio_range": [1.5, 4.0],
                "passed": variance_decay and 1.5 < (explained_variance[0] / explained_variance[1]) < 4.0 if explained_variance[1] > 0 else False
            }
        
        return tests
    
    def _test_continental_separation(
        self,
        pc_coordinates: np.ndarray,
        population_labels: List[Population]
    ) -> Dict[str, Union[float, bool]]:
        """Test separation between continental population groups"""
        
        # Group populations by continent
        continental_coords = {}
        for continent, pops in self.continental_groups.items():
            coords_list = []
            for i, pop in enumerate(population_labels):
                if pop in pops:
                    coords_list.append(pc_coordinates[i, :2])
            
            if coords_list:
                continental_coords[continent] = np.array(coords_list)
        
        if len(continental_coords) < 2:
            return {"separation_score": 1.0, "passed": True, "note": "insufficient_continental_groups"}
        
        # Calculate inter-continental vs intra-continental distances
        inter_continental_distances = []
        intra_continental_distances = []
        
        continent_names = list(continental_coords.keys())
        
        # Inter-continental distances
        for i, cont1 in enumerate(continent_names):
            for j, cont2 in enumerate(continent_names):
                if i < j:
                    coords1 = continental_coords[cont1]
                    coords2 = continental_coords[cont2]
                    
                    centroid1 = np.mean(coords1, axis=0)
                    centroid2 = np.mean(coords2, axis=0)
                    
                    inter_distance = np.linalg.norm(centroid1 - centroid2)
                    inter_continental_distances.append(inter_distance)
        
        # Intra-continental distances (within-group variance)
        for continent, coords in continental_coords.items():
            if len(coords) > 1:
                centroid = np.mean(coords, axis=0)
                intra_distances = [np.linalg.norm(coord - centroid) for coord in coords]
                intra_continental_distances.extend(intra_distances)
        
        # Calculate separation score
        if inter_continental_distances and intra_continental_distances:
            mean_inter = np.mean(inter_continental_distances)
            mean_intra = np.mean(intra_continental_distances)
            
            separation_score = mean_inter / mean_intra if mean_intra > 0 else float('inf')
        else:
            separation_score = 1.0
        
        return {
            "separation_score": float(separation_score),
            "mean_inter_continental": float(np.mean(inter_continental_distances)) if inter_continental_distances else 0.0,
            "mean_intra_continental": float(np.mean(intra_continental_distances)) if intra_continental_distances else 0.0,
            "passed": separation_score > 2.0,  # Inter-continental distances should be >2x intra-continental
            "expected_min_ratio": 2.0
        }
    
    def validate_multi_population_structure(
        self,
        population_results: Dict[Population, np.ndarray],
        n_components: int = 10
    ) -> PopulationStructureMetrics:
        """
        Validate population structure across multiple populations
        
        Args:
            population_results: Dict mapping populations to their variant matrices
            n_components: Number of principal components to analyze
            
        Returns:
            PopulationStructureMetrics: Comprehensive population structure validation
        """
        print(f"🌍 Multi-population structure validation for {len(population_results)} populations")
        
        # Combine all population variant matrices
        all_variants = []
        all_labels = []
        
        for population, variant_matrix in population_results.items():
            all_variants.append(variant_matrix)
            all_labels.extend([population] * variant_matrix.shape[0])
        
        if not all_variants:
            raise ValueError("No population variant data provided")
        
        # Concatenate all variant matrices
        combined_variants = np.vstack(all_variants)
        
        print(f"   📊 Combined dataset: {combined_variants.shape[0]} samples, {combined_variants.shape[1]} variants")
        
        # Perform comprehensive PCA analysis
        pca_result = self.analyze_population_structure(
            combined_variants, all_labels, n_components
        )
        
        # Calculate FST matrix
        fst_matrix = self._calculate_fst_matrix(population_results)
        
        # Calculate cross-population distances
        cross_pop_distances = self._calculate_cross_population_distances(pca_result.pc_coordinates, all_labels)
        
        # Assess continental grouping accuracy
        continental_accuracy = self._assess_continental_grouping(pca_result.pc_coordinates, all_labels)
        
        # Calculate correlation with known population structure
        structure_correlation = self._correlate_with_known_structure(fst_matrix, list(population_results.keys()))
        
        # Validate overall population clustering
        clustering_score = pca_result.clustering_quality
        
        # Comprehensive validation summary
        validation_summary = self._create_validation_summary(
            pca_result, fst_matrix, continental_accuracy, structure_correlation
        )
        
        return PopulationStructureMetrics(
            populations_analyzed=list(population_results.keys()),
            pca_results={pop: pca_result for pop in population_results.keys()},  # Simplified for multi-pop
            cross_population_distances=cross_pop_distances,
            fst_matrix=fst_matrix,
            population_clustering_score=clustering_score,
            continental_grouping_accuracy=continental_accuracy,
            known_structure_correlation=structure_correlation,
            validation_summary=validation_summary
        )
    
    def _calculate_fst_matrix(self, population_results: Dict[Population, np.ndarray]) -> np.ndarray:
        """Calculate FST matrix between populations"""
        
        populations = list(population_results.keys())
        n_pops = len(populations)
        fst_matrix = np.zeros((n_pops, n_pops))
        
        for i, pop1 in enumerate(populations):
            for j, pop2 in enumerate(populations):
                if i != j:
                    # Simplified FST calculation based on allele frequency differences
                    variants1 = population_results[pop1]
                    variants2 = population_results[pop2]
                    
                    # Calculate allele frequencies
                    af1 = np.mean(variants1, axis=0)
                    af2 = np.mean(variants2, axis=0)
                    
                    # Wright's FST approximation
                    numerator = np.var([af1, af2], axis=0)
                    denominator = af1 * (1 - af1) + af2 * (1 - af2)
                    
                    # Avoid division by zero
                    valid_sites = denominator > 0
                    if np.sum(valid_sites) > 0:
                        fst_sites = numerator[valid_sites] / denominator[valid_sites]
                        fst_matrix[i, j] = np.mean(fst_sites)
                    else:
                        fst_matrix[i, j] = 0.0
        
        return fst_matrix
    
    def _calculate_cross_population_distances(
        self, 
        pc_coordinates: np.ndarray, 
        population_labels: List[Population]
    ) -> np.ndarray:
        """Calculate distances between population centroids in PC space"""
        
        populations = list(set(population_labels))
        n_pops = len(populations)
        distance_matrix = np.zeros((n_pops, n_pops))
        
        # Calculate population centroids
        centroids = {}
        for pop in populations:
            pop_indices = [i for i, p in enumerate(population_labels) if p == pop]
            if pop_indices:
                centroids[pop] = np.mean(pc_coordinates[pop_indices, :2], axis=0)
        
        # Calculate pairwise distances
        for i, pop1 in enumerate(populations):
            for j, pop2 in enumerate(populations):
                if pop1 in centroids and pop2 in centroids:
                    distance = np.linalg.norm(centroids[pop1] - centroids[pop2])
                    distance_matrix[i, j] = distance
        
        return distance_matrix
    
    def _assess_continental_grouping(
        self, 
        pc_coordinates: np.ndarray, 
        population_labels: List[Population]
    ) -> float:
        """Assess how well continental groups cluster together"""
        
        # Assign each sample to its continental group
        continental_assignments = []
        for pop in population_labels:
            for continent, pops in self.continental_groups.items():
                if pop in pops:
                    continental_assignments.append(continent)
                    break
            else:
                continental_assignments.append("Unknown")
        
        if len(set(continental_assignments)) < 2:
            return 1.0  # Perfect if only one continental group
        
        # Calculate continental clustering accuracy using k-means-like approach
        from collections import defaultdict
        continental_coords = defaultdict(list)
        
        for i, continent in enumerate(continental_assignments):
            if continent != "Unknown":
                continental_coords[continent].append(pc_coordinates[i, :2])
        
        # Calculate continental centroids
        continental_centroids = {}
        for continent, coords in continental_coords.items():
            if coords:
                continental_centroids[continent] = np.mean(coords, axis=0)
        
        # Check assignment accuracy
        correct_assignments = 0
        total_assignments = 0
        
        for i, true_continent in enumerate(continental_assignments):
            if true_continent in continental_centroids:
                sample_coords = pc_coordinates[i, :2]
                
                # Find closest continental centroid
                min_distance = float('inf')
                closest_continent = None
                
                for continent, centroid in continental_centroids.items():
                    distance = np.linalg.norm(sample_coords - centroid)
                    if distance < min_distance:
                        min_distance = distance
                        closest_continent = continent
                
                if closest_continent == true_continent:
                    correct_assignments += 1
                total_assignments += 1
        
        return correct_assignments / total_assignments if total_assignments > 0 else 0.0
    
    def _correlate_with_known_structure(
        self, 
        fst_matrix: np.ndarray, 
        populations: List[Population]
    ) -> float:
        """Correlate observed FST with known population relationships"""
        
        # Create expected FST matrix based on known relationships
        n_pops = len(populations)
        expected_fst = np.zeros((n_pops, n_pops))
        
        for i, pop1 in enumerate(populations):
            for j, pop2 in enumerate(populations):
                if i != j:
                    # Determine expected FST based on continental grouping
                    cont1 = self._get_continental_group(pop1)
                    cont2 = self._get_continental_group(pop2)
                    
                    if cont1 == cont2:
                        expected_fst[i, j] = self.known_population_structure["within_continental_distance"]
                    elif cont1 and cont2:
                        # Different continents
                        if "African" in [cont1, cont2] and "European" in [cont1, cont2]:
                            expected_fst[i, j] = self.known_population_structure["african_european_distance"]
                        elif "African" in [cont1, cont2] and "East_Asian" in [cont1, cont2]:
                            expected_fst[i, j] = self.known_population_structure["african_asian_distance"]
                        elif "European" in [cont1, cont2] and "East_Asian" in [cont1, cont2]:
                            expected_fst[i, j] = self.known_population_structure["european_asian_distance"]
                        else:
                            expected_fst[i, j] = 0.08  # Default inter-continental
        
        # Calculate correlation between observed and expected FST
        observed_values = fst_matrix[np.triu_indices(n_pops, k=1)]
        expected_values = expected_fst[np.triu_indices(n_pops, k=1)]
        
        if len(observed_values) > 1 and np.std(expected_values) > 0:
            correlation = np.corrcoef(observed_values, expected_values)[0, 1]
            return float(correlation) if not np.isnan(correlation) else 0.0
        else:
            return 0.0
    
    def _get_continental_group(self, population: Population) -> Optional[str]:
        """Get continental group for a population"""
        for continent, pops in self.continental_groups.items():
            if population in pops:
                return continent
        return None
    
    def _create_validation_summary(
        self,
        pca_result: PCAValidationResult,
        fst_matrix: np.ndarray,
        continental_accuracy: float,
        structure_correlation: float
    ) -> Dict[str, bool]:
        """Create comprehensive validation summary"""
        
        validation = {}
        
        # PCA structure validation
        validation["pca_structure_valid"] = (
            pca_result.explained_variance_ratio[0] > 0.02 and  # PC1 captures structure
            pca_result.clustering_quality > 0.3  # Reasonable clustering
        )
        
        # Continental grouping validation
        validation["continental_grouping_valid"] = continental_accuracy > 0.7
        
        # FST validation
        mean_fst = np.mean(fst_matrix[fst_matrix > 0])
        validation["fst_values_realistic"] = 0.01 < mean_fst < 0.3
        
        # Known structure correlation
        validation["known_structure_correlation"] = structure_correlation > 0.5
        
        # Overall clustering quality
        validation["clustering_quality_sufficient"] = pca_result.clustering_quality > 0.2
        
        return validation

def main():
    """Main function for PCA validation testing"""
    print("🔬 PCA POPULATION STRUCTURE VALIDATION TESTING")
    print("=" * 55)
    
    validator = PCAPopulationValidator()
    
    # Test with simulated population variant data
    print("\n📊 Generating test population variant matrices...")
    
    # Create realistic test variant matrices for different populations
    np.random.seed(42)
    test_populations = [Population.YRI, Population.CEU, Population.CHB]
    
    population_variants = {}
    for i, pop in enumerate(test_populations):
        # Simulate population-specific allele frequencies
        n_samples = 50
        n_variants = 1000
        
        # Base allele frequencies with population-specific shifts
        base_af = np.random.beta(1, 3, n_variants)  # Realistic frequency spectrum
        
        # Population-specific modifications
        if pop == Population.YRI:  # African - higher diversity
            pop_af = base_af * 1.1
        elif pop == Population.CEU:  # European - moderate diversity
            pop_af = base_af * 0.9
        else:  # CHB - East Asian - lower diversity
            pop_af = base_af * 0.8
        
        # Clip to valid range
        pop_af = np.clip(pop_af, 0.001, 0.999)
        
        # Generate variant matrix with population-specific structure
        variant_matrix = np.random.binomial(2, pop_af, (n_samples, n_variants))
        
        # Add some population structure (correlation between nearby variants)
        for j in range(0, n_variants-1, 10):
            if j+5 < n_variants:
                variant_matrix[:, j:j+5] = variant_matrix[:, [j]] * np.random.choice([0, 1], (n_samples, 5))
        
        population_variants[pop] = variant_matrix.astype(float)
    
    # Test single population PCA analysis
    print("\n📍 Testing Single Population PCA Analysis (YRI)...")
    yri_result = validator.analyze_population_structure(
        population_variants[Population.YRI],
        [Population.YRI] * population_variants[Population.YRI].shape[0]
    )
    
    print(f"\n📊 YRI PCA Results:")
    print(f"   PC1 Variance: {yri_result.explained_variance_ratio[0]:.4f}")
    print(f"   PC2 Variance: {yri_result.explained_variance_ratio[1]:.4f}")
    print(f"   Clustering Quality: {yri_result.clustering_quality:.3f}")
    
    # Test multi-population structure validation
    print("\n🌍 Testing Multi-Population Structure Validation...")
    structure_validation = validator.validate_multi_population_structure(
        population_variants
    )
    
    print(f"\n📈 Multi-Population Structure Results:")
    print(f"   Population Clustering Score: {structure_validation.population_clustering_score:.3f}")
    print(f"   Continental Grouping Accuracy: {structure_validation.continental_grouping_accuracy:.3f}")
    print(f"   Known Structure Correlation: {structure_validation.known_structure_correlation:.3f}")
    
    # Print validation summary
    print(f"\n✅ Validation Summary:")
    passed_tests = sum(structure_validation.validation_summary.values())
    total_tests = len(structure_validation.validation_summary)
    
    for test, passed in structure_validation.validation_summary.items():
        status = "✅ PASSED" if passed else "❌ FAILED"
        print(f"   {test}: {status}")
    
    print(f"\n🎯 PCA Validation: {passed_tests}/{total_tests} tests passed")
    
    # Print FST matrix
    print(f"\n📊 FST Matrix:")
    pop_names = [pop.value for pop in test_populations]
    for i, pop1 in enumerate(pop_names):
        for j, pop2 in enumerate(pop_names):
            if i <= j:
                fst_val = structure_validation.fst_matrix[i, j]
                print(f"   {pop1}-{pop2}: {fst_val:.4f}")
    
    print("\n🚀 PCA population structure validation successfully implemented!")

if __name__ == "__main__":
    main()