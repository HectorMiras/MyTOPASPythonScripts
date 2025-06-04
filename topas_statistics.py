"""Module for handling TOPAS simulation statistics and measurements."""

from dataclasses import dataclass
from typing import Optional
import numpy as np


@dataclass
class TOPASStatistics:
    """Class to hold statistics from TOPAS output files.
    
    Attributes:
        sum_value (float): Sum of scored quantity
        mean (float): Mean value of scored quantity per history
        count_in_bin (int): Number of counts in the bin
        second_moment (float): Second moment of the distribution
        variance (float): Variance of the distribution per history
        standard_deviation (float): Standard deviation of the distribution per history
        histories_with_scorer_active (int): Number of histories that contributed to the scorer
        measurement_type (str): Type of measurement (e.g., 'DoseToMedium', 'EnergyDeposit')
        units (str): Units of the measurement
        run_variance (float): Variance between different runs (optional)
        n_runs (int): Number of runs combined (optional)
    """
    
    sum_value: float
    mean: float
    count_in_bin: int
    second_moment: float
    variance: float
    standard_deviation: float
    histories_with_scorer_active: int
    measurement_type: str
    units: str
    run_variance: float = 0.0
    n_runs: int = 1
    sum_squares: float = 0.0  # For calculating run variance
    
    @classmethod
    def from_dict(cls, stats: dict, measurement_type: str, units: str) -> 'TOPASStatistics':
        """Create a TOPASStatistics instance from a dictionary of statistics.
        
        Args:
            stats: Dictionary containing statistics from TOPAS output
            measurement_type: Type of measurement (e.g., 'DoseToMedium')
            units: Units of the measurement
        
        Returns:
            TOPASStatistics instance
        """
        return cls(
            sum_value=stats.get('Sum', 0.0),
            mean=stats.get('Mean', 0.0),
            count_in_bin=int(stats.get('Count_in_Bin', 0)),
            second_moment=stats.get('Second_Moment', 0.0),
            variance=stats.get('Variance', 0.0),
            standard_deviation=stats.get('Standard_Deviation', 0.0),
            histories_with_scorer_active=int(stats.get('Histories_with_Scorer_Active', 0)),
            measurement_type=measurement_type,
            units=units
        )
    
    def __str__(self) -> str:
        """Return a string representation of the statistics."""
        return (f"{self.measurement_type}: {self.sum_value:.6e} {self.units} "
                f"(mean: {self.mean:.6e} {self.units}/hist, "
                f"histories: {self.histories_with_scorer_active})")
    
    def __add__(self, other: 'TOPASStatistics') -> 'TOPASStatistics':
        """Add two TOPASStatistics instances.
        
        Used for combining results from multiple runs.
        
        Args:
            other: Another TOPASStatistics instance
            
        Returns:
            A new TOPASStatistics instance with combined values
            
        Raises:
            ValueError: If measurement types or units don't match
        """
        if self.measurement_type != other.measurement_type or self.units != other.units:
            raise ValueError("Cannot add statistics with different measurement types or units")
            
        total_histories = self.histories_with_scorer_active + other.histories_with_scorer_active
        if total_histories == 0:
            return self  # Return copy of self if no histories
            
        # Combine the statistics
        combined_sum = self.sum_value + other.sum_value
        combined_count = self.count_in_bin + other.count_in_bin
        combined_second_moment = self.second_moment + other.second_moment
        
        # Calculate new mean
        combined_mean = combined_sum / total_histories
        
        # Calculate new variance and standard deviation (per history)
        combined_variance = (combined_second_moment / total_histories - combined_mean ** 2)
        combined_std_dev = (combined_variance ** 0.5) if combined_variance > 0 else 0.0
        
        # Calculate variance between runs
        n_runs = self.n_runs + other.n_runs
        combined_sum_squares = self.sum_squares + other.sum_squares + (self.sum_value ** 2) + (other.sum_value ** 2)
        run_variance = np.sqrt(combined_sum_squares / n_runs - (combined_sum / n_runs) ** 2)
        
        return TOPASStatistics(
            sum_value=combined_sum,
            mean=combined_mean,
            count_in_bin=combined_count,
            second_moment=combined_second_moment,
            variance=combined_variance,
            standard_deviation=combined_std_dev,
            histories_with_scorer_active=total_histories,
            measurement_type=self.measurement_type,
            units=self.units,
            run_variance=run_variance,
            n_runs=n_runs,
            sum_squares=combined_sum_squares
        )
