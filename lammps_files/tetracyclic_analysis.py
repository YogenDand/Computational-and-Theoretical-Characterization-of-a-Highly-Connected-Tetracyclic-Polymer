#!/usr/bin/env python3
"""
Theoretical Graph Analysis for the Tetracyclic "Alpha" Polymer
Calculates the theoretical g-factor based on Cantarella et al. (2022)
"""

import numpy as np
import networkx as nx
from scipy.linalg import pinv
import matplotlib.pyplot as plt

class AlphaGraphAnalysis:
    """
    Performs theoretical analysis of the 6-vertex, 9-edge tetracyclic
    "Alpha" polymer graph from Figure 3 of Cantarella et al. (2022).
    """

    def __init__(self):
        """Initializes the analysis by creating the graph."""
        self.G = self.create_alpha_graph()
        self.expected_g_from_paper = 17/49

    def create_alpha_graph(self):
        """Creates the specific 6-vertex, 9-edge 'Alpha' graph."""
        G = nx.Graph()
        # Define the 6 junction vertices and 9 edges (chains)
        # as depicted in Figure 3 of the paper.
        nodes = [1, 2, 3, 4, 5, 6]
        edges = [
            (1, 2), (1, 4), (1, 5),
            (3, 2), (3, 4), (3, 6)
            (5, 2), (5, 4), (5, 6)
        ]
        G.add_nodes_from(nodes)
        G.add_edges_from(edges)
        return G

    def calculate_theoretical_g_factor(self):
        """
        Calculates the asymptotic g-factor using Theorem 5 from
        Cantarella et al. (2022).
        g = (3 / e²) * (Tr(L⁺) + Loops/3 - 1/6)
        """
        v = self.G.number_of_nodes()
        e = self.G.number_of_edges()
        
        # Cycle rank (Loops) = e - v + 1
        cycle_rank = e - v + 1

        # Get the normalized graph Laplacian
        L_norm = nx.normalized_laplacian_matrix(self.G).toarray()

        # Compute the Moore-Penrose pseudoinverse
        L_norm_plus = pinv(L_norm)

        # Get the trace of the pseudoinverse
        trace_L_plus = np.trace(L_norm_plus)
        
        # Apply the formula from Theorem 5
        g_factor = (3 / (e**2)) * (trace_L_plus + (cycle_rank / 3.0) - (1.0 / 6.0))
        
        print("--- Theoretical g-factor Calculation ---")
        print(f"Vertices (v): {v}")
        print(f"Edges (e):    {e}")
        print(f"Cycle Rank (Loops): {cycle_rank}")
        print(f"Trace of L_norm_plus: {trace_L_plus:.4f}")
        print(f"Calculated g-factor:  {g_factor:.6f}")
        print(f"Expected from paper:  {self.expected_g_from_paper:.6f} (109/245)")
        print("-" * 40)
        
        return g_factor

    def visualize_graph(self):
        """Visualize the graph structure."""
        plt.figure(figsize=(8, 6))
        pos = nx.spring_layout(self.G, seed=42)
        nx.draw(self.G, pos, with_labels=True, node_color='skyblue',
                node_size=1200, font_size=14, font_weight='bold',
                width=2.0, edge_color='gray')
        plt.title(f'Tetracyclic "Alpha" Polymer Graph (g ≈ 0.445)')
        plt.savefig('alpha_graph_structure.png', dpi=300)
        print("Graph visualization saved as alpha_graph_structure.png")

def main():
    """Run the full theoretical analysis."""
    print("=" * 60)
    print("Theoretical Analysis of the 'Alpha' Polymer Graph")
    print("=" * 60)
    
    analyzer = AlphaGraphAnalysis()
    analyzer.visualize_graph()
    theoretical_g = analyzer.calculate_theoretical_g_factor()

    print("\n✅ Analysis complete.")
    print("This theoretical g-factor is the benchmark for your LAMMPS simulation.")

if __name__ == "__main__":
    main()
